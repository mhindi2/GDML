# -*- coding: utf-8 -*-
# Fri Feb 11 01:11:14 PM PST 2022
# **************************************************************************
# *                                                                        *
# *   Copyright (c) 2021 Keith Sloan <keith@sloan-home.co.uk>              *
# *                                                                        *
# *   This program is free software; you can redistribute it and/or modify *
# *   it under the terms of the GNU Lesser General Public License (LGPL)   *
# *   as published by the Free Software Foundation; either version 2 of    *
# *   the License, or (at your option) any later version.                  *
# *   for detail see the LICENCE text file.                                *
# *                                                                        *
# *   This program is distributed in the hope that it will be useful,      *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of       *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
# *   GNU Library General Public License for more details.                 *
# *                                                                        *
# *   You should have received a copy of the GNU Library General Public    *
# *   License along with this program; if not, write to the Free Software  *
# *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 *
# *   USA                                                                  *
# *                                                                        *
# *   Acknowledgements :
# *                                                                        *
# *                                                                        *
# **************************************************************************
__title__ = "FreeCAD - OBJ -> GDML Tessellated importer"
__author__ = "Keith Sloan <keith@sloan-home.co.uk>"
__url__ = ["https://github.com/KeithSloan/FreeCAD_GDML"]

import FreeCAD, FreeCADGui
import os, io, sys, re
import Part, Draft

from PySide import QtGui, QtCore

def joinDir(path):
    import os
    __dirname__ = os.path.dirname(__file__)
    return(os.path.join(__dirname__, path))

# Save the native open function to avoid collisions
if open.__module__ in ['__builtin__', 'io']:
    pythonopen = open  # to distinguish python built-in open function from the one declared her


def open(filename):
    "called when freecad opens a file."

    print('Open : '+filename)
    docName = os.path.splitext(os.path.basename(filename))[0]
    print('path : '+filename)
    if filename.lower().endswith('.obj'):
        try:
            doc = FreeCAD.ActiveDocument()
            print('Active Doc')

        except:
            print('New Doc')
            doc = FreeCAD.newDocument(docName)

        processOBJ(doc, filename)
        return doc


def insert(filename, docname):
    "called when freecad imports a file"
    print('Insert filename : '+filename+' docname : '+docname)
    global doc
    groupname = os.path.splitext(os.path.basename(filename))[0]
    try:
        doc = FreeCAD.getDocument(docname)
    except NameError:
        doc = FreeCAD.newDocument(docname)
    if filename.lower().endswith('.obj'):
        processOBJ(doc, filename)


class switch(object):
    value = None

    def __new__(class_, value):
        class_.value = value
        return True


def case(*args):
    return any((arg == switch.value for arg in args))


class ColourWidget(QtGui.QLineEdit):

    def __init__(self, colour):
        super().__init__()
        palette = self.palette()
        # palette.setColor(QtGui.QPalette.Button, QtGui.QColor(colour))
        palette.setColor(QtGui.QPalette.Base, QtGui.QColor(colour))
        # palette.setColor(QtGui.QPalette.Window, QtGui.QColor(colour))
        self.setAutoFillBackground(True)
        # palette.setColor(QtGui.QPalette.Window, QtGui.QColor(QtGui.qRed))
        # palette.setColor(QtGui.QPalette.Window, QtGui.qRed)
        self.setStyleSheet("QPushButton {border-color: black; border: 2px;}")
        self.setPalette(palette)
        self.setReadOnly(True)
        self.setMaxLength = 5
        # self.setFlat(True)
        self.update()

class TextWidget(QtGui.QLineEdit):

    def __init__(self, text):
        super().__init__()
        #print(f"Text Widget {text}")
        self.insert(text)
        self.setReadOnly(True)

class MapMaterialObj2GDML(QtGui.QWidget):

    def __init__(self, objName, objMat, gdmlWidget, colour):
        super().__init__()
        #print('Map Material Entry Obj ->  GDML: ')
        #print(f'Name {objName} Material {objMat}')
        self.hbox = QtGui.QHBoxLayout()
        self.hbox.addWidget(TextWidget(objName))
        self.hbox.addWidget(TextWidget(objMat))
        self.hbox.addWidget(gdmlWidget)
        self.hbox.addWidget(ColourWidget(colour))
        self.setLayout(self.hbox)


class MaterialMapList(QtGui.QScrollArea):

    def __init__(self):
        super().__init__()
        from .GDMLMaterials import getMaterialsList
        # Scroll Area which contains the widgets, set as the centralWidget
        # Widget that contains the collection of Vertical Box
        self.widget = QtGui.QWidget()
        self.matList = getMaterialsList()
        #print(f"Material List {self.matList}")
        # The Vertical Box that contains the Horizontal Boxes of  labels and buttons
        self.vbox = QtGui.QVBoxLayout()
        self.widget.setLayout(self.vbox)

        # Scroll Area Properties
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.setWidgetResizable(True)
        self.setWidget(self.widget)

    def addEntry(self, objName, objMat, gdmlMat, colour):
        from .GDMLMaterials import GDMLMaterial
        #print('Add Entry')
        matWidget = GDMLMaterial(self.matList, gdmlMat)
        self.vbox.addWidget(MapMaterialObj2GDML(objName, objMat,  matWidget, colour))

class MapObjmat2GDMLmatDialog(QtGui.QDialog):
    def __init__(self, *args):
        super(MapObjmat2GDMLmatDialog, self).__init__()
        self.setupUi()
        self.initUI()

    def initUI(self):
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.setMouseTracking(True)
        self.show()

    def setupUi(self):
        self.setObjectName("Dialog")
        self.resize(400, 362)
        mainLayout = QtGui.QVBoxLayout()
        self.buttonBox = QtGui.QDialogButtonBox(self)
        self.buttonBox.setGeometry(QtCore.QRect(30, 30, 541, 130))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(
            QtGui.QDialogButtonBox.Cancel | QtGui.QDialogButtonBox.Ok
        )
        self.buttonBox.setObjectName("buttonBox")
        self.buttonBox.accepted.connect(self.action)
        self.buttonBox.rejected.connect(self.onCancel)

        #self.textEdit = QtGui.QTextEdit(self)
        #self.textEdit.setGeometry(QtCore.QRect(10, 10, 381, 141))
        #self.textEdit.setLocale(
        #    QtCore.QLocale(
        #       QtCore.QLocale.English, QtCore.QLocale.UnitedKingdom
        #    )
        #)
        #self.textEdit.setReadOnly(True)
        #self.textEdit.setObjectName("textEdit")
        #self.verticalLayoutWidget = QtGui.QWidget(self)
        #self.verticalLayoutWidget.setGeometry(QtCore.QRect(60, 150, 271, 151))
        #self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        #self.verticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        #self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        #self.verticalLayout.setObjectName("verticalLayout")
        #self.setLayout(self.verticalLayout)
        #self.groupBox = QtGui.QGroupBox(self.verticalLayoutWidget)
        #self.groupBox.setObjectName("groupBox")
        #self.fullDisplayRadioButton = QtGui.QRadioButton(self.groupBox)
        #self.fullDisplayRadioButton.setGeometry(QtCore.QRect(10, 30, 105, 22))
        #self.fullDisplayRadioButton.setChecked(True)
        #self.fullDisplayRadioButton.setObjectName("fullDisplayRadioButton")
        #self.samplesRadioButton = QtGui.QRadioButton(self.groupBox)
        self.mapList = MaterialMapList()
        self.mapLayout = QtGui.QVBoxLayout()
        self.mapLayout.addWidget(self.mapList)
        mainLayout.addLayout(self.mapLayout)
        mainLayout.addWidget(self.buttonBox)
        self.setLayout(mainLayout)
        self.setGeometry(30, 30, 800, 250)
        self.setWindowTitle("Map Materials Obj -> GDML")
        #self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.retStatus = 0

    def initMaterials(self):
        from .GDMLMaterials import getMaterialsList, GDMLMaterial
        self.materialsList = getMaterialsList()

    def addMaterialMapping(self, name, objMat, gdmlMat):
        #print(f"Add Material Map Obj {name} Material {objMat} to GDML mat {gdmlMat}")        
        self.mapList.addEntry(name, objMat, gdmlMat, None)

    def fullDisplayRadioButtonToggled(self):
        self.fullDisplayRadioButton.blockSignals(True)

        if self.fullDisplayRadioButton.isChecked():
            self.fractionSpinBox.setEnabled(False)
            self.fractionsLabel.setEnabled(False)
        else:
            self.fractionSpinBox.setEnabled(True)
            self.fractionsLabel.setEnabled(True)

        self.fullDisplayRadioButton.blockSignals(False)

    def action(self):
        print(f"Process FC Meshes to GDML Tessellations")
        self.accept()

    def onCancel(self):
        print(f"reject")
        self.reject()

    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("Mesh2TessellateDialog", "Mesh2Tess"))
        self.textEdit.setHtml(
            _translate(
                "Mesh2TessellateDialog",
                '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0//EN" "http://www.w3.org/TR/REC-html40/strict.dtd">\n'
                '<html><head><meta name="qrichtext" content="1" /><style type="text/css">\n'
                "p, li { white-space: pre-wrap; }\n"
                "</style></head><body style=\" font-family:'Noto Sans'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
                '<p style=" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;">Convert a mesh to a GDML Tessellated solid. For very large meshes, the building of a full solid in python might take a <span style=" font-weight:600;">very long</span> time. To speed up conversion, a sampling of the facets can be displayed, instead of creating a full solid. <span style=" font-weight:600;">On export to GDML, the full solid will be exported</span>. The fraction of faces displayed can be later changed in the Properties paenl.</p></body></html>',
            )
        )

        self.groupBox.setTitle(
            _translate("Mesh2TessellateDialog", "Tessellation display")
        )
        self.fullDisplayRadioButton.setToolTip(
            _translate("Mesh2TessellateDialog", "Display full solid")
        )
        self.fullDisplayRadioButton.setText(
            _translate("Mesh2TessellateDialog", "Full solid")
        )
        self.samplesRadioButton.setToolTip(
            _translate("Mesh2TessellateDialog", "Sample facets")
        )
        self.samplesRadioButton.setText(
            _translate("Mesh2TessellateDialog", "Samples only")
        )
        self.fractionsLabel.setText(
            _translate("Mesh2TessellateDialog", "Sampled fraction (%)")
        )
        self.fractionSpinBox.setToolTip(
            _translate("Mesh2TessellateDialog", "Percent of facets sampled")
        )

def getSelectedMaterial():
    from .exportGDML import nameFromLabel
    from .GDMLObjects import GDMLmaterial

    list = FreeCADGui.Selection.getSelection()
    if list is not None:
        for obj in list:
            if hasattr(obj, 'Proxy'):
                if isinstance(obj.Proxy, GDMLmaterial) is True:
                    return nameFromLabel(obj.Label)

    return 0


def getVert(s):
    if '/' in s:
        ret = int(s[:s.index('/')])-1
        # print(ret)
    else:
        ret = int(s)-1
    return(ret)


def findNearestPart(obj):
    if obj.TypeId == "App::Part":
        return obj
    for o in obj.InList:
        return findNearestPart(o)   

def processOBJ(doc, filename):

    import FreeCADGui, Mesh, re
    from .GDMLObjects import checkMaterialDefinitionsExist
    from datetime import datetime

    from .GDMLObjects import GDMLTessellated, ViewProvider
    from .GDMLCommands import Mesh2TessDialog

    print("Check Materials definitions exist")
    checkMaterialDefinitionsExist()
    print("import OBJ as GDML Tessellated")
    startTime = datetime.now()
    mapDialog = MapObjmat2GDMLmatDialog()
    #mapDialog.initMaterials()
    #mapDialog.exec_()
    #return
    # Preprocess file collecting Object and Material definitions
    fp = pythonopen(filename)
    data = fp.read()
    pattern = re.compile(r"^(?:[0g]|usemtl)\s.*", re.MULTILINE)
    objMatList = pattern.findall(data)
    #print(f"Obj Mat List {objMatList}")
    objDict = {}
    name = ""
    objMat = ""
    for i in objMatList:
        if i[:2] == 'g ':
            name = i.lstrip('g ')
            #print(f"Name {name}")
        if i[:7] == 'usemtl ':
            objMat = i.lstrip('usemtl ')
            #print(f"Material {objMat}")
            objDict[name] = objMat
            mapDialog.addMaterialMapping(name, objMat, "G4_A-150_TISSUE")
    #print(f"Obj Dict {objDict}")
    mapDialog.exec_()
    preTime = datetime.now()
    print(f"Time for preprocess objects materials {preTime - startTime}")
    # Read OBJ file using FC mesh
    meshDoc = FreeCAD.newDocument("TempObj")
    #print(f"Active document {FreeCADGui.ActiveDocument.Document.Name}")
    #Mesh.open(filename)
    Mesh.insert(filename)
    fcMeshTime = datetime.now()
    print(f"Time for FC mesh load of OBJ file {fcMeshTime - preTime}")
    FreeCADGui.Selection.clearSelection()
    #for obj in doc.Objects:
    #    FreeCADGui.Selection.addSelection(obj)
    FreeCADGui.Selection.addSelection(meshDoc.Objects[0])
    sel = FreeCADGui.Selection.getSelection()
    print(f"Selection {sel}")
    dialog = Mesh2TessDialog(sel, doc)
    dialog.exec_()
    FreeCADGui.setActiveDocument(doc)
    FreeCAD.ActiveDocument.recompute()
    if FreeCAD.GuiUp:
        FreeCADGui.SendMsgToActiveView("ViewFit")
    FreeCAD.Console.PrintMessage("End import OBJ file\n")
