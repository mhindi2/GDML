# Sun Feb 20 12:50:19 PM PST 2022
# Fri Feb 11 12:09:03 PM PST 2022
# **************************************************************************
# *                                                                        *
# *   Copyright (c) 2019 Keith Sloan <keith@sloan-home.co.uk>              *
# *             (c) 2020 Dam Lambert                                       *
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
# *   Acknowledgements : Ideas & code copied from                          *
# *                      https://github.com/ignamv/geanTipi                *
# *                                                                        *
# ***************************************************************************
__title__ = "FreeCAD - GDML exporter Version"
__author__ = "Keith Sloan <keith@sloan-home.co.uk>"
__url__ = ["https://github.com/KeithSloan/FreeCAD_Geant4"]

import FreeCAD, os, Part, math
from FreeCAD import Vector
from .GDMLObjects import GDMLcommon, GDMLBox, GDMLTube
# modif add
# from .GDMLObjects import getMult, convertionlisteCharToLunit

import sys
try:
    import lxml.etree as ET
    FreeCAD.Console.PrintMessage("running with lxml.etree\n")
    XML_IO_VERSION = 'lxml'
except ImportError:
    try:
        import xml.etree.ElementTree as ET
        FreeCAD.Console.PrintMessage("running with xml.etree.ElementTree\n")
        XML_IO_VERSION = 'xml'
    except ImportError:
        FreeCAD.Console.PrintMessage('pb xml lib not found\n')
        sys.exit()
# xml handling
# import argparse
# from   xml.etree.ElementTree import XML
#################################

global zOrder

from .GDMLObjects import GDMLQuadrangular, GDMLTriangular, \
                        GDML2dVertex, GDMLSection, \
                        GDMLmaterial, GDMLfraction, \
                        GDMLcomposite, GDMLisotope, \
                        GDMLelement, GDMLconstant, GDMLvariable, \
                        GDMLquantity

from . import GDMLShared

# ***************************************************************************
# Tailor following to your requirements ( Should all be strings )          *
# no doubt there will be a problem when they do implement Value
if open.__module__ in ['__builtin__', 'io']:
    pythonopen = open  # to distinguish python built-in open function from the one declared here

# ## modifs lambda


def verifNameUnique(name):
    # need to be done!!
    return True

# ## end modifs lambda

#################################
# Switch functions
################################


class switch(object):
    value = None

    def __new__(class_, value):
        class_.value = value
        return True


def case(*args):
    return any((arg == switch.value for arg in args))


class MultiPlacer:
    def __init__(self, obj):
        self.obj = obj

    def place(self, volRef):
        print("Can't place base class MultiPlace")

    def xml(self):
        print("Can't place base class MultiPlace")

    def name(self):
        return self.obj.Label

    @staticmethod
    def getPlacer(obj):
        if obj.TypeId == 'Part::Mirroring':
            return MirrorPlacer(obj)
        else:
            print(f'{obj.Label} is not a placer')
            return None


class MirrorPlacer(MultiPlacer):
    def __init__(self, obj):
        super().__init__(obj)

    def place(self, volRef):
        global structure
        assembly = ET.Element('assembly', {'name': self.obj.Label})
        # structure.insert(0, assembly)
        # insert just before worlVol, which should be last
        worldIndex = len(structure) - 1
        structure.insert(worldIndex, assembly)
        pos = self.obj.Source.Placement.Base
        name = volRef
        pvol = ET.SubElement(assembly, 'physvol')
        ET.SubElement(pvol, 'volumeref', {'ref': volRef})
        exportPosition(name, pvol, pos)

        name = volRef+'_mirror'
        pvol = ET.SubElement(assembly, 'physvol')
        ET.SubElement(pvol, 'volumeref', {'ref': volRef})
        normal = self.obj.Normal
        # reflect the position about the reflection plane
        unitVec = normal.normalize()
        posAlongNormal = pos.dot(unitVec)*unitVec
        posParalelToPlane = pos - posAlongNormal
        newPos = posParalelToPlane - posAlongNormal
        exportPosition(name, pvol, newPos)
        # first reflect about x-axis
        # then rotate to bring x-axis to direction of normal
        rotX = False
        if normal.x == 1:
            scl = Vector(-1, 1, 1)
        elif normal.y == 1:
            scl = Vector(1, -1, 1)
        elif normal.z == 1:
            scl = Vector(1, 1, 1-1)
        else:
            scl = Vector(-1, 1, 1)
            rotX = True
        exportScaling(name, pvol, scl)
        if rotX is True:
            # rotation to bring normal to x-axis (might have to reverse)
            rot = FreeCAD.Rotation(normal, Vector(1, 0, 0))
            exportRotation(name, pvol, rot)


#########################################################
# Pretty format GDML                                    #
#########################################################


def indent(elem, level=0):
    i = "\n" + level*"  "
    j = "\n" + (level-1)*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for subelem in elem:
            indent(subelem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = j
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = j
    return elem

#########################################


def nameFromLabel(label):
    if ' ' not in label:
        return label
    else:
        return(label.split(' ')[0])


def initGDML():
    NS = 'http://www.w3.org/2001/XMLSchema-instance'
    location_attribute = '{%s}noNamespaceSchemaLocation' % NS
    gdml = ET.Element('gdml', attrib={location_attribute:
                                      'http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd'})
    # print(gdml.tag)

    return gdml

#################################
#  Setup GDML environment
#################################


def GDMLstructure():
    # print("Setup GDML structure")
    #################################
    # globals
    ################################
    global gdml, constants, variables, define, materials, solids, \
           structure, setup
    global WorldVOL
    global defineCnt, LVcount, PVcount, POScount, ROTcount, SCLcount
    global centerDefined
    global identityDefined
    global gxml

    centerDefined = False
    identityDefined = False
    defineCnt = LVcount = PVcount = POScount = ROTcount = SCLcount = 1

    gdml = initGDML()
    define = ET.SubElement(gdml, 'define')
    materials = ET.SubElement(gdml, 'materials')
    solids = ET.SubElement(gdml, 'solids')
    solids.clear()
    structure = ET.SubElement(gdml, 'structure')
    setup = ET.SubElement(gdml, 'setup', {'name': 'Default', 'version': '1.0'})
    gxml = ET.Element('gxml')

    return structure


def defineMaterials():
    # Replaced by loading Default
    # print("Define Materials")
    global materials


def exportDefine(name, v):
    global define
    ET.SubElement(define, 'position', {'name': name, 'unit': 'mm',
                                       'x': str(v[0]),
                                       'y': str(v[1]),
                                       'z': str(v[2])})


def exportDefineVertex(name, v, index):
    global define
    ET.SubElement(define, 'position', {'name': name + str(index),
                                       'unit': 'mm',
                                       'x': str(v.X),
                                       'y': str(v.Y),
                                       'z': str(v.Z)})


def defineWorldBox(bbox):
    global solids
    for obj in FreeCAD.ActiveDocument.Objects:
        # print("{} + {} = ".format(bbox, obj.Shape.BoundBox))
        if hasattr(obj, "Shape"):
            bbox.add(obj.Shape.BoundBox)
        if hasattr(obj, "Mesh"):
            bbox.add(obj.Mesh.BoundBox)
        if hasattr(obj, "Points"):
            bbox.add(obj.Points.BoundBox)
    #   print(bbox)
    # Solids get added to solids section of gdml ( solids is a global )
    name = 'WorldBox'
    ET.SubElement(solids, 'box', {'name': name,
                                  'x': str(1000),
                                  'y': str(1000),
                                  'z': str(1000),
                                  # 'x': str(2*max(abs(bbox.XMin), abs(bbox.XMax))), \
                                  # 'y': str(2*max(abs(bbox.YMin), abs(bbox.YMax))), \
                                  # 'z': str(2*max(abs(bbox.ZMin), abs(bbox.ZMax))), \
                                  'lunit': 'mm'})
    return(name)


def quaternion2XYZ(rot):
    '''
    convert a quaternion rotation to a sequence of rotations around X, Y, Z
    Here is my (Munther Hindi) derivation:
    First, the rotation matrices for rotations around X, Y, Z axes.

    Rx = [      1       0       0]
        [      0  cos(a) -sin(a)]
        [      0  sin(a)  cos(a)]

    Ry= [ cos(b)       0  sin(b)]
       [      0       1       0]
       [-sin(b)       0   cos(b)]

    Rz = [ cos(g) -sin(g)       0]
        [ sin(g)  cos(g)       0]
        [      0       0       1]

    R = Rz Ry Rx = [
    [cos(b)*cos(g), cos(g)*sin(a)*sin(b) - cos(a)*sin(g), cos(a)*cos(g)*sin(b) + sin(a)*sin(g)]
    [cos(b)*sin(g), sin(a)*sin(b)*sin(g) + cos(a)*cos(g), cos(a)*sin(b)*sin(g) - cos(g)*sin(a)]
    [-sin(b),        cos(b)*sin(a),                       ,  cos(a)*cos(b)]]

    To get the angles b(eta), g(amma) for rotations around y, z axes, transform the unit vector (1,0,0)
    [x,y,z] = Q*(1,0,0) = R*(1,0,0) ==>
    x = cos(b)*cos(g)
    y = cos(b)*sin(g)
    z = -sin(b)

    ==>   g = atan2(y, x) = atan2(sin(b)*sin(g), cos(g)*sin(b)) = atan2(sin(g), cos(g))
    then  b = atan2(-z*cos(g), x) = atan2(sin(b)*cos(g), cos(b)*cos(g)] = atan2(sin(b), cos(b))

    Once b, g are found, a(lpha) can be found by transforming (0, 0, 1), or (0,1,0)
    Since now g, b are known, one can form the inverses of Ry, Rz:
    Rz^-1 = Rz(-g)
    Ry^-1 = Ry(-b)

    Now R*(0,0,1) = Rz*Ry*Rz(0,1,0) = (x, y, z)
    multiply both sides by Ry^-1 Rz^-1:
    Ry^-1 Rz^-1 Rz Ry Rx (0,1,0) = Rx (0,1,0) = Ry(-b) Rz(-g) (x, y, z) = (xp, yp, zp)
    ==>
    xp = 0
    yp = cos(a)
    zp = sin(a)

    and a = atan2(zp, yp)
    '''
    v = rot*Vector(1, 0, 0)
    g = math.atan2(v.y, v.x)
    b = math.atan2(-v.z*math.cos(g), v.x)
    Ryinv = FreeCAD.Rotation(Vector(0, 1, 0), math.degrees(-b))
    Rzinv = FreeCAD.Rotation(Vector(0, 0, 1), math.degrees(-g))
    vp = Ryinv*Rzinv*rot*Vector(0, 1, 0)
    a = math.atan2(vp.z, vp.y)

    print([math.degrees(a), math.degrees(b), math.degrees(g)])
    return [math.degrees(a), math.degrees(b), math.degrees(g)]


def createLVandPV(obj, name, solidName):
    #
    # Logical & Physical Volumes get added to structure section of gdml
    #
    # Need to update so that use export of Rotation & position
    # rather than this as well i.e one Place
    #
    global PVcount, POScount, ROTcount
    pvName = 'PV'+name+str(PVcount)
    PVcount += 1
    pos = obj.Placement.Base
    lvol = ET.SubElement(structure, 'volume', {'name': pvName})
    material = getMaterial(obj)
    ET.SubElement(lvol, 'materialref', {'ref': material})
    ET.SubElement(lvol, 'solidref', {'ref': solidName})
    # Place child physical volume in World Volume
    phys = ET.SubElement(lvol, 'physvol', {'name': 'PV-'+name})
    ET.SubElement(phys, 'volumeref', {'ref': pvName})
    x = pos[0]
    y = pos[1]
    z = pos[2]
    if x != 0 and y != 0 and z != 0:
        posName = 'Pos'+name+str(POScount)
        POScount += 1
        ET.SubElement(phys, 'positionref', {'name': posName})
        ET.SubElement(define, 'position', {'name': posName, 'unit': 'mm',
                                           'x': str(x), 'y': str(y), 'z': str(z)})
    # Realthunders enhancement to toEuler ixyz is intrinsic
    rot = obj.Placement.Rotation
    if hasattr(rot, 'toEulerAngles'):
        angles = rot.toEulerAngles('ixyz')
        angles = (angles[2], angles[1], angles[0])
    else:
        print('Export of rotation probably wrong')
        print('Needs toEulerAngles function - Use LinkStage 3')
        angles = rot.toEuler()
    GDMLShared.trace("Angles")
    GDMLShared.trace(angles)
    angles = quaternion2XYZ(rot)
    a0 = angles[0]
    # print(a0)
    a1 = angles[1]
    # print(a1)
    a2 = angles[2]
    # print(a2)
    if a0 != 0 and a1 != 0 and a2 != 0:
        rotName = 'Rot'+name+str(ROTcount)
        ROTcount += 1
        ET.SubElement(phys, 'rotationref', {'name': rotName})
        ET.SubElement(define, 'rotation', {'name': rotName, 'unit': 'deg',
                                           'x': str(-a0), 'y': str(-a1), 'z': str(-a2)})


def createAdjustedLVandPV(obj, name, solidName, delta):
    # Allow for difference in placement between FreeCAD and GDML
    adjObj = obj
    rot = FreeCAD.Rotation(obj.Placement.Rotation)
    adjObj.Placement.move(rot.multVec(delta))  # .negative()
    createLVandPV(adjObj, name, solidName)


def reportObject(obj):

    GDMLShared.trace("Report Object")
    GDMLShared.trace(obj)
    GDMLShared.trace("Name : "+obj.Name)
    GDMLShared.trace("Type : "+obj.TypeId)
    if hasattr(obj, 'Placement'):
        print("Placement")
        print("Pos   : "+str(obj.Placement.Base))
        print("axis  : "+str(obj.Placement.Rotation.Axis))
        print("angle : "+str(obj.Placement.Rotation.Angle))

    while switch(obj.TypeId):

        ###########################################
        # FreeCAD GDML Parts                      #
        ###########################################
        if case("Part::FeaturePython"):
            GDMLShared.trace("Part::FeaturePython")
            #
            # if hasattr(obj.Proxy,'Type'):
            #   print (obj.Proxy.Type)
            #   print (obj.Name)
            # else :
            #   print("Not a GDML Feature")

            # print dir(obj)
            # print dir(obj.Proxy)
            # print("cylinder : Height "+str(obj.Height)+ " Radius "+str(obj.Radius))
            break
        ###########################################
        # FreeCAD Parts                           #
        ###########################################
        if case("Part::Sphere"):
            print("Sphere Radius : "+str(obj.Radius))
            break

        if case("Part::Box"):
            print("cube : (" + str(obj.Length)+"," + str(obj.Width)+"," +
                  str(obj.Height) + ")")
            break

        if case("Part::Cylinder"):
            print("cylinder : Height " + str(obj.Height) +
                  " Radius "+str(obj.Radius))
            break

        if case("Part::Cone"):
            print("cone : Height "+str(obj.Height) +
                  " Radius1 "+str(obj.Radius1)+" Radius2 "+str(obj.Radius2))
            break

        if case("Part::Torus"):
            print("Torus")
            print(obj.Radius1)
            print(obj.Radius2)
            break

        if case("Part::Prism"):
            print("Prism")
            break

        if case("Part::RegularPolygon"):
            print("RegularPolygon")
            break

        if case("Part::Extrusion"):
            print("Extrusion")
            break

        if case("Circle"):
            print("Circle")
            break

        if case("Extrusion"):
            print("Wire extrusion")
            break

        if case("Mesh::Feature"):
            print("Mesh")
            # print dir(obj.Mesh)
            break

        print("Other")
        print(obj.TypeId)
        break


def processPlanar(obj, shape, name):
    print('Polyhedron ????')
    global defineCnt
    #
    # print("Add tessellated Solid")
    tess = ET.SubElement(solids, 'tessellated', {'name': name})
    # print("Add Vertex positions")
    for f in shape.Faces:
        baseVrt = defineCnt
        for vrt in f.Vertexes:
            vnum = 'v'+str(defineCnt)
            ET.SubElement(define, 'position', {'name': vnum,
                                               'x': str(vrt.Point.x),
                                               'y': str(vrt.Point.y),
                                               'z': str(vrt.Point.z),
                                               'unit': 'mm'})
            defineCnt += 1
        # print("Add vertex to tessellated Solid")
        vrt1 = 'v'+str(baseVrt)
        vrt2 = 'v'+str(baseVrt+1)
        vrt3 = 'v'+str(baseVrt+2)
        vrt4 = 'v'+str(baseVrt+3)
        NumVrt = len(f.Vertexes)
        if NumVrt == 3:
            ET.SubElement(tess, 'triangular', {
               'vertex1': vrt1,
               'vertex2': vrt2,
               'vertex3': vrt3,
               'type': 'ABSOLUTE'})
        elif NumVrt == 4:
            ET.SubElement(tess, 'quadrangular', {
               'vertex1': vrt1,
               'vertex2': vrt2,
               'vertex3': vrt3,
               'vertex4': vrt4,
               'type': 'ABSOLUTE'})


def checkShapeAllPlanar(Shape):
    for f in Shape.Faces:
        if f.Surface.isPlanar() is False:
            return False
    return True


#    Add XML for TessellateSolid
def mesh2Tessellate(mesh, name):
    global defineCnt

    baseVrt = defineCnt
    # print ("mesh")
    # print (mesh)
    # print ("Facets")
    # print (mesh.Facets)
    # print ("mesh topology")
    # print (dir(mesh.Topology))
    # print (mesh.Topology)
    #
    #    mesh.Topology[0] = points
    #    mesh.Topology[1] = faces
    #
    #    First setup vertex in define section vetexs (points)
    # print("Add Vertex positions")
    for fc_points in mesh.Topology[0]:
        # print(fc_points)
        v = 'v'+str(defineCnt)
        ET.SubElement(define, 'position', {'name': v,
                                           'x': str(fc_points[0]),
                                           'y': str(fc_points[1]),
                                           'z': str(fc_points[2]),
                                           'unit': 'mm'})
        defineCnt += 1
    #
    #     Add faces
    #
    # print("Add Triangular vertex")
    tess = ET.SubElement(solids, 'tessellated', {'name': name})
    for fc_facet in mesh.Topology[1]:
        # print(fc_facet)
        vrt1 = 'v'+str(baseVrt+fc_facet[0])
        vrt2 = 'v'+str(baseVrt+fc_facet[1])
        vrt3 = 'v'+str(baseVrt+fc_facet[2])
        ET.SubElement(tess, 'triangular', {
           'vertex1': vrt1, 'vertex2': vrt2,
           'vertex3': vrt3, 'type': 'ABSOLUTE'})


def processMesh(obj, Mesh, Name):
    #  obj needed for Volune names
    #  object maynot have Mesh as part of Obj
    #  Name - allows control over name
    print("Create Tessellate Logical Volume")
    createLVandPV(obj, Name, 'Tessellated')
    mesh2Tessellate(Mesh, Name)
    return(Name)


def shape2Mesh(shape):
    import MeshPart
    return (MeshPart.meshFromShape(Shape=shape, Deflection=0.0))
#            Deflection= params.GetFloat('meshdeflection',0.0))


def processObjectShape(obj):
    # Check if Planar
    # If plannar create Tessellated Solid with 3 & 4 vertex as appropriate
    # If not planar create a mesh and the a Tessellated Solid with 3 vertex
    # print("Process Object Shape")
    # print(obj)
    # print(obj.PropertiesList)
    if not hasattr(obj, 'Shape'):
        return
    shape = obj.Shape
    # print (shape)
    # print(shape.ShapeType)
    while switch(shape.ShapeType):
        if case("Mesh::Feature"):
            print("Mesh - Should not occur should have been handled")
            # print("Mesh")
            # tessellate = mesh2Tessellate(mesh)
            # return(tessellate)
            # break

            print("ShapeType Not handled")
            print(shape.ShapeType)
            break

#   Dropped through to here
#   Need to check has Shape

    # print('Check if All planar')
    planar = checkShapeAllPlanar(shape)
    # print(planar)

    if planar:
        return(processPlanar(obj, shape, obj.Name))

    else:
        # Create Mesh from shape & then Process Mesh
        # to create Tessellated Solid in Geant4
        return(processMesh(obj, shape2Mesh(shape), obj.Name))


def processSection(obj):
    # print("Process Section")
    ET.SubElement(solids, 'section', {
       'vertex1': obj.v1,
       'vertex2': obj.v2,
       'vertex3': obj.v3, 'vertex4': obj.v4,
       'type': obj.vtype})


def addPhysVol(xmlVol, volName):
    GDMLShared.trace("Add PhysVol to Vol : " + volName)
    # print(ET.tostring(xmlVol))
    pvol = ET.SubElement(xmlVol, 'physvol', {'name': 'PV-'+volName})
    ET.SubElement(pvol, 'volumeref', {'ref': volName})
    return pvol


def cleanVolName(obj, volName):
    # Get proper Volume Name
    # print('clean name : '+volName)
    if hasattr(obj, 'Copynumber'):
        # print('Has copynumber')
        i = len(volName)
        if '_' in volName and i > 2:
            volName = volName[:-2]
    # print('returning name : '+volName)
    return volName


def addPhysVolPlacement(obj, xmlVol, volName, placement):
    # obj: App:Part to be placed.
    # xmlVol: the xml that the <physvol is a subelement of.
    # It may be a <volume, or an <assembly
    # volref: the name of the volume being placed
    # placement: the placement of the <physvol
    # For most situations, the placement (pos, rot) should be that
    # of the obj (obj.Placement.Base, obj.Placement.Rotation), but
    # if the user specifies a placement for the solid, then the palcement
    # has to be a product of both placements. Here we don't try to figure
    # that out, so we demand the placement be given explicitly

    # returns xml of of reated <physvol element

    # Get proper Volume Name
    # I am commenting this out I don't know why it's needed.
    # the <volume or <assembly name is ceated withoutout any cleanup,m so the
    # reference to it musl also not have any cleanup
    # refName = cleanVolName(obj, volName)
    # GDMLShared.setTrace(True)
    GDMLShared.trace("Add PhysVol to Vol : "+volName)
    # print(ET.tostring(xmlVol))
    if xmlVol is not None:
        if not hasattr(obj, 'CopyNumber'):
            pvol = ET.SubElement(xmlVol, 'physvol', {'name': 'PV-' + volName})
        else:
            cpyNum = str(obj.CopyNumber)
            GDMLShared.trace('CopyNumber : '+cpyNum)
            pvol = ET.SubElement(xmlVol, 'physvol', {'copynumber': cpyNum})

        ET.SubElement(pvol, 'volumeref', {'ref': volName})
        processPlacement(volName, pvol, placement)
        if hasattr(obj, 'GDMLscale'):
            scaleName = volName+'scl'
            ET.SubElement(pvol, 'scale', {'name': scaleName,
                                          'x': str(obj.GDMLscale[0]),
                                          'y': str(obj.GDMLscale[1]),
                                          'z': str(obj.GDMLscale[2])})

        return pvol


def exportPosition(name, xml, pos):
    global POScount
    global centerDefined
    GDMLShared.trace('export Position')
    GDMLShared.trace(pos)
    x = pos[0]
    y = pos[1]
    z = pos[2]
    if x == 0 and y == 0 and z == 0:
        if not centerDefined:
            centerDefined = True
            ET.SubElement(define, 'position', {'name': 'center',
                                               'x': '0', 'y': '0', 'z': '0',
                                               'unit': 'mm'})
        ET.SubElement(xml, 'positionref', {'ref': 'center'})

    else:
        posName = 'P-' + name + str(POScount)
        POScount += 1
        posxml = ET.SubElement(define, 'position', {'name': posName,
                                                    'unit': 'mm'})
        if x != 0:
            posxml.attrib['x'] = str(x)
        if y != 0:
            posxml.attrib['y'] = str(y)
        if z != 0:
            posxml.attrib['z'] = str(z)
        ET.SubElement(xml, 'positionref', {'ref': posName})


def exportRotation(name, xml, rot):
    print('Export Rotation')
    global ROTcount
    global identityDefined
    if rot.Angle == 0:
        if not identityDefined:
            identityDefined = True
            ET.SubElement(define, 'rotation', {'name': 'identity',
                                               'x': '0', 'y': '0', 'z': '0'})
        ET.SubElement(xml, 'rotationref', {'ref': 'identity'})

    else:

        # Realthunders enhancement to toEuler ixyz is intrinsic
        if hasattr(rot, 'toEulerAngles'):
            angles = rot.toEulerAngles('ixyz')
            angles = (angles[2], angles[1], angles[0])
        else:
            print('Export of rotation probably wrong')
            print('Needs toEulerAngles function - Use LinkStage 3')
            angles = rot.toEuler()
        GDMLShared.trace("Angles")
        GDMLShared.trace(angles)
        a0 = angles[0]
        print(a0)
        a1 = angles[1]
        print(a1)
        a2 = angles[2]
        print(a2)
        if a0 != 0 or a1 != 0 or a2 != 0:
            rotName = 'R-'+name+str(ROTcount)
            ROTcount += 1
            rotxml = ET.SubElement(define, 'rotation', {'name': rotName,
                                                        'unit': 'deg'})
            if abs(a2) != 0:
                rotxml.attrib['x'] = str(-a2)
            if abs(a1) != 0:
                rotxml.attrib['y'] = str(-a1)
            if abs(a0) != 0:
                rotxml.attrib['z'] = str(-a0)
            ET.SubElement(xml, 'rotationref', {'ref': rotName})


def exportScaling(name, xml, scl):
    global SCLcount
    global centerDefined
    GDMLShared.trace('export Scaling')
    GDMLShared.trace(scl)
    x = scl[0]
    y = scl[1]
    z = scl[2]
    sclName = 'S-' + name + str(SCLcount)
    SCLcount += 1
    ET.SubElement(define, 'scale', {'name': sclName,
                                    'x': str(x),
                                    'y': str(y),
                                    'z': str(z)})
    ET.SubElement(xml, 'scaleref', {'ref': sclName})


def processPlacement(name, xml, placement):
    exportPosition(name, xml, placement.Base)
    exportRotation(name, xml, placement.Rotation)


def processPosition(obj, solid):
    if obj.Placement.Base == FreeCAD.Vector(0, 0, 0):
        return
    GDMLShared.trace("Define position & references to Solid")
    exportPosition(obj.Name, solid, obj.Placement.Base)


def processRotation(obj, solid):
    if obj.Placement.Rotation.Angle == 0:
        return
    GDMLShared.trace('Deal with Rotation')
    exportRotation(obj.Name, solid, obj.Placement.Rotation)


def testDefaultPlacement(obj):
    # print(dir(obj.Placement.Rotation))
    # print('Test Default Placement : '+obj.Name)
    # print(obj.Placement.Base)
    # print(obj.Placement.Rotation.Angle)
    if obj.Placement.Base == FreeCAD.Vector(0, 0, 0) and \
       obj.Placement.Rotation.Angle == 0:
        return True
    else:
        return False


def testAddPhysVol(obj, xmlParent, volName):
    if testDefaultPlacement(obj) is False:
        if xmlParent is not None:
            pvol = addPhysVol(xmlParent, volName)
            processPlacement(obj, pvol)
        else:
            print('Root/World Volume')


def addVolRef(volxml, volName, obj, solidName=None):
    # Pass material as Boolean
    material = getMaterial(obj)
    if solidName is None:
        solidName = nameOfGDMLobject(obj)
    ET.SubElement(volxml, 'materialref', {'ref': material})
    ET.SubElement(volxml, 'solidref', {'ref': solidName})

    ET.SubElement(gxml, 'volume', {'name': volName, 'material': material})

    if hasattr(obj.ViewObject, 'ShapeColor') and volName != WorldVOL:
        colour = obj.ViewObject.ShapeColor
        colStr = '#'+''.join('{:02x}'.format(round(v*255)) for v in colour)
        ET.SubElement(volxml, 'auxiliary', {'auxtype': 'Color',
                                            'auxvalue': colStr})
    # print(ET.tostring(volxml))


def nameOfGDMLobject(obj):
    name = obj.Label
    if len(name) > 4:
        if name[0:4] == 'GDML':
            if '_' in name:
                return(name.split('_', 1)[1])
    return name


def processIsotope(obj, item):  # maybe part of material or element (common code)
    if hasattr(obj, 'Z'):
        # print(dir(obj))
        item.set('Z', str(obj.Z))

    if hasattr(obj, 'N'):
        # print(dir(obj))
        item.set('N', str(obj.N))

    if hasattr(obj, 'formula'):
        # print(dir(obj))
        item.set('formula', str(obj.formula))

    if hasattr(obj, 'unit') or hasattr(obj, 'atom_value') or \
       hasattr(obj, 'value'):
        atom = ET.SubElement(item, 'atom')

        if hasattr(obj, 'unit'):
            atom.set('unit', str(obj.unit))

        if hasattr(obj, 'type'):
            atom.set('unit', str(obj.type))

        if hasattr(obj, 'atom_value'):
            atom.set('value', str(obj.atom_value))

        if hasattr(obj, 'value'):
            atom.set('value', str(obj.value))


def processMaterials():
    print("\nProcess Materials")
    global materials

    for GName in ['Constants', 'Variables',
                  'Isotopes', 'Elements', 'Materials']:
        Grp = FreeCAD.ActiveDocument.getObject(GName)
        if Grp is not None:
            # print(Grp.TypeId+" : "+Grp.Name)
            print(Grp.Label)
            if processGroup(Grp) is False:
                break


def processFractionsComposites(obj, item):
    # Fractions are used in Material and Elements
    if isinstance(obj.Proxy, GDMLfraction):
        # print("GDML fraction :" + obj.Label)
        # need to strip number making it unique
        ET.SubElement(item, 'fraction', {'n': str(obj.n),
                                         'ref': nameFromLabel(obj.Label)})

    if isinstance(obj.Proxy, GDMLcomposite):
        # print("GDML Composite")
        ET.SubElement(item, 'composite', {'n': str(obj.n),
                                          'ref': nameFromLabel(obj.Label)})


def createMaterials(group):
    global materials
    for obj in group:
        if obj.Label != 'Geant4':
            item = ET.SubElement(materials, 'material', {
               'name':
               nameFromLabel(obj.Label)})

            # Dunit & Dvalue must be first for Geant4
            if hasattr(obj, 'Dunit') or hasattr(obj, 'Dvalue'):
                # print("Dunit or DValue")
                D = ET.SubElement(item, 'D')
                if hasattr(obj, 'Dunit'):
                    D.set('unit', str(obj.Dunit))

                if hasattr(obj, 'Dvalue'):
                    D.set('value', str(obj.Dvalue))

                if hasattr(obj, 'Tunit') and hasattr(obj, 'Tvalue'):
                    ET.SubElement(item, 'T', {'unit': obj.Tunit,
                                              'value': str(obj.Tvalue)})

                if hasattr(obj, 'MEEunit'):
                    ET.SubElement(item, 'MEE', {'unit': obj.MEEunit,
                                                'value': str(obj.MEEvalue)})
            # process common options material / element
            processIsotope(obj, item)
            if len(obj.Group) > 0:
                for o in obj.Group:
                    processFractionsComposites(o, item)


def createElements(group):
    global materials
    for obj in group:
        # print(f'Element : {obj.Label}')
        item = ET.SubElement(materials, 'element', {
           'name': nameFromLabel(obj.Label)})
        # Common code IsoTope and Elements1
        processIsotope(obj, item)

        if len(obj.Group) > 0:
            for o in obj.Group:
                processFractionsComposites(o, item)


def createConstants(group):
    global define
    for obj in group:
        if isinstance(obj.Proxy, GDMLconstant):
            # print("GDML constant")
            # print(dir(obj))

            ET.SubElement(define, 'constant', {'name': obj.Label,
                                               'value': obj.value})


def createVariables(group):
    global define
    for obj in group:
        if isinstance(obj.Proxy, GDMLvariable):
            # print("GDML variable")
            # print(dir(obj))

            ET.SubElement(define, 'variable', {'name': obj.Label,
                                               'value': obj.value})


def createQuantities(group):
    global define
    for obj in group:
        if isinstance(obj.Proxy, GDMLquantity):
            # print("GDML quantity")
            # print(dir(obj))

            ET.SubElement(define, 'quantity', {'name': obj.Label,
                                               'type': obj.type,
                                               'unit': obj.unit,
                                               'value': obj.value})


def createIsotopes(group):
    global materials
    for obj in group:
        if isinstance(obj.Proxy, GDMLisotope):
            # print("GDML isotope")
            # item = ET.SubElement(materials,'isotope',{'N': str(obj.N), \
            #                                           'Z': str(obj.Z), \
            #                                           'name' : obj.Label})
            # ET.SubElement(item,'atom',{'unit': obj.unit, \
            #                           'value': str(obj.value)})
            item = ET.SubElement(materials, 'isotope', {'name': obj.Label})
            processIsotope(obj, item)


def processGroup(obj):
    print('Process Group '+obj.Label)
    # print(obj.TypeId)
    # print(obj.Group)
    #      if hasattr(obj,'Group') :
    # return
    if hasattr(obj, 'Group'):
        # print("   Object List : "+obj.Label)
        # print(obj)
        while switch(obj.Name):
            if case("Constants"):
                # print("Constants")
                createConstants(obj.Group)
                break

            if case("Variables"):
                # print("Variables")
                createVariables(obj.Group)
                break

            if case("Quantities"):
                # print("Quantities")
                createQuantities(obj.Group)
                break

            if case("Isotopes"):
                # print("Isotopes")
                createIsotopes(obj.Group)
                break

            if case("Elements"):
                # print("Elements")
                createElements(obj.Group)
                break

            if case("Materials"):
                print("Materials")
                createMaterials(obj.Group)
                break

            if case("Geant4"):
                # Do not export predefine in Geant4
                print("Geant4")
                break


from itertools import islice


def consume(iterator):
    next(islice(iterator, 2, 2), None)


def getXmlVolume(volObj):
    global structure
    if volObj is None:
        return None
    xmlvol = structure.find("volume[@name='%s']" % volObj.Label)
    if xmlvol is None:
        print(volObj.Label+' Not Found')
    return xmlvol


def getDefaultMaterial():
    # should get this from GDML settings under "settings"
    # for now this does not exist, so simply put steel
    return 'G4_STAINLESS-STEEL'


def getBooleanCount(obj):
    GDMLShared.trace('get Count : ' + obj.Name)
    if hasattr(obj, 'Tool'):
        GDMLShared.trace('Has tool - check Base')
        baseCnt = getBooleanCount(obj.Base)
        toolCnt = getBooleanCount(obj.Tool)
        GDMLShared.trace('Count is : ' + str(baseCnt + toolCnt))
        return (baseCnt + toolCnt)
    else:
        return 0


def getMaterial(obj):
    # Temporary fix until the SetMaterials works
    # Somehow (now Feb 20) if a new gdml object is added
    # the defalut material is Geant4, and SetMaterials fails to change it
    from .GDMLMaterials import getMaterialsList
    GDMLShared.trace('get Material : '+obj.Label)
    if hasattr(obj, 'material'):
        material = obj.material
    elif hasattr(obj, 'Tool'):
        GDMLShared.trace('Has tool - check Base')
        material = getMaterial(obj.Base)
    else:
        return getDefaultMaterial()

    if material not in getMaterialsList():
        material = getDefaultMaterial()

    return material


'''
def printObjectInfo(xmlVol, volName, xmlParent, parentName):
    print("Process Object : "+obj.Label+' Type '+obj.TypeId)
    if xmlVol is not None :
       xmlstr = ET.tostring(xmlVol) 
    else :
       xmlstr = 'None'
    print('Volume : '+volName+' : '+str(xmlstr))
    if xmlParent is not None :
       xmlstr = ET.tostring(xmlParent) 
    else :
       xmlstr = 'None'
    print('Parent : '+str(parentName)+' : '+str(xmlstr))
'''


def exportCone(name, radius, height):
    cylEl = ET.SubElement(solids, 'cone', {'name': name,
                                           'rmin1': '0',
                                           'rmax1': str(radius),
                                           'rmin2': '0',
                                           'rmax2': str(radius),
                                           'z': str(height),
                                           'startphi': '0',
                                           'deltaphi': '360',
                                           'aunit': 'deg', 'lunit': 'mm'})
    return cylEl


def insertXMLvolume(name):
    # Insert at beginning for sub volumes
    GDMLShared.trace('insert xml volume : ' + name)
    elem = ET.Element('volume', {'name': name})
    global structure
    structure.insert(0, elem)
    return elem


def insertXMLvolObj(obj):
    # name = cleanVolName(obj, obj.Label)
    name = obj.Label
    return insertXMLvolume(name)


def insertXMLassembly(name):
    # Insert at beginning for sub volumes
    GDMLShared.trace('insert xml assembly : ' + name)
    elem = ET.Element('assembly', {'name': name})
    global structure
    structure.insert(0, elem)
    return elem


def insertXMLassemObj(obj):
    # name = cleanVolName(obj, obj.Label)
    name = obj.Label
    return insertXMLassembly(name)


def createXMLvol(name):
    return ET.SubElement(structure, 'volume', {'name': name})


def processAssembly(vol, xmlVol, xmlParent, parentName):
    # vol - Volume Object
    # xmlVol - xml of this assembly
    # xmlParent - xml of this volumes Paretnt
    # App::Part will have Booleans & Multifuse objects also in the list
    # So for s in list is not so good
    # type 1 straight GDML type = 2 for GEMC
    # xmlVol could be created dummy volume
    GDMLShared.setTrace(True)
    volName = vol.Label
    # volName = cleanVolName(vol, vol.Label)
    GDMLShared.trace('Process Assembly : '+volName)
    # if GDMLShared.getTrace() == True :
    #   printVolumeInfo(vol, xmlVol, xmlParent, parentName)
    assemObjs = assemblyHeads(vol)
    print(f'processAssembly: vol.TypeId {vol.TypeId}')
    for obj in assemObjs:
        if obj.TypeId == 'App::Part':
            processVolAssem(obj, xmlVol, volName)
        elif obj.TypeId == 'App::Link':
            print('Process Link')
            # objName = cleanVolName(obj, obj.Label)
            addPhysVolPlacement(obj, xmlVol, obj.LinkedObject.Label,
                                obj.Placement)
        else:
            _ = processVolume(obj, xmlVol)

    addPhysVolPlacement(vol, xmlParent, volName, vol.Placement)


def processVolume(vol, xmlParent, volName=None):
    from .solidsExporter import SolidExporter

    # vol - Volume Object
    # xmlVol - xml of this volume
    # xmlParent - xml of this volumes Paretnt
    # App::Part will have Booleans & Multifuse objects also in the list
    # So for s in list is not so good
    # type 1 straight GDML type = 2 for GEMC
    # xmlVol could be created dummy volume
    if vol.TypeId == 'App::Link':
        print('Volume is Link')
        # objName = cleanVolName(obj, obj.Label)
        addPhysVolPlacement(vol, xmlParent, vol.LinkedObject.Label,
                            vol.Placement)
        return

    if volName is None:
        volName = vol.Label
    if vol.TypeId == 'App::Part':
        topObject = topObj(vol)
    else:
        topObject = vol
    if topObject is None:
        return

    if isMultiPlacement(topObject):
        xmlVol, volName = processMultiPlacement(topObject)
        partPlacement = topObject.Placement

    else:
        solidExporter = SolidExporter.getExporter(topObject)
        if solidExporter is None:
            return
        solidExporter.export()
        print(f'solids count {len(list(solids))}')
        # 1- adds a <volume element to <structure with name volName
        xmlVol = insertXMLvolume(volName)
        # 2- add material info to the generated <volume pointerd to by xmlVol
        addVolRef(xmlVol, volName, topObject, solidExporter.name())
        # 3- add a <physvol. A <physvol, can go under the <worlVol, or under
        #    a <assembly
        # first we need to convolve the solids placement, with the vols placement
        partPlacement = solidExporter.placement()

    volPlace = vol.Placement

    # if the container is NOT an assembly, we convolve the placement of
    # the solid with that of its container. An assembly gets its own
    # placement in thr GDML, so don't do the convolution here
    if xmlParent is not None:
        if xmlParent.tag == "assembly":
            print("xmplParent tag: "+xmlParent.tag)
            placement = partPlacement
        else:
            placement = volPlace*partPlacement
    else:
        placement = partPlacement  # addphysVolPlacement does not add a <physvol
                                # if xmlParent is None, so this has no effect
    addPhysVolPlacement(vol, xmlParent, volName, placement)
    if hasattr(vol, 'SensDet'):
        if vol.SensDet is not None:
            print('Volume : ' + volName)
            print('SensDet : ' + vol.SensDet)
            ET.SubElement(xmlVol, 'auxiliary', {'auxtype': 'SensDet',
                                                'auxvalue': vol.SensDet})
    print(f'Processed Volume : {volName}')

    return xmlVol


def processVolAssem(vol, xmlParent, parentName):
    # vol - Volume Object
    # xmlVol - xml of this volume
    # xmlParent - xml of this volumes Paretnt
    # xmlVol could be created dummy volume
    print('process volasm '+vol.Label)
    volName = vol.Label
    if isAssembly(vol):
        newXmlVol = insertXMLassembly(volName)
        processAssembly(vol, newXmlVol, xmlParent, parentName)
    else:
        processVolume(vol, xmlParent)


def printVolumeInfo(vol, xmlVol, xmlParent, parentName):
    if xmlVol is not None:
        xmlstr = ET.tostring(xmlVol)
    else:
        xmlstr = 'None'
    print(xmlstr)
    GDMLShared.trace('     '+vol.Label + ' - ' + str(xmlstr))
    if xmlParent is not None:
        xmlstr = ET.tostring(xmlParent)
    else:
        xmlstr = 'None'
    GDMLShared.trace('     Parent : ' + str(parentName) + ' : ' + str(xmlstr))


def processMultiPlacement(obj):
    from .solidsExporter import SolidExporter

    print(f'procesMultiPlacement {obj.Label}')

    def getChildren(obj):
        children = []
        for o in obj.OutList:
            if o.TypeId != 'App::Origin':
                children.append(o)
                children += getChildren(o)

        return children

    children = [obj] + getChildren(obj)
    # export first solid in solids (booleans etc)
    for i, s in enumerate(children):
        if SolidExporter.isSolid(s):
            exporter = SolidExporter.getExporter(s)
            exporter.export()
            solidName = exporter.name()
            volName = 'LV-'+solidName
            volXML = insertXMLvolume(volName)
            addVolRef(volXML, obj.Label, s, solidName)
            break
    placers = children[:i]  # placers without the solids
    for pl in reversed(placers):
        placer = MultiPlacer.getPlacer(pl)
        placer.place(volName)
        volName = placer.name()
        volXML = placer.xml()

    return volXML, volName  # name of last placer (an assembly)


def createWorldVol(volName):
    print("Need to create Dummy Volume and World Box ")
    bbox = FreeCAD.BoundBox()
    boxName = defineWorldBox(bbox)
    worldVol = ET.SubElement(structure, 'volume', {'name': volName})
    ET.SubElement(worldVol, 'materialref', {'ref': 'G4_AIR'})
    ET.SubElement(worldVol, 'solidref', {'ref': boxName})
    ET.SubElement(gxml, 'volume', {'name': volName, 'material': 'G4_AIR'})
    return worldVol


def isAssembly(obj):
    # return True if obj is an assembly
    # to be an assembly the obj must be:
    # (1) and App::Part or an App::Link and
    # (2) it has either (1) At least one App::Part as a subpart or
    #                   (2) more than one "terminal" object
    # A terminal object is one that has associated with it ONE volume
    # A volume refers to ONE solid
    # A terminal item CAN be a boolean, or an extrusion (and in the future
    # a chamfer or a fillet. So a terminal element need NOT have an empty
    # OutList
    # N.B. App::Link is treated as a non-assembly, eventhough it might be linked
    # to an assembly, because all we need to place it is the volref of its link

    subObjs = []
    if obj.TypeId != 'App::Part':
        return False
    for ob in obj.OutList:
        if ob.TypeId == 'App::Part' or ob.TypeId == 'App::Link':
            return True  # Yes, even if ONE App::Part is under this, we treat it as an assembly
        else:
            if ob.TypeId != 'App::Origin':
                subObjs.append(ob)

    # now remove any OutList objects from the subObjs
    for subObj in subObjs[:]:  # the slice is a COPY of the list, not the list itself
        if hasattr(subObj, 'OutList'):
            for o in subObj.OutList:
                if o in subObjs:
                    subObjs.remove(o)

    if len(subObjs) > 1:
        return True
    else:
        return False


def assemblyHeads(obj):
    # return a list of subassembly heads for this object
    # Subassembly heads are themselves either assemblies
    # or terminal objects (those that produce a <volume and <physvol)
    assemblyHeads = []
    if isAssembly(obj):
        for ob in obj.OutList:
            if ob.TypeId == 'App::Part' or ob.TypeId == 'App::Link':
                print(f'adding {ob.Label}')
                assemblyHeads.append(ob)
            else:
                if ob.TypeId != 'App::Origin':
                    print(f'adding {ob.Label}')
                    assemblyHeads.append(ob)

        # now remove any OutList objects from from the subObjs
        for subObj in assemblyHeads[:]:  # the slice is a COPY of the list, not the list itself
            if hasattr(subObj, 'OutList'):
                # App:Links has the object they are linked to in their OutList
                # We do not want to remove the link!
                if subObj.TypeId == 'App::Link':
                    print(f'skipping {subObj.Label}')
                    continue
                for o in subObj.OutList:
                    if o in assemblyHeads:
                        print(f'removing {ob.Label}')
                        assemblyHeads.remove(o)

    return assemblyHeads


def topObj(obj):
    # The topmost object in an App::Part
    # The App::Part is assumed NOT to be an assembly
    if isAssembly(obj):
        print(f'***** Non Assembly expected  for {obj.Label}')
        return

    if not hasattr(obj, 'OutList'):
        return obj

    if len(obj.OutList) == 0:
        return obj

    sublist = []
    for ob in obj.OutList:
        if ob.TypeId != 'App::Origin':
            sublist.append(ob)

    for subObj in sublist[:]:
        if hasattr(subObj, 'OutList'):
            for o in subObj.OutList:
                if o in sublist:
                    sublist.remove(o)

    if len(sublist) > 1:
        print(f'Found more than one top object in {obj.Label}. \n Returning first only')
    elif len(sublist) == 0:
        return None

    return sublist[0]


def isMultiPlacement(obj):
    return obj.TypeId == 'Part::Mirroring'


def countGDMLObj(objList):
    # Return counts GDML objects exportables
    # #rint('countGDMLObj')
    GDMLShared.trace('countGDMLObj')
    gcount = 0
    # print(range(len(objList)))
    for idx in range(len(objList)):
        # print('idx : '+str(idx))
        obj = objList[idx]
        if obj.TypeId == 'Part::FeaturePython':
            gcount += 1
        if obj.TypeId == 'Part::Cut' \
           or obj.TypeId == 'Part::Fuse' \
           or obj.TypeId == 'Part::Common':
            gcount -= 1
    # print('countGDMLObj - Count : '+str(gcount))
    GDMLShared.trace('countGDMLObj - gdml : ' + str(gcount))
    return gcount


def checkGDMLstructure(objList):
    # Should be
    # World Vol - App::Part
    # App::Origin
    # GDML Object
    GDMLShared.trace('check GDML structure')
    GDMLShared.trace(objList)
    # print(objList)
    cnt = countGDMLObj(objList)
    if cnt > 1:  # More than one GDML Object need to insert Dummy
        return False
    if cnt == 1 and len(objList) == 2:  # Just a single GDML obj insert Dummy
        return False
    return True


def locateXMLvol(vol):
    global structure
    xmlVol = structure.find("volume[@name='%s']" % vol.Label)
    return xmlVol


def exportWorldVol(vol, fileExt):

    global WorldVOL
    WorldVOL = vol.Label
    if fileExt != '.xml':
        print('Export World Process Volume : ' + vol.Label)
        GDMLShared.trace('Export Word Process Volume' + vol.Label)
        ET.SubElement(setup, 'world', {'ref': vol.Label})

        if checkGDMLstructure(vol.OutList) is False:
            GDMLShared.trace('Insert Dummy Volume')
            createXMLvol('dummy')
            xmlParent = createWorldVol(vol.Label)
            parentName = vol.Label
            addPhysVol(xmlParent, 'dummy')
        else:
            GDMLShared.trace('Valid Structure')
            xmlParent = None
            parentName = None
    else:
        xmlParent = None
        parentName = None

    if hasattr(vol, 'OutList'):
        # print(vol.OutList)
        cnt = countGDMLObj(vol.OutList)
    print('Root GDML Count ' + str(cnt))

    if cnt > 0:  # one GDML defining world volume
        if isAssembly(vol):
            heads = assemblyHeads(vol)
            worlSolid = heads[0]
            xmlVol = processVolume(worlSolid, xmlParent, volName=WorldVOL)
            for obj in heads[1:]:  # skip first volume (done above)
                processVolAssem(obj, xmlVol, WorldVOL)
        else:
            xmlVol = processVolume(vol, xmlParent)
    else:  # no volume defining world
        xmlVol = insertXMLassembly(vol.Label)
        processAssembly(vol, xmlVol, xmlParent, parentName)


def exportElementAsXML(dirPath, fileName, flag, elemName, elem):
    # gdml is a global
    global gdml, docString, importStr
    if elem is not None:
        # xmlElem = ET.Element('xml')
        # xmlElem.append(elem)
        # indent(xmlElem)
        if flag is True:
            filename = fileName+'-' + elemName + '.xml'
        else:
            filename = elemName + '.xml'
        # ET.ElementTree(xmlElem).write(os.path.join(dirPath,filename))
        ET.ElementTree(elem).write(os.path.join(dirPath, filename))
        docString += '<!ENTITY ' + elemName + ' SYSTEM "' + filename+'">\n'
        gdml.append(ET.Entity(elemName))


def exportGDMLstructure(dirPath, fileName):
    global gdml, docString, importStr
    print("Write GDML structure to Directory")
    gdml = initGDML()
    docString = '\n<!DOCTYPE gdml [\n'
    # exportElementAsXML(dirPath, fileName, False, 'constants',constants)
    exportElementAsXML(dirPath, fileName, False, 'define', define)
    exportElementAsXML(dirPath, fileName, False, 'materials', materials)
    exportElementAsXML(dirPath, fileName, True, 'solids', solids)
    exportElementAsXML(dirPath, fileName, True, 'structure', structure)
    exportElementAsXML(dirPath, fileName, False, 'setup', setup)
    print(f"setup : {setup}")
    docString += ']>\n'
    # print(docString)
    # print(len(docString))
    # gdml = ET.fromstring(docString.encode("UTF-8"))
    indent(gdml)
    ET.ElementTree(gdml).write(os.path.join(dirPath, fileName+'.gdml'),
                               doctype=docString.encode('UTF-8'))
    print("GDML file structure written")


def exportGDML(first, filepath, fileExt):
    from . import GDMLShared
    global zOrder

    # GDMLShared.setTrace(True)
    GDMLShared.trace('exportGDML')
    print("====> Start GDML Export 1.6")
    print('File extension : '+fileExt)

    GDMLstructure()
    zOrder = 1
    processMaterials()
    exportWorldVol(first, fileExt)
    # format & write GDML file
    # xmlstr = ET.tostring(structure)
    # print('Structure : '+str(xmlstr))
    if fileExt == '.gdml':
        # indent(gdml)
        print(len(list(solids)))
        print("Write to gdml file")
        # ET.ElementTree(gdml).write(filepath, 'utf-8', True)
        # ET.ElementTree(gdml).write(filepath, xml_declaration=True)
        ET.ElementTree(gdml).write(filepath, pretty_print=True,
                                   xml_declaration=True)
        print("GDML file written")

    if fileExt == '.GDML':
        filePath = os.path.split(filepath)
        print('Input File Path : '+filepath)
        fileName = os.path.splitext(filePath[1])[0]
        print('File Name : '+fileName)
        dirPath = os.path.join(filePath[0], fileName)
        print('Directory Path : '+dirPath)
        if os.path.exists(dirPath) is False:
            if os.path.isdir(dirPath) is False:
                os.makedirs(dirPath)
        if os.path.isdir(dirPath) is True:
            exportGDMLstructure(dirPath, fileName)
        else:
            print('Invalid Path')
            # change to Qt Warning

    if fileExt == '.xml':
        xmlElem = ET.Element('xml')
        xmlElem.append(solids)
        xmlElem.append(structure)
        indent(xmlElem)
        ET.ElementTree(xmlElem).write(filepath)
        print("XML file written")


def exportGDMLworld(first, filepath, fileExt):
    if filepath.lower().endswith('.gdml'):
        # GDML Export
        print('GDML Export')
        # if hasattr(first,'InList') :
        #   print(len(first.InList))

        if hasattr(first, 'OutList'):
            cnt = countGDMLObj(first.OutList)
            GDMLShared.trace('Count : ' + str(cnt))
            if cnt > 1:
                from .GDMLQtDialogs import showInvalidWorldVol
                showInvalidWorldVol()
            else:
                exportGDML(first, filepath, fileExt)


def hexInt(f):
    return hex(int(f*255))[2:].zfill(2)


def formatPosition(pos):
    s = str(pos[0]) + '*mm ' + str(pos[1]) + '*mm ' +str(pos[2]) + '*mm'
    print(s)
    return s


def scanForStl(first, gxml, path, flag):

    from .GDMLColourMap import lookupColour

    # if flag == True ignore Parts that convert
    print('scanForStl')
    print(first.Name+' : '+first.Label+' : '+first.TypeId)
    while switch(first.TypeId):

        if case("App::Origin"):
            # print("App Origin")
            return

        if case("App::GeoFeature"):
            # print("App GeoFeature")
            return

        if case("App::Line"):
            # print("App Line")
            return

        if case("App::Plane"):
            # print("App Plane")
            return

        break

    if flag is True:
        #
        #  Now deal with objects that map to GDML solids
        #
        while switch(first.TypeId):
            if case("Part::FeaturePython"):
                return

            if case("Part::Box"):
                print("    Box")
                return

            if case("Part::Cylinder"):
                print("    Cylinder")
                return

            if case("Part::Cone"):
                print("    Cone")
                return

            if case("Part::Sphere"):
                print("    Sphere")
                return

            break

    # Deal with Booleans which will have Tool
    if hasattr(first, 'Tool'):
        print(first.TypeId)
        scanForStl(first.Base, gxml, path, flag)
        scanForStl(first.Tool, gxml, path, flag)

    if hasattr(first, 'OutList'):
        for obj in first.OutList:
            scanForStl(obj, gxml, path, flag)

    if first.TypeId != 'App::Part':
        if hasattr(first, 'Shape'):
            print('Write out stl')
            print('===> Name : '+first.Name+' Label : '+first.Label+' \
            Type :'+first.TypeId+' : '+str(hasattr(first, 'Shape')))
            newpath = os.path.join(path, first.Label + '.stl')
            print('Exporting : ' + newpath)
            first.Shape.exportStl(newpath)
            # Set Defaults
            colHex = 'ff0000'
            mat = 'G4Si'
            if hasattr(first.ViewObject, 'ShapeColor'):
                # print(dir(first))
                col = first.ViewObject.ShapeColor
                colHex = hexInt(col[0]) + hexInt(col[1]) + hexInt(col[2])
                print('===> Colour '+str(col) + ' '+colHex)
                mat = lookupColour(col)
                print('Material : '+mat)
                if hasattr(first, 'Placement'):
                    print(first.Placement.Base)
                    pos = formatPosition(first.Placement.Base)
                    ET.SubElement(gxml, 'volume', {'name': first.Label,
                                                   'color': colHex,
                                                   'material': mat,
                                                   'position': pos})


def exportGXML(first, path, flag):
    print('Path : '+path)
    # basename = 'target_'+os.path.basename(path)
    gxml = ET.Element('gxml')
    print('ScanForStl')
    scanForStl(first, gxml, path, flag)
    # format & write gxml file
    indent(gxml)
    print("Write to gxml file")
    # ET.ElementTree(gxml).write(os.path.join(path,basename+'.gxml'))
    ET.ElementTree(gxml).write(os.path.join(path, 'target_cad.gxml'))
    print("gxml file written")


def exportMaterials(first, filename):
    if filename.lower().endswith('.xml'):
        print('Export Materials to XML file : '+filename)
        xml = ET.Element('xml')
        global define
        define = ET.SubElement(xml, 'define')
        global materials
        materials = ET.SubElement(xml, 'materials')
        processMaterials()
        indent(xml)
        ET.ElementTree(xml).write(filename)
    else:
        print('File extension must be xml')


def create_gcard(path, flag):
    basename = os.path.basename(path)
    print('Create gcard : '+basename)
    print('Path : '+path)
    gcard = ET.Element('gcard')
    ET.SubElement(gcard, 'detector', {'name': 'target_cad', 'factory': 'CAD'})
    if flag is True:
        ET.SubElement(gcard, 'detector', {
           'name': 'target_gdml', 'factory': 'GDML'})
    indent(gcard)
    path = os.path.join(path, basename + '.gcard')
    ET.ElementTree(gcard).write(path)


def checkDirectory(path):
    if not os.path.exists(path):
        print('Creating Directory : ' + path)
        os.mkdir(path)


def exportGEMC(first, path, flag):
    # flag = True  GEMC - GDML
    # flag = False just CAD
    global gxml

    print('Export GEMC')
    # basename = os.path.basename(path)
    print(path)
    print(flag)
    checkDirectory(path)
    # Create CAD directory
    cadPath = os.path.join(path, 'cad')
    checkDirectory(cadPath)
    # Create gcard
    create_gcard(path, flag)
    exportGXML(first, cadPath, flag)
    if flag is True:
        print('Create GDML directory')
        gdmlPath = os.path.join(path, 'gdml')
        checkDirectory(gdmlPath)
        # gdmlFilePath  = os.path.join(gdmlPath,basename+'.gdml')
        gdmlFilePath = os.path.join(gdmlPath, 'target_gdml.gdml')
        exportGDML(first, gdmlFilePath, 'gdml')
        # newpath = os.path.join(gdmlPath,basename+'.gxml')
        newpath = os.path.join(gdmlPath, 'target_gdml.gxml')
        indent(gxml)
        ET.ElementTree(gxml).write(newpath)


def export(exportList, filepath):
    "called when FreeCAD exports a file"
    global refPlacement

    refPlacement = {}  # a dictionary of name as key, and placement as value
                       # the name could that of <solid, an <assembly, or a <physvol

    first = exportList[0]
    print(f'Export Volume: {first.Label}')

    import os
    path, fileExt = os.path.splitext(filepath)
    print('filepath : '+path)
    print('file extension : '+fileExt)

    if fileExt.lower() == '.gdml':
        if first.TypeId == "App::Part":
            exportGDMLworld(first, filepath, fileExt)

        elif first.Label == "Materials":
            exportMaterials(first, filepath)

        else:
            print("Needs to be a Part for export")
            from PySide import QtGui
            QtGui.QMessageBox.critical(None,
                                       'Need to select a Part for export',
                                       'Press OK')

    elif fileExt.lower == '.xml':
        print('Export XML structure & solids')
        exportGDML(first, filepath, '.xml')

    if fileExt == '.gemc':
        exportGEMC(first, path, False)

    elif fileExt == '.GEMC':
        exportGEMC(first, path, True)
