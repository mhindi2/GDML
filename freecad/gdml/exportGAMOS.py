#**************************************************************************
#*                                                                        * 
#*   Copyright (c) 2019 Keith Sloan <keith@sloan-home.co.uk>              *
#*             (c) 2020 Dam Lambert                                       *
#*                                                                        *
#*   This program is free software; you can redistribute it and/or modify *
#*   it under the terms of the GNU Lesser General Public License (LGPL)   *
#*   as published by the Free Software Foundation; either version 2 of    *
#*   the License, or (at your option) any later version.                  *
#*   for detail see the LICENCE text file.                                *
#*                                                                        *
#*   This program is distributed in the hope that it will be useful,      *
#*   but WITHOUT ANY WARRANTY; without even the implied warranty of       *
#*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
#*   GNU Library General Public License for more details.                 *
#*                                                                        *
#*   You should have received a copy of the GNU Library General Public    *
#*   License along with this program; if not, write to the Free Software  *
#*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 *
#*   USA                                                                  *
#*                                                                        *
#*   Acknowledgements : Ideas & code copied from                          *
#*                      https://github.com/ignamv/geanTipi                *
#*                                                                        *
#***************************************************************************
__title__="FreeCAD - GAMOS exporter Version"
__author__ = "Keith Sloan <keith@sloan-home.co.uk>"
__url__ = ["https://github.com/KeithSloan/FreeCAD_Geant4"]

if open.__module__ == '__builtin__':
    pythonopen = open # to distinguish python built-in open function from the one declared here

import FreeCAD, os, Part, math
from FreeCAD import Vector
from .GDMLObjects import GDMLcommon, GDMLBox, GDMLTube

# modif add 
#from .GDMLObjects import getMult, convertionlisteCharToLunit

import sys

try: import FreeCADGui
except ValueError: gui = False
else: gui = True

global zOrder

from .GDMLObjects import GDMLQuadrangular, GDMLTriangular, \
                        GDML2dVertex, GDMLSection, \
                        GDMLmaterial, GDMLfraction, \
                        GDMLcomposite, GDMLisotope, \
                        GDMLelement, GDMLconstant

from . import GDMLShared

#***************************************************************************
# Tailor following to your requirements ( Should all be strings )          *
# no doubt there will be a problem when they do implement Value
if open.__module__ in ['__builtin__', 'io']:
    pythonopen = open # to distinguish python built-in open function from the one declared here

def verifNameUnique(name):
   # need to be done!!
   return True

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

def nameFromLabel(label) :
    if ' ' not in label :
       return label
    else :
       return(label.split(' ')[0])

def exportDefine(name, v) :
    global define
    print('define : '+name)
    #print(v)
    #print(v[0])
    #ET.SubElement(define,'position',{'name' : name, 'unit': 'mm',   \
    #            'x': str(v[0]), 'y': str(v[1]), 'z': str(v[2]) })

def defineWorldBox(bbox,fp):
    for obj in FreeCAD.ActiveDocument.Objects :
        # print("{} + {} = ".format(bbox, obj.Shape.BoundBox))
        if hasattr(obj,"Shape"):
           bbox.add(obj.Shape.BoundBox)
        if hasattr(obj,"Mesh"):
           bbox.add(obj.Mesh.BoundBox)
        if hasattr(obj,"Points"):
           bbox.add(obj.Points.BoundBox)
    #   print(bbox)
    # Solids get added to solids section of gdml ( solids is a global )
    name = 'WorldBox'
    fp.write(':SOLID '+name+' BOX '+'1000'+'1000'+'1000')
    #ET.SubElement(solids, 'box', {'name': name,
    #                         'x': str(1000), \
    #                         'y': str(1000), \
    #                         'z': str(1000), \
    #                 #'x': str(2*max(abs(bbox.XMin), abs(bbox.XMax))), \
    #                 #'y': str(2*max(abs(bbox.YMin), abs(bbox.YMax))), \
    #                 #'z': str(2*max(abs(bbox.ZMin), abs(bbox.ZMax))), \
    #                 'lunit': 'mm'})
    return(name)

def processPlanar(obj, shape, name ) :
    print ('Polyhedron ????')
    global defineCnt
    #
    #print("Add tessellated Solid")
    #tess = ET.SubElement(solids,'tessellated',{'name': name})
    #print("Add Vertex positions")
    for f in shape.Faces :
       baseVrt = defineCnt
       for vrt in f.Vertexes :
           vnum = 'v'+str(defineCnt)
           #ET.SubElement(define, 'position', {'name': vnum, \
           #   'x': str(vrt.Point.x), \
           #   'y': str(vrt.Point.y), \
           #   'z': str(vrt.Point.z), \
           #   'unit': 'mm'})
           #defineCnt += 1
       #print("Add vertex to tessellated Solid")
       vrt1 = 'v'+str(baseVrt)
       vrt2 = 'v'+str(baseVrt+1)
       vrt3 = 'v'+str(baseVrt+2)
       vrt4 = 'v'+str(baseVrt+3)
       NumVrt = len(f.Vertexes)
       if NumVrt == 3 :
          pass
          #ET.SubElement(tess,'triangular',{ \
          #            'vertex1': vrt1, \
          #            'vertex2': vrt2, \
          #            'vertex3': vrt3, \
          #            'type': 'ABSOLUTE'})
       elif NumVrt == 4 :
          pass
          #ET.SubElement(tess,'quadrangular',{ \
          #            'vertex1': vrt1, \
          #            'vertex2': vrt2, \
          #            'vertex3': vrt3, \
          #            'vertex4': vrt4, \
          #            'type': 'ABSOLUTE'})

def checkShapeAllPlanar(Shape) :
    for f in Shape.Faces :
        if f.Surface.isPlanar() == False :
           return False
        break
    return True

#    Add XML for TessellateSolid
def mesh2Tessellate(mesh, name) :
     global defineCnt

     baseVrt = defineCnt
     print ("mesh")
     print (mesh)
     print (dir(mesh))
     print ("Facets")
     print (mesh.Facets)
     print ("mesh topology")
     print (dir(mesh.Topology))
     print (mesh.Topology)
#
#    mesh.Topology[0] = points
#    mesh.Topology[1] = faces
#
#    First setup vertex in define section vetexs (points) 
     print("Add Vertex positions")
     for fc_points in mesh.Topology[0] : 
         print(fc_points)
         v = 'v'+str(defineCnt)
         #ET.SubElement(define, 'position', {'name': v, \
         #         'x': str(fc_points[0]), \
         #         'y': str(fc_points[1]), \
         #         'z': str(fc_points[2]), \
         #         'unit': 'mm'})
         #defineCnt += 1         
#                  
#     Add faces
#
     print("Add Triangular vertex")
     #tess = ET.SubElement(solids,'tessellated',{'name': name})
     for fc_facet in mesh.Topology[1] : 
       print(fc_facet)
       vrt1 = 'v'+str(baseVrt+fc_facet[0])
       vrt2 = 'v'+str(baseVrt+fc_facet[1])
       vrt3 = 'v'+str(baseVrt+fc_facet[2])
       #ET.SubElement(tess,'triangular',{ \
       #  'vertex1': vrt1, 'vertex2': vrt2 ,'vertex3': vrt3, 'type': 'ABSOLUTE'})


def processMesh(obj, Mesh, Name) :
    #  obj needed for Volune names
    #  object maynot have Mesh as part of Obj
    #  Name - allows control over name
    print("Create Tessellate Logical Volume")
    createLVandPV(obj, Name, 'Tessellated')
    mesh2Tessellate(Mesh, Name)
    return(Name)

def shape2Mesh(shape) :
     import MeshPart
     return (MeshPart.meshFromShape(Shape=shape, Deflection = 0.0))
#            Deflection= params.GetFloat('meshdeflection',0.0)) 

def processObjectShape(obj) :
    # Check if Planar
    # If plannar create Tessellated Solid with 3 & 4 vertex as appropriate
    # If not planar create a mesh and the a Tessellated Solid with 3 vertex
    print("Process Object Shape")
    print(obj)
    print(obj.PropertiesList)
    if not hasattr(obj,'Shape') :
       return 
    shape = obj.Shape
    #print (shape)
    #print(shape.ShapeType)
    while switch(shape.ShapeType) : 
      if case("Mesh::Feature") :
         print("Mesh - Should not occur should have been handled")
         #print("Mesh")
         #tessellate = mesh2Tessellate(mesh) 
         #return(tessellate)
         #break

         print("ShapeType Not handled")
         print(shape.ShapeType)
         break

#   Dropped through to here
#   Need to check has Shape

    print('Check if All planar')
    planar = checkShapeAllPlanar(shape)
    print(planar)

    if planar :
       return(processPlanar(obj,shape,obj.Name))

    else :
       # Create Mesh from shape & then Process Mesh
       #to create Tessellated Solid in Geant4
       return(processMesh(obj,shape2Mesh(shape),obj.Name))

def BoxObjectVariables(obj) :
    # This for non GDML Box
 

    vars = 'x' + str(obj.Length.Value) + 'y' + str(obj.Width.Value) + \
                          'z' + str(obj.Height.Value)
    delta = FreeCAD.Vector(obj.Length.Value / 2, \
                           obj.Width.Value / 2,  \
                           obj.Height.Value / 2)

    return(vars, delta)

def CylinderObjectVariables(obj) :
    vars = 'rmax' + str(obj.Radius.Value) + \
           'deltaphi' + str(float(obj.Angle)) + \
           'z' + str(obj.Height.Value)
    delta = FreeCAD.Vector(0, 0, obj.Height.Value / 2)
    return(vars, delta)

def ConeObjectVariable(obj) :
    # Needs unique Name
    coneName = obj.Name
    #ET.SubElement(solids, 'cone',{'name': coneName, \
    #                       'rmax1': str(obj.Radius1.Value),  \
    #                       'rmax2': str(obj.Radius2.Value),  \
    #                       'deltaphi': str(float(obj.Angle)), \
    #                       'aunit': obj.aunit,
    #                       'z': str(obj.Height.Value),  \
    #                       'lunit' : 'mm'})
    if addVolsFlag :
       # Adjustment for position in GDML
       delta = FreeCAD.Vector(0, 0, obj.Height.Value / 2)
       createAdjustedLVandPV(obj, obj.Name, coneName, delta)
    return(coneName)

def processSection(obj, addVolsflag) :
    print("Process Section")
    #ET.SubElement(solids, 'section',{'vertex1': obj.v1, \
    #        'vertex2': obj.v2, 'vertex3': obj.v3, 'vertex4': obj.v4, \
    #        'type': obj.vtype})


def processSphereObject(obj, addVolsFlag) :
    # Needs unique Name
    #modif lambda (if we change the name here, each time we import and export the file, the name will be change 
    #sphereName = 'Sphere' + obj.Name
    sphereName = obj.Name

    #ET.SubElement(solids, 'sphere',{'name': sphereName, \
    #                       'rmax': str(obj.Radius.Value), \
    #                       'starttheta': str(90.-float(obj.Angle2)), \
    #                       'deltatheta': str(float(obj.Angle2-obj.Angle1)), \
    #                       'deltaphi': str(float(obj.Angle3)), \
    #                       'aunit': obj.aunit,
    #                       'lunit' : 'mm'})
    if addVolsFlag :
       createLVandPV(obj,obj.Name,sphereName)
    return(sphereName)

def exportPosition(name, xml, pos) :
    global POScount
    GDMLShared.trace('export Position')
    GDMLShared.trace(pos)
    x = pos[0]
    y = pos[1]
    z = pos[2]
    posName = 'P-'+name+str(POScount)
    POScount += 1
    #posxml = ET.SubElement(define,'position',{'name' : posName, \
    #                      'unit': 'mm'})
    if x != 0 :
       posxml.attrib['x'] = str(x)
    if y != 0 :
       posxml.attrib['y'] = str(y)
    if z != 0 :
       posxml.attrib['z'] = str(z)
    #ET.SubElement(xml,'positionref',{'ref' : posName})

def exportRotation(name, xml, Rotation) :
    global ROTcount
    if Rotation.Angle != 0 :
        angles = Rotation.toEuler()
        GDMLShared.trace("Angles")
        GDMLShared.trace(angles)
        a0 = angles[0]
        a1 = angles[1]
        a2 = angles[2]
        if a0!=0 or a1!=0 or a2!=0 :
            rotName = 'R-'+name+str(ROTcount)
            ROTcount += 1
            #rotxml = ET.SubElement(define, 'rotation', {'name': rotName, \
            #        'unit': 'deg'})
            if abs(a2) != 0 :
                rotxml.attrib['x']=str(a2)
            if abs(a1) != 0 :
                rotxml.attrib['y']=str(a1)
            if abs(a0) != 0 :
                rotxml.attrib['z']=str(a0)
            #ET.SubElement(xml, 'rotationref', {'ref': rotName})

def processPosition(obj, solid) :
    if obj.Placement.Base == FreeCAD.Vector(0,0,0) :
        return
    GDMLShared.trace("Define position & references to Solid")
    exportPosition(obj.Name, solid, obj.Placement.Base)

def processRotation(obj, solid) :
    if obj.Placement.Rotation.Angle == 0 :
       return
    GDMLShared.trace('Deal with Rotation')
    exportRotation(obj.Name,solid,obj.Placement.Rotation)


def testDefaultPlacement(obj) :
    #print(dir(obj.Placement.Rotation))
    print('Test Default Placement : '+obj.Name)
    print('No longer used ??')
    print(obj.Placement.Base)
    print(obj.Placement.Rotation.Angle)
    if obj.Placement.Base == FreeCAD.Vector(0,0,0) and \
       obj.Placement.Rotation.Angle == 0 :
       return True
    else :
       return False

def processMaterialsGroupObjects(fp) :
    print("\nProcess Materials Group Objects")
    fp.write('Process Materials')
   
    for obj in FreeCAD.ActiveDocument.Objects:
       print(obj.TypeId+" : "+obj.Name)

       if obj.TypeId == "App::DocumentObjectGroupPython":
          print("   Object List : "+obj.Name)
          #print(obj)
          while switch(obj.Name) :
             if case("Constants") : 
                print("Constants")
                processConstants(obj, fp)
                fp.flush()
                break

             if case("Isotopes") :
                print("Isotopes")
                processIsotopes(obj, fp)
                fp.flush()
                break
            
             if case("Elements") :
                print("Elements")
                processElements(obj, fp)
                fp.flush()
                break

             if case("Materials") : 
                print("Materials")
                processMaterials(obj, fp)
                fp.flush()
                break
         
             break
       else :
          return

def processConstants(obj,fp):
    print('Process Constants')
    for subobj in obj.OutList :
        #if isinstance(subobj.Proxy.GDMLconstant) :
        fp.write(':P '+subobj.Name+' '+ str(subobj.Value)+'\n')

def processIsotopes(obj,fp):
    print('Process Isotopes')
    fp.write('Process IsoTopes')
    #if isinstance(obj.Proxy,GDMLisotope) :
    for subobj in obj.OutList :
        #if isinstance(subobj.Proxy,GDMLisotope) :
        print('Write ISOT')
        fp.write(':ISOT '+subobj.Name+' '+str(subobj.Z)+' '+str(subobj.N) + \
                  ' '+str(subobj.value)+'\n')

def processElement(obj,fp) :
    #if isinstance(obj.Proxy,GDMLelement) :
    #if isinstance(obj,GDMLelement) :
       if hasattr(obj,'OutList') :
          num = len(obj.OutList)
          string = ':ELEM_FROM_ISOT '+obj.Name+' '+str(num)
          for eobj in obj.OutList :
              string = string +' '+eobj.Name+' '+ str(eobj.n)
          string = string + '\n'
          fp.write(string)
       else :
          fp.write(':ELEM '+obj.Name+' '+str(obj.Z)+'\n')
 
def processElements(obj,fp):
    print('Process Elements')
    if hasattr(obj,'OutList') :
       for subobj in obj.OutList :
           processElement(subobj,fp)

def processMaterial(obj,fp) :
    #if isinstance(obj.Proxy,GDMLmaterial) :
       if hasattr(obj,'OutList') :
          num = len(obj.OutList)
          string = ':MIXT '+obj.Name+' '+str(obj.density)+str(num)
          for eobj in obj.OutList :
              string=string+' '+eobj.Name+' '+str(eobj.n)
          string = string + '\n'
          fp.write(string)
       else :
          fp.write('MATE '+obj.Name,' '+str(obj.Z) + \
                  ' '+str(obj.A)+' '+str(obj.density)+'\n')

def processMaterials(obj,fp):
    print('Process Materials')
    if hasattr(obj,'OutList') :
       for subobj in obj.OutList :
           processMaterial(subobj,fp)

def tobesorted():

          if isinstance(obj.Proxy,GDMLcomposite) :
             print("GDML Composite")
             fp.write('ELEM_FROM_ISOT '+ nameFromLabel(obj.Label))

             #ET.SubElement(item,'composite',{'n': str(obj.n), \

def processGDMLname(obj) :
    return(obj.Label.split('_',1)[1])

def appendVector(s,v) :
    return(s+str(v[0])+' '+str(v[1])+' '+str(v[2]))

def GDMLArb8ObjectVariables(obj) :
    print('Arb8 - Not Supported')

    #fp.write(*    ET.SubElement(solids, 'arb8',{'name': arb8Name, \
    #                      'v1x': str(obj.v1x),  \
    #                      'v1y': str(obj.v1y),  \
    #                      'v2x': str(obj.v2x),  \
    #                      'v2y': str(obj.v2y),  \
    #                      'v3x': str(obj.v3x),  \
    #                      'v3y': str(obj.v3y),  \
    #                      'v4x': str(obj.v4x),  \
    #                      'v4y': str(obj.v4y),  \
    #                      'v5x': str(obj.v5x),  \
    #                      'v5y': str(obj.v5y),  \
    #                      'v6x': str(obj.v6x),  \
    #                      'v6y': str(obj.v6y),  \
    #                      'v7x': str(obj.v7x),  \
    #                      'v7y': str(obj.v7y),  \
    #                      'v8x': str(obj.v8x),  \
    #                      'v8y': str(obj.v8y),  \
    #                      'dz': str(obj.dz),  \
    #                      'lunit' : obj.lunit})
    return (arb8Name)

def GDMLBoxObjectVariables(obj) :

    vars = 'BOX: '+name+' ' + str(obj.x) +  \
                        ' ' + str(obj.y) +  \
                        ' ' + str(obj.z)
    return vars

def GDMLConeObjectVariables(obj) :
   
    # add check 360 delta for shorter output 
    vars = 'CONE: '+name+ \
                          ' '+ str(obj.rmin1),  \
                          ' '+ str(obj.rmin2),  \
                          ' '+ str(obj.rmax1),  \
                          ' '+ str(obj.rmax2),  \
                          ' '+ str(obj.z),  \
                          ' '+ str(obj.startphi), \
                          ' '+ str(obj.deltaphi)
    return vars

def GDMLCutTubeObjectVariables(obj) :
    print('Cut Tube')
    #fp.write('    ET.SubElement(solids, 'cutTube',{'name': cTubeName, \
    #                      'rmin': str(obj.rmin),  \
    #                      'rmax': str(obj.rmax),  \
    #                      'startphi': str(obj.startphi), \
    #                      'deltaphi': str(obj.deltaphi), \
    #                      'aunit': obj.aunit, \
    #                      'z': str(obj.z),  \
    #                      'highX':str(obj.highX), \
    #                      'highY':str(obj.highY), \
    #                      'highZ':str(obj.highZ), \
    #                      'lowX':str(obj.lowX), \
    #                      'lowY':str(obj.lowY), \
    #                      'lowZ':str(obj.lowZ), \
    #                      'lunit' : obj.lunit})
    return(name)

def GDMLElConeObjectVariables(obj) :
    vars = 'ELLIPTICALCONE: '+name+ \
                ' ' + str(obj.dx), \
                ' ' + str(obj.dy), \
                ' ' + str(obj.zcut), \
                ' ' + str(obj.zmax)
    return vars

def GDMLEllipsoidObjectVariables(obj) :
    vars = 'ELLIPSOID: '+name+ \
                          ' ' + str(obj.ax),  \
                          ' ' + str(obj.by),  \
                          ' ' + str(obj.cz),  \
                          ' ' + str(obj.zcut1),  \
                          ' ' + str(obj.zcut2)
    return vars

def GDMLElTubeObjectVariables(obj) :
    vars = 'ELLIPTICALTUBE: ' \
                          ' ' + str(obj.dx),  \
                          ' ' + str(obj.dy),  \
                          ' ' + str(obj.dz)
    return vars

def GDMLOrbObjectVariables(obj) :
    vars = 'ORB: '+name, \
                   ' ' + str(obj.r)
    return vars

def GDMLParaObjectVariables(obj) :
    vars = 'PARA: '+name+ \
                          ' ' + str(obj.x),  \
                          ' ' + str(obj.y),  \
                          ' ' + str(obj.z),  \
                          ' ' + str(obj.alpha), \
                          ' ' + str(obj.theta), \
                          ' ' + str(obj.phi)
    return vars


def GDMLPolyconeObjectVariables(obj) :
    vars = 'POLYCONE: '+name+ \
                          ' ' + str(obj.startphi),  \
                          ' ' + str(obj.deltaphi),  \
                          ' ' + str(len(obj.OutList))
    print('Add position, Tangent distances')
    for zplane in obj.OutList :
        vars = vars + str(zplane.rmin) + str(zplane.rmax) + str(zplane.z)
    
    return vars

def GDMLSphereObjectVariables(obj) :
    vars = 'SPHERE: '+name+ \
             ' ' + str(obj.rmin)+ \
             ' ' + str(obj.rmax)+ \
             ' ' + str(obj.startphi)+ \
             ' ' + str(obj.deltaphi)+ \
             ' ' + str(obj.deltatheta)
    return vars

def GDMLTessellatedObjectVariables(obj) :
    vertex = obj.Proxy.Vertex
    facets = obj.Proxy.Facets
    vars = 'TESSElATED: '+name+' '+str(len(obj.OutList))+' '
    print(len(obj.OutList))
    for f in facets :
        vars = vars + len(f)
        vars = vars + vertex[f[0]] +' '+vertex[f[1]]+' '+vertex[f[2]]
        if len(f) == 4 :
           vars = vars + vertex[f[3]]
    return vars

def GDMLTetraObjectVariables(obj) :
    vars = 'TET: '+name
    appendVector(vars,obj,v1)
    appendVector(vars,obj,v2)
    appendVector(vars,obj,v3)
    return vars

def GDMLTorusObjectVariables(obj) :
    vars = 'TORUS: '+name+ \
             ' ' + str(obj.rmin)+ \
             ' ' + str(obj.rmax)+ \
             ' ' + str(obj.rtor)+ \
             ' ' + str(obj.startphi)+ \
             ' ' + str(obj.deltaphi)
    return vars


def GDMLTrapObjectVariables(obj) :
    vars = 'TRAP: '+name+ \
             ' ' + str(obj.z)+  \
             ' ' + str(obj.theta)+  \
             ' ' + str(obj.phi)+ \
             ' ' + str(obj.x1)+  \
             ' ' + str(obj.x2)+  \
             ' ' + str(obj.x3)+  \
             ' ' + str(obj.x4)+  \
             ' ' + str(obj.y1)+  \
             ' ' + str(obj.y2)+  \
             ' ' + str(obj.alpha)+ \
             ' ' + str(obj.alpha)
    return vars

def GDMLTrdObjectVariables(obj) :
    vars = 'TRD: '+name+ \
             ' ' + str(obj.x1)+  \
             ' ' + str(obj.x2)+  \
             ' ' + str(obj.y1)+  \
             ' ' + str(obj.y2)+  \
             ' ' + str(obj.z)
    return vars

def GDMLTubeObjectVariables(obj) :
    vars = 'TUBS: '+name+ \
             ' ' + str(obj.rmin)+  \
             ' ' + str(obj.rmax)+  \
             ' ' + str(obj.z)+  \
             ' ' + str(obj.startphi)+ \
             ' ' + str(obj.deltaphi)
    return vars

def GDMLXtruObjectVariables(obj) :
    # Needs unique Name
    xtruName = obj.Label.split('_',1)[1]

    if flag == True :
        #xtru = ET.SubElement(solids, 'xtru',{'name': xtruName, \
        #                  'lunit' : obj.lunit})
        for items in obj.OutList :
            #if items.Type == 'twoDimVertex' :
            #    #ET.SubElement(xtru, 'twoDimVertex',{'x': str(items.x), \
            #    #                   'y': str(items.y)})
            #if items.Type == 'section' :
            #    ET.SubElement(xtru, 'section',{'zOrder': str(items.zOrder), \
            #                      'zPosition': str(items.zPosition), \
            #                      'xOffset' : str(items.xOffset), \
            #                      'yOffset' : str(items.yOffset), \
            #                      'scalingFactor' : str(items.scalingFactor)})
            pass
    return(xtruName)

def GDML2dVertexVariables(obj) :
    if flag == True :
        print("Process 2d Vertex")
        #ET.SubElement(solids, 'twoDimVertex',{'x': obj.x, 'y': obj.y})

def GDMLSolidVariables(obj) :
    # Deal with GDML Solids first
    # Deal with FC Objects that convert
    #print(obj.Proxy.Type)
    name = processGDMLname(obj)
    print(fp)
    while switch(obj.Proxy.Type) :
       if case("GDMLArb8") :
          print("      GDMLArb8") 
          return(GDMLArb8ObjectVariables(obj))
          break

       if case("GDMLBox") :
          #print("      GDMLBox") 
          return(GDMLBoxObjectVariables(obj))
          break

       if case("GDMLCone") :
          #print("      GDMLCone") 
          return(GDMLConeObjectVariables(obj))
          break

       if case("GDMLcutTube") :
          #print("      GDMLcutTube") 
          return(GDMLCutTubeObjectVariables(obj))
          break
       
       if case("GDMLElCone") :
          #print("      GDMLElCone") 
          return(GDMLElConeObjectVariables(obj))
          break

       if case("GDMLEllipsoid") :
          #print("      GDMLEllipsoid") 
          return(GDMLEllipsoidObjectVariables(obj))
          break

       if case("GDMLElTube") :
          #print("      GDMLElTube") 
          return(GDMLElTubeObjectVariables(obj))
          break

       if case("GDMLOrb") :
          #print("      GDMLOrb") 
          return(GDMLOrbObjectVariables(obj))
          break

       if case("GDMLPara") :
          #print("      GDMLPara") 
          return(GDMLParaObjectVariables(obj))
          break
             
       if case("GDMLPolycone") :
          #print("      GDMLPolycone") 
          return(GDMLPolyconeObjectVariables(obj))
          break
             
       if case("GDMLPolyhedra") :
          #print("      GDMLPolyhedra") 
          return(GDMLPolyhedraObjectVariables(obj))
          break
             
       if case("GDMLSphere") :
          #print("      GDMLSphere") 
          return(GDMLSphereObjectVariables(obj))
          break

       if case("GDMLTessellated") :
          #print("      GDMLTessellated") 
          return(GDMLTessellatedObjectVariables(obj))
          break

       if case("GDMLTetra") :
          #print("      GDMLTetra") 
          return(GDMLTetraObjectVariables(obj))
          break

       if case("GDMLTorus") :
          print("      GDMLTorus") 
          return(GDMLTorusObjectVariables(obj))
          break

       if case("GDMLTrap") :
          #print("      GDMLTrap") 
          return(GDMLTrapObjectVariables(obj))
          break

       if case("GDMLTrd") :
          #print("      GDMLTrd") 
          return(GDMLTrdObjectVariables(obj))
          break

       if case("GDMLTube") :
          #print("      GDMLTube") 
          return(GDMLTubeObjectVariables(obj))
          break

       if case("GDMLXtru") :
          #print("      GDMLXtru") 
          return(GDMLXtruObjectVariables(obj))
          break

       print("Not yet Handled")
       break  

def processSolid(obj, fp) :
    # export solid & return Name
    while switch(obj.TypeId) :

        if case("Part::FeaturePython"):
            #print("   Python Feature")
            if hasattr(obj.Proxy, 'Type') :
                #print(obj.Proxy.Type) 
                return(GDMLSolidVariables(obj))

        # Test for Booleans
        if case("Part::Fuse") :
           print("   Union")
           unionName = 'Union'+obj.Name
           ref1 = processSolid(obj.Base, fp)
           ref2 = processSolid(obj.Tool, fp)
           #union = ET.SubElement(solids,'union',{'name': unionName })
           #ET.SubElement(union,'first', {'ref': ref1})
           #ET.SubElement(union,'second',{'ref': ref2})
           processPosition(obj,pvol)
           processRotation(obj,pvol)
           processPosition(obj.Tool,union)
           processRotation(obj.Tool,union)
           return True, unionName
           break

        if case("Part::Common") :
           print("   Intersection")
           intersectName = 'Intersect'+obj.Name
           ref1 = processSolid(obj.Base, fp)
           ref2 = processSolid(obj.Tool, fp)
           #intersect = ET.SubElement(solids,'intersection',{'name': \
           #          intersectName })
           #ET.SubElement(intersect,'first', {'ref': ref1})
           #ET.SubElement(intersect,'second',{'ref': ref2})
           processPosition(obj,pvol)
           processRotation(obj,pvol)
           processPosition(obj.Tool,intersect)
           processRotation(obj.Tool,intersect)
           return True, intersectName
           break

        if case("Part::MultiFuse") :
           print("   Multifuse") 
           # test and fix
           multName = 'MultiFuse'+obj.Name
           print('Output Solids')
           for sub in obj.OutList:
               processSolid(sub,fp)
           print('Output Solids Complete')
           #multUnion = ET.SubElement(solids,'multiUnion',{'name': multName })
           for sub in obj.OutList:
               print(sub.Name)
               node = processMuNod(multUnion, 'node-'+str(num))
               #ET.SubElement(node,'solid',{'ref':sub.Name})
               processPosition(sub,node)
               processRotation(sub,node)
           return multName
           break

        if case("Part::MultiCommon") :
           print("   Multi Common / intersection")
           print("   Not available in GDML")
           exit(-3)
           break

        if case("Mesh::Feature") :
           print("   Mesh Feature") 
           # test and Fix
           return(processMesh(obj, obj.Mesh, obj.Name))
           break

        if case("Part::FeaturePython"):
           #print("   Python Feature")
           if hasattr(obj.Proxy, 'Type') :
              #print(obj.Proxy.Type) 
              solidName = processSolid(obj, fp)
              processPosition(obj,pvol)
              processRotation(obj,pvol)
           break  

        # Same as Part::Feature but no position
        if case("App::FeaturePython") :
           print("App::FeaturePython") 
        
        #  
        #  Now deal with objects that map to GDML solids
        #
        if case("Part::Box") :
            print("    Box")
            return(processBoxObject(obj))
            break

        if case("Part::Cylinder") :
            print("    Cylinder")
            return(processCylinderObject(obj))
            break

        if case("Part::Cone") :
            print("    Cone")
            return(processConeObject(obj))
            break

        if case("Part::Sphere") :
            print("    Sphere")
            return(processSphereObject(obj))
            break

        # Not a Solid that translated to GDML solid
        # Dropped through so treat object as a shape
        # Need to check obj has attribute Shape
        # Create tessellated solid
        #
        #return(processObjectShape(obj, addVolsFlag))
        print("Convert FreeCAD shape to GDML Tessellated")
        print(obj.TypeId)
        if hasattr(obj,'Shape') :
           if obj.Shape.isValid() : 
              return(processObjectShape(obj))
        break

import collections
from itertools import islice

def consume(iterator) :
    next(islice(iterator,2,2), None)

def createWorldVol(volName,fp) :
    print("Need to create Dummy Volume and World Box ")
    bbox = FreeCAD.BoundBox()
    boxName = defineWorldBox(bbox,fp)
    fp.write(':VOLU '+volName+' '+boxname+' '+'G4_Galactic') 
    print("Need to FIX !!!! To use defined gas")

def checkGDMLstructure(objList) :
    # Should be 
    # World Vol - App::Part
    # App::Origin
    # GDML Object
    if objList[0].TypeId != 'App::Origin' \
        or objList[2].TypeId != 'App::Part' :
            return False
    return True

def solidDefinition(obj) :
    return('Define '+obj.Name, 0)

def processVolPart(obj, fp ) :
    count = 0
    if hasattr(obj,'OutList') :
       for subobj in obj.OutList :
           if subobj.TypeId == 'Part::FeaturePython' :
              count = count + 1
       print('Number of Solids : '+str(count))
       for subobj in obj.OutList :
           while switch(subobj.TypeId) :
              if case('App::Part') :
                 print('Sub Vol/Part')
                 processVolPart(subobj, fp)
                 break
             
              if case('App::Origin') :
                 pass
                 break

              solidDef, delta = solidDefinition(subobj)
              if count == 1 :
                 # Does Volume Name match Object Name
                 if subobj.Name == obj.Name :
                    # Output as Combined Vol Object
                    fp.write('VOLU: '+ solidDef)
                 else :
                    fp.write('VOL: ' + obj.Name)
                    fp.write('SOL: ' + solidDef) 
              
              else :
                 fp.write('VOLU: LV'+ subobj.Name + ' ' + solidDef)

              break

def exportGAMOS(first,filename) :
       # GAMOS Export
       print("\nStart GDML GAMOS Export 0.1")

       fp = pythonopen(filename,'w')
       print('Process Materials')
       processMaterialsGroupObjects(fp)
       print('Process Structure')
       processVolPart(first,fp)
       print("GAMOS file written")

def export(exportList,filename) :
    "called when FreeCAD exports a file"
    
    first = exportList[0]
    if first.TypeId == "App::Part" :
       exportGAMOS(first,filename)

       #from PyQt4 import QtGui
       #QtGui.QMessageBox.critical(None,'Need to select a Part for export','Press OK')

