import FreeCAD
from .exportGDML import exportDefineVertex, exportDefine
from . import GDMLShared
from .exportGDML import solids
from .exportGDML import exportPosition, exportRotation
from FreeCAD import Vector

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


class SolidExporter:
    # Abstract class to export object as gdml
    solidExporters = {
        "GDMLArb8": "GDMLArb8Exporter",
        "GDMLBox": "GDMLBoxExporter",
        "GDMLCone": "GDMLConeExporter",
        "GDMLcutTube": "GDMLcutTubeExporter",
        "GDMLElCone": "GDMLElConeExporter",
        "GDMLEllipsoid": "GDMLEllipsoidExporter",
        "GDMLElTube": "GDMLElTubeExporter",
        "GDMLHype": "GDMLHypeExporter",
        "GDMLOrb": "GDMLOrbExporter",
        "GDMLPara": "GDMLParaExporter",
        "GDMLParaboloid": "GDMLParaboloidExporter",
        "GDMLPolycone": "GDMLPolyconeExporter",
        "GDMLGenericPolycone": "GDMLGenericPolyconeExporter",
        "GDMLPolyhedra": "GDMLPolyhedraExporter",
        "GDMLGenericPolyhedra": "GDMLGenericPolyhedraExporter",
        "GDMLSphere": "GDMLSphereExporter",
        "GDMLTessellated": "GDMLTessellatedExporter",
        "GDMLGmshTessellated": "GDMLGmshTessellatedExporter",
        "GDMLTetra": "GDMLTetraExporter",
        "GDMLTetrahedron": "GDMLTetrahedronExporter",
        "GDMLTorus": "GDMLTorusExporter",
        "GDMLTrap": "GDMLTrapExporter",
        "GDMLTrd": "GDMLTrdExporter",
        "GDMLTube": "GDMLTubeExporter",
        "GDMLTwistedbox": "GDMLTwistedboxExporter",
        "GDMLTwistedtrap": "GDMLTwistedtrapExporter",
        "GDMLTwistedtrd": "GDMLTwistedtrdExporter",
        "GDMLTwistedtubs": "GDMLTwistedtubsExporter",
        "GDMLXtru": "GDMLXtruExporter",
        "Part::MultiFuse": "MultiFuseExporter",
        "Part::Extrusion": "ExtrusionExporter",
        "Part::Revolution": "RevolutionExporter",
        "Part::Box": "BoxExporter",
        "Part::Cylinder": "CylinderExporter",
        "Part::Cone": "ConeExporter",
        "Part::Sphere": "SphereExporter",
        "Part::Cut": "BooleanExporter",
        "Part::Fuse": "BooleanExporter",
        "Part::Common": "BooleanExporter"}

    @staticmethod
    def isSolid(obj):
        print(f'isSolid {obj.Label}')
        if obj.TypeId == "Part::FeaturePython":
            typeId = obj.Proxy.Type
            if typeId == 'Array':
                if obj.ArrayType == 'ortho':
                    return True
                elif obj.ArrayType == 'polar':
                    return True
            else:
                return obj.Proxy.Type in SolidExporter.solidExporters
        else:
            return obj.TypeId in SolidExporter.solidExporters

    @staticmethod
    def getExporter(obj):
        from .revolveExporter import RevolutionExporter
        from .extrusionExporter import ExtrusionExporter
        if obj.TypeId == "Part::FeaturePython":
            typeId = obj.Proxy.Type
            if typeId == 'Array':
                if obj.ArrayType == 'ortho':
                    return OrthoArrayExporter(obj)
                elif obj.ArrayType == 'polar':
                    return PolarArrayExporter(obj)
        else:
            typeId = obj.TypeId

        if typeId in SolidExporter.solidExporters:
            classname = SolidExporter.solidExporters[typeId]
            # kludge for classes imported from another module
            # globals["RevolutionExporter"] returns key error
            if classname == "ExtrusionExporter":
                return ExtrusionExporter(obj)
            elif classname == "RevolutionExporter":
                return RevolutionExporter(obj)
            else:
                klass = globals()[classname]
                return klass(obj)
        else:
            print(f'{obj.Label} does not have a Solid Exporter')
            return None

    def __init__(self, obj):
        self.obj = obj

    def name(self):
        return self.obj.Label

    def position(self):
        return self.obj.Placement.Base

    def rotation(self):
        return self.obj.Placement.Rotation

    def placement(self):
        return FreeCAD.Placement(self.position(), self.rotation())

    def export(self):
        print('This abstract base')
        return


class BoxExporter(SolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'box', {'name': self.name(),
                                      'x': str(self.obj.Length.Value),
                                      'y': str(self.obj.Width.Value),
                                      'z': str(self.obj.Height.Value),
                                      'lunit': 'mm'})

    def position(self):
        delta = FreeCAD.Vector(self.obj.Length.Value / 2,
                               self.obj.Width.Value / 2,
                               self.obj.Height.Value / 2)
        # The interpretation of placement for non GDML Part::xxxx objects
        # (box, tube, etc)
        # is different from that of their GDML counterparts.
        # For GDML solids, the rotation is applied first (about the center of
        # the solid)
        # then the rotated solid is translated.
        # For Part::xxxx solids, the part is first translated, then rotated.
        # To reproduce the GDML behavior we have to produce a rotated position
        # to export to gdml
        unrotatedpos = self.obj.Placement.Base + delta
        pos = self.obj.Placement.Rotation*unrotatedpos
        return pos


class CylinderExporter(SolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        # Needs unique Name
        # This is for non GDML cylinder/tube
        ET.SubElement(solids, 'tube', {'name': self.name(),
                                       'rmax': str(self.obj.Radius.Value),
                                       'deltaphi': str(float(self.obj.Angle.Value)),
                                       'aunit': 'deg',
                                       'z': str(self.obj.Height.Value),
                                       'lunit': 'mm'})

    def position(self):
        delta = FreeCAD.Vector(0, 0, self.obj.Height.Value / 2)
        # see comments in BoxExporter
        unrotatedpos = self.obj.Placement.Base + delta
        pos = self.obj.Placement.Rotation*unrotatedpos
        return pos


class ConeExporter(SolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'cone', {
            'name': self.name(),
            'rmax1': str(self.obj.Radius1.Value),
            'rmax2': str(self.obj.Radius2.Value),
            'deltaphi': str(float(self.obj.Angle.Value)),
            'aunit': 'deg',
            'z': str(self.obj.Height.Value),
            'lunit': 'mm'})

    def position(self):
        # Adjustment for position in GDML
        delta = FreeCAD.Vector(0, 0, self.obj.Height.Value / 2)
        # see comments in BoxExporter
        unrotatedpos = self.obj.Placement.Base + delta
        pos = self.obj.Placement.Rotation*unrotatedpos
        return pos


class SphereExporter(SolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'sphere', {
            'name': self.name(),
            'rmax': str(self.obj.Radius.Value),
            'starttheta': str(90.-float(self.obj.Angle2.Value)),
            'deltatheta': str(float(self.obj.Angle2.Value - self.obj.Angle1.Value)),
            'deltaphi': str(float(self.obj.Angle3.Value)),
            'aunit': 'deg',
            'lunit': 'mm'})

    def position(self):
        # see comments in processBoxObject
        unrotatedpos = self.obj.Placement.Base
        pos = self.obj.Placement.Rotation*unrotatedpos
        return pos


class BooleanExporter(SolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def isBoolean(self, obj):
        id = obj.TypeId
        return (id == "Part::Cut" or id == "Part::Fuse" or
                id == "Part::Common")

    def boolOperation(self, obj):
        opsDict = {"Part::Cut": 'subtraction',
                   "Part::Fuse": 'union',
                   "Part::Common": 'intersection'}
        if obj.TypeId in opsDict:
            return opsDict[obj.TypeId]
        else:
            print(f'Boolean type {obj.TypId} not handled yet')
            return None

    def export(self):
        '''
        In FreeCAD doc booleans that are themselves composed of other booleans
        are listed in sequence, eg:
                  topBool:
                      Base: Nonbool_0
                      Tool: bool1:
                            Base: bool2:
                                  Base: Nonbool_1
                                  Tool: Nonbool_2
                            Tool: bool3:
                                  Base: Nonbool_3
                                  Tool: Nonbool_4
        In the gdml file, boolean solids must always refer to PREVIOUSLY
        defined solids. So the last booleans must be written first:
        <Nonbool_0 />
        <Nonbool_1 />
        <Nonbool_2 />
        <Nonbool_3 />
        <Nonbool_4 />
        <bool3: 1st=Nonbool_3, 2nd=Nonbool_4 />
        <bool2: 1st=Nonbool_1, 2nd=Nonbool_2 />
        <bool1: 1st=bool2, 2nd=bool3 />
        <TopBool: 1st=Nonbool_0, 2nd=bool1 />

        The code below first builds the list of booleans in order:
        [topBool, bool1, bool2, bool3]

        Then outputs them to gdml in reverse order.
        In the process of scanning for booleans, the Nonbooleans are exported
        '''
        GDMLShared.trace('Process Boolean Object')

        obj = self.obj
        boolsList = [obj]  # list of booleans that are part of obj
        # dynamic list the is used to figure out when we've iterated over all
        #  subobjects that are booleans
        tmpList = [obj]
        ref1 = {}  # first solid exporter
        ref2 = {}  # second solid exporter
        while(len(tmpList) > 0):
            obj1 = tmpList.pop()
            solidExporter = SolidExporter.getExporter(obj1.Base)
            ref1[obj1] = solidExporter
            if self.isBoolean(obj1.Base):
                tmpList.append(obj1.Base)
                boolsList.append(obj1.Base)
            else:
                solidExporter.export()

            solidExporter = SolidExporter.getExporter(obj1.Tool)
            ref2[obj1] = solidExporter
            if self.isBoolean(obj1.Tool):
                tmpList.append(obj1.Tool)
                boolsList.append(obj1.Tool)
            else:
                solidExporter.export()

        # Now tmpList is empty and boolsList has list of all booleans
        for obj1 in reversed(boolsList):
            operation = self.boolOperation(obj1)
            if operation is None:
                continue
            solidName = obj1.Label
            boolXML = ET.SubElement(solids, str(operation), {
                'name': solidName})
            ET.SubElement(boolXML, 'first', {'ref': ref1[obj1].name()})
            ET.SubElement(boolXML, 'second', {'ref': ref2[obj1].name()})
            # process position & rotation
            pos = ref2[obj1].position()
            exportPosition(ref2[obj1].name(), boolXML, pos)
            # For booleans, gdml want actual rotation, not reverse
            # processRotation export negative of rotation angle(s)
            # This is ugly way of NOT reversing angle:

            rot = FreeCAD.Rotation()
            rot.Axis = ref2[obj1].rotation().Axis
            angle = ref2[obj1].rotation().Angle
            rot.Angle = -angle
            exportRotation(ref2[obj1].name(), boolXML, rot)


class GDMLSolidExporter(SolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def name(self):
        from .exportGDML import nameOfGDMLobject
        return nameOfGDMLobject(self.obj)


class GDMLArb8Exporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'arb8', {'name': self.name(),
                                       'v1x': str(self.obj.v1x),
                                       'v1y': str(self.obj.v1y),
                                       'v2x': str(self.obj.v2x),
                                       'v2y': str(self.obj.v2y),
                                       'v3x': str(self.obj.v3x),
                                       'v3y': str(self.obj.v3y),
                                       'v4x': str(self.obj.v4x),
                                       'v4y': str(self.obj.v4y),
                                       'v5x': str(self.obj.v5x),
                                       'v5y': str(self.obj.v5y),
                                       'v6x': str(self.obj.v6x),
                                       'v6y': str(self.obj.v6y),
                                       'v7x': str(self.obj.v7x),
                                       'v7y': str(self.obj.v7y),
                                       'v8x': str(self.obj.v8x),
                                       'v8y': str(self.obj.v8y),
                                       'dz': str(self.obj.dz),
                                       'lunit': self.obj.lunit})


class GDMLBoxExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'box', {'name': self.name(),
                                      'x': str(self.obj.x),
                                      'y': str(self.obj.y),
                                      'z': str(self.obj.z),
                                      'lunit': self.obj.lunit})


class GDMLConeExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'cone', {'name': self.name(),
                                       'rmin1': str(self.obj.rmin1),
                                       'rmin2': str(self.obj.rmin2),
                                       'rmax1': str(self.obj.rmax1),
                                       'rmax2': str(self.obj.rmax2),
                                       'startphi': str(self.obj.startphi),
                                       'deltaphi': str(self.obj.deltaphi),
                                       'aunit': self.obj.aunit,
                                       'z': str(self.obj.z),
                                       'lunit': self.obj.lunit})


class GDMLCutTubeExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'cutTube', {'name': self.name(),
                                          'rmin': str(self.obj.rmin),
                                          'rmax': str(self.obj.rmax),
                                          'startphi': str(self.obj.startphi),
                                          'deltaphi': str(self.obj.deltaphi),
                                          'aunit': self.obj.aunit,
                                          'z': str(self.obj.z),
                                          'highX': str(self.obj.highX),
                                          'highY': str(self.obj.highY),
                                          'highZ': str(self.obj.highZ),
                                          'lowX': str(self.obj.lowX),
                                          'lowY': str(self.obj.lowY),
                                          'lowZ': str(self.obj.lowZ),
                                          'lunit': self.obj.lunit})


class GDMLElConeExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'elcone', {'name': self.name(),
                                         'dx': str(self.obj.dx),
                                         'dy': str(self.obj.dy),
                                         'zcut': str(self.obj.zcut),
                                         'zmax': str(self.obj.zmax),
                                         'lunit': str(self.obj.lunit)})


class GDMLEllipsoidExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'ellipsoid', {'name': self.name(),
                                            'ax': str(self.obj.ax),
                                            'by': str(self.obj.by),
                                            'cz': str(self.obj.cz),
                                            'zcut1': str(self.obj.zcut1),
                                            'zcut2': str(self.obj.zcut2),
                                            'lunit': self.obj.lunit})


class GDMLElTubeExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'eltube', {'name': self.name(),
                                         'dx': str(self.obj.dx),
                                         'dy': str(self.obj.dy),
                                         'dz': str(self.obj.dz),
                                         'lunit': self.obj.lunit})


class GDMLHypeExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'hype', {'name': self.name(),
                                       'rmin': str(self.obj.rmin),
                                       'rmax': str(self.obj.rmax),
                                       'z': str(self.obj.z),
                                       'inst': str(self.obj.inst),
                                       'outst': str(self.obj.outst),
                                       'aunit': self.obj.aunit,
                                       'lunit': self.obj.lunit})


class GDMLParaboloidExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'paraboloid', {'name': self.name(),
                                             'rlo': str(self.obj.rlo),
                                             'rhi': str(self.obj.rhi),
                                             'dz': str(self.obj.dz),
                                             'lunit': self.obj.lunit})


class GDMLOrbExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'orb', {'name': self.name(),
                                      'r': str(self.obj.r),
                                      'lunit': self.obj.lunit})


class GDMLParaExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'para', {'name': self.name(),
                                       'x': str(self.obj.x),
                                       'y': str(self.obj.y),
                                       'z': str(self.obj.z),
                                       'alpha': str(self.obj.alpha),
                                       'theta': str(self.obj.theta),
                                       'phi': str(self.obj.phi),
                                       'aunit': str(self.obj.aunit),
                                       'lunit': self.obj.lunit})


class GDMLPolyconeExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        cone = ET.SubElement(solids, 'polycone', {
            'name': self.name(),
            'startphi': str(self.obj.startphi),
            'deltaphi': str(self.obj.deltaphi),
            'aunit': self.obj.aunit,
            'lunit': self.obj.lunit})

        for zplane in self.obj.OutList:
            ET.SubElement(cone, 'zplane', {'rmin': str(zplane.rmin),
                                           'rmax': str(zplane.rmax),
                                           'z': str(zplane.z)})


class GDMLGenericPolyconeExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        cone = ET.SubElement(solids, 'genericPolycone', {
            'name': self.name(),
            'startphi': str(self.obj.startphi),
            'deltaphi': str(self.obj.deltaphi),
            'aunit': self.obj.aunit,
            'lunit': self.obj.lunit})
        for rzpoint in self.obj.OutList:
            ET.SubElement(cone, 'rzpoint', {'r': str(rzpoint.r),
                                            'z': str(rzpoint.z)})


class GDMLGenericPolyhedraExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        polyhedra = ET.SubElement(solids, 'genericPolyhedra', {
            'name': self.name(),
            'startphi': str(self.obj.startphi),
            'deltaphi': str(self.obj.deltaphi),
            'numsides': str(self.obj.numsides),
            'aunit': self.obj.aunit,
            'lunit': self.obj.lunit})
        for rzpoint in self.obj.OutList:
            ET.SubElement(polyhedra, 'rzpoint', {'r': str(rzpoint.r),
                                                 'z': str(rzpoint.z)})


class GDMLPolyhedraExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        poly = ET.SubElement(solids, 'polyhedra', {
            'name': self.name(),
            'startphi': str(self.obj.startphi),
            'deltaphi': str(self.obj.deltaphi),
            'numsides': str(self.obj.numsides),
            'aunit': self.obj.aunit,
            'lunit': self.obj.lunit})

        for zplane in self.obj.OutList:
            ET.SubElement(poly, 'zplane', {'rmin': str(zplane.rmin),
                                           'rmax': str(zplane.rmax),
                                           'z': str(zplane.z)})


class GDMLSphereExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'sphere', {
            'name': self.name(),
            'rmin': str(self.obj.rmin),
            'rmax': str(self.obj.rmax),
            'startphi': str(self.obj.startphi),
            'deltaphi': str(self.obj.deltaphi),
            'starttheta': str(self.obj.starttheta),
            'deltatheta': str(self.obj.deltatheta),
            'aunit': self.obj.aunit,
            'lunit': self.obj.lunit})


class GDMLTessellatedExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        tessName = self.name()
        # Use more readable version
        tessVname = tessName + '_'
        # print(dir(obj))
        vertexHashcodeDict = {}

        tess = ET.SubElement(solids, 'tessellated', {'name': tessName})
        for i, v in enumerate(self.obj.Shape.Vertexes):
            vertexHashcodeDict[v.hashCode()] = i
            exportDefineVertex(tessVname, v, i)

        for f in self.obj.Shape.Faces:
            # print(f'Normal at : {n} dot {dot} {clockWise}')
            vertexes = f.OuterWire.OrderedVertexes
            if len(f.Edges) == 3:
                i0 = vertexHashcodeDict[vertexes[0].hashCode()]
                i1 = vertexHashcodeDict[vertexes[1].hashCode()]
                i2 = vertexHashcodeDict[vertexes[2].hashCode()]
                ET.SubElement(tess, 'triangular', {
                    'vertex1': tessVname+str(i0),
                    'vertex2': tessVname+str(i1),
                    'vertex3': tessVname+str(i2),
                    'type': 'ABSOLUTE'})
            elif len(f.Edges) == 4:
                i3 = vertexHashcodeDict[vertexes[3].hashCode()]
                ET.SubElement(tess, 'quadrangular', {
                    'vertex1': tessVname+str(i0),
                    'vertex2': tessVname+str(i1),
                    'vertex3': tessVname+str(i2),
                    'vertex4': tessVname+str(i3),
                    'type': 'ABSOLUTE'})


class GDMLTetraExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        tetraName = self.name()
        v1Name = tetraName + 'v1'
        v2Name = tetraName + 'v2'
        v3Name = tetraName + 'v3'
        v4Name = tetraName + 'v4'
        exportDefine(v1Name, self.obj.v1)
        exportDefine(v2Name, self.obj.v2)
        exportDefine(v3Name, self.obj.v3)
        exportDefine(v4Name, self.obj.v4)

        ET.SubElement(solids, 'tet', {'name': tetraName,
                                      'vertex1': v1Name,
                                      'vertex2': v2Name,
                                      'vertex3': v3Name,
                                      'vertex4': v4Name})


class GDMLTetrahedronExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        global structure
        global solids
        tetrahedronName = self.name()
        print('Len Tet' + str(len(self.obj.Proxy.Tetra)))
        count = 0
        for t in self.obj.Proxy.Tetra:
            tetraName = 'Tetra_' + str(count)
            v1Name = tetraName + 'v1'
            v2Name = tetraName + 'v2'
            v3Name = tetraName + 'v3'
            v4Name = tetraName + 'v4'
            exportDefine(v1Name, t[0])
            exportDefine(v2Name, t[1])
            exportDefine(v3Name, t[2])
            exportDefine(v4Name, t[3])
            ET.SubElement(solids, 'tet', {'name': tetraName,
                                          'vertex1': v1Name,
                                          'vertex2': v2Name,
                                          'vertex3': v3Name,
                                          'vertex4': v4Name})
            lvName = 'LVtetra' + str(count)
            lvol = ET.SubElement(structure, 'volume', {'name': lvName})
            ET.SubElement(lvol, 'materialref', {'ref': self.obj.material})
            ET.SubElement(lvol, 'solidref', {'ref': tetraName})
            count += 1

        # Now put out Assembly
        assembly = ET.SubElement(structure, 'assembly', {
            'name': tetrahedronName})
        count = 0
        for t in self.obj.Proxy.Tetra:
            lvName = 'LVtetra' + str(count)
            physvol = ET.SubElement(assembly, 'physvol')
            ET.SubElement(physvol, 'volumeref', {'ref': lvName})
            # ET.SubElement(physvol, 'position')
            # ET.SubElement(physvol, 'rotation')
            count += 1


class GDMLTorusExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'torus', {'name': self.name(),
                                        'rmin': str(self.obj.rmin),
                                        'rmax': str(self.obj.rmax),
                                        'rtor': str(self.obj.rtor),
                                        'startphi': str(self.obj.startphi),
                                        'deltaphi': str(self.obj.deltaphi),
                                        'aunit': self.obj.aunit,
                                        'lunit': self.obj.lunit})


class GDMLTrapExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'trap', {'name': self.name(),
                                       'z': str(self.obj.z),
                                       'theta': str(self.obj.theta),
                                       'phi': str(self.obj.phi),
                                       'x1': str(self.obj.x1),
                                       'x2': str(self.obj.x2),
                                       'x3': str(self.obj.x3),
                                       'x4': str(self.obj.x4),
                                       'y1': str(self.obj.y1),
                                       'y2': str(self.obj.y2),
                                       'alpha1': str(self.obj.alpha),
                                       'alpha2': str(self.obj.alpha),
                                       'aunit': self.obj.aunit,
                                       'lunit': self.obj.lunit})


class GDMLTrdExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'trd', {'name': self.name(),
                                      'z': str(self.obj.z),
                                      'x1': str(self.obj.x1),
                                      'x2': str(self.obj.x2),
                                      'y1': str(self.obj.y1),
                                      'y2': str(self.obj.y2),
                                      'lunit': self.obj.lunit})


class GDMLTubeExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'tube', {'name': self.name(),
                                       'rmin': str(self.obj.rmin),
                                       'rmax': str(self.obj.rmax),
                                       'startphi': str(self.obj.startphi),
                                       'deltaphi': str(self.obj.deltaphi),
                                       'aunit': self.obj.aunit,
                                       'z': str(self.obj.z),
                                       'lunit': self.obj.lunit})


class GDMLTwistedboxExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'twistedbox', {
            'name': self.name(),
            'PhiTwist': str(self.obj.PhiTwist),
            'x': str(self.obj.x),
            'y': str(self.obj.y),
            'z': str(self.obj.z),
            'aunit': str(self.obj.aunit),
            'lunit': self.obj.lunit})


class GDMLTwistedtrdExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'twistedtrd', {
            'name': self.name(),
            'PhiTwist': str(self.obj.PhiTwist),
            'x1': str(self.obj.x1),
            'x2': str(self.obj.x2),
            'y1': str(self.obj.y1),
            'y2': str(self.obj.y2),
            'z': str(self.obj.z),
            'aunit': str(self.obj.aunit),
            'lunit': self.obj.lunit})


class GDMLTwistedtrapExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'twistedtrap', {
            'name': self.name(),
            'PhiTwist': str(self.obj.PhiTwist),
            'x1': str(self.obj.x1),
            'x2': str(self.obj.x2),
            'y1': str(self.obj.y1),
            'y2': str(self.obj.y2),
            'x3': str(self.obj.x3),
            'x4': str(self.obj.x4),
            'z': str(self.obj.z),
            'Theta': str(self.obj.Theta),
            'Phi': str(self.obj.Phi),
            'Alph': str(self.obj.Alph),
            'aunit': str(self.obj.aunit),
            'lunit': self.obj.lunit})


class GDMLTwistedtubsExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'twistedtubs', {
            'name': self.name(),
            'twistedangle': str(self.obj.twistedangle),
            'endinnerrad': str(self.obj.endinnerrad),
            'endouterrad': str(self.obj.endouterrad),
            'zlen': str(self.obj.zlen),
            'phi': str(self.obj.phi),
            'aunit': str(self.obj.aunit),
            'lunit': self.obj.lunit})


class GDMLXtruExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        xtru = ET.SubElement(solids, 'xtru', {'name': self.name(),
                                              'lunit': self.obj.lunit})
        for items in self.obj.OutList:
            if items.Type == 'twoDimVertex':
                ET.SubElement(xtru, 'twoDimVertex', {'x': str(items.x),
                                                     'y': str(items.y)})
            if items.Type == 'section':
                ET.SubElement(xtru, 'section', {
                    'zOrder': str(items.zOrder),
                    'zPosition': str(items.zPosition),
                    'xOffset': str(items.xOffset),
                    'yOffset': str(items.yOffset),
                    'scalingFactor': str(items.scalingFactor)})


class GDML2dVertexExporter(GDMLSolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def export(self):
        ET.SubElement(solids, 'twoDimVertex', {'x': self.obj.x,
                                               'y': self.obj.y})


class MultiFuseExporter(SolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def name(self):
        solidName = 'MultiFuse' + self.obj.Label
        return solidName

    def export(self):
        from .exportGDML import processPlacement
        GDMLShared.trace("Multifuse - multiunion")
        # test and fix
        # First add solids in list before reference
        print('Output Solids')
        exporters = []
        for sub in self.obj.OutList:
            exporter = SolidExporter.getExporter(sub)
            if exporter is not None:
                exporter.export()
                exporters.append(exporter)

        GDMLShared.trace('Output Solids Complete')
        multUnion = ET.SubElement(solids, 'multiUnion', {
            'name': self.name()})

        num = 1
        for exp in exporters:
            GDMLShared.trace(exp.name())
            node = ET.SubElement(multUnion, 'multiUnionNode', {
                'name': 'node-'+str(num)})
            ET.SubElement(node, 'solid', {'ref': exp.name()})
            processPlacement(exp.name(), node, exp.placement)
            num += 1

        GDMLShared.trace('Return MultiFuse')


class OrthoArrayExporter(SolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def name(self):
        solidName = 'MultiUnion-' + self.obj.Label
        return solidName

    def export(self):
        base = self.obj.OutList[0]
        baseExporter = SolidExporter.getExporter(base)
        baseExporter.export()
        volRef = baseExporter.name()
        unionXML = ET.SubElement(solids, 'multiUnion', {'name': self.name()})
        # structure.insert(0, assembly)
        translate = Vector(0, 0, 0)
        for ix in range(self.obj.NumberX):
            translate = ix*self.obj.IntervalX
            for iy in range(self.obj.NumberY):
                translate += iy*self.obj.IntervalY
                for iz in range(self.obj.NumberZ):
                    nodeName = f'{self.name()}_{ix}_{iy}_{iz}'
                    translate += iz*self.obj.IntervalZ
                    nodeXML = ET.SubElement(unionXML, 'multiUnionNode', {'name': nodeName})
                    ET.SubElement(nodeXML, 'solid', {'ref': volRef})
                    ET.SubElement(nodeXML, 'position', {
                        'x': str(translate.x),
                        'y': str(translate.y),
                        'z': str(translate.z),
                        'unit': 'mm'})


class PolarArrayExporter(SolidExporter):
    def __init__(self, obj):
        super().__init__(obj)

    def name(self):
        solidName = 'MultiUnion-' + self.obj.Label
        return solidName

    def export(self):
        base = self.obj.OutList[0]
        baseExporter = SolidExporter.getExporter(base)
        baseExporter.export()
        volRef = baseExporter.name()
        unionXML = ET.SubElement(solids, 'multiUnion', {'name': self.name()})
        dthet = self.obj.Angle/self.obj.NumberPolar
        positionVector = self.obj.Base.Placement.Base
        axis = self.obj.Axis
        # TODO adjust for center of rotation != origin
        for i in range(self.obj.NumberPolar):
            rot = FreeCAD.Rotation(axis, i*dthet)
            pos = rot*positionVector     # position has to be roated too!
            rot.Angle = -rot.Angle   # undo angle reversal by exportRotation
            nodeName = f'{self.name()}_{i}'
            nodeXML = ET.SubElement(unionXML, 'multiUnionNode', {'name': nodeName})
            ET.SubElement(nodeXML, 'solid', {'ref': volRef})
            exportPosition(nodeName, nodeXML, pos)
            exportRotation(nodeName, nodeXML, rot)
