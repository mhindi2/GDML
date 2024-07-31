#
# First select the World volume and then execute this macro
#

import sys
import random
import time

import FreeCAD as App
from FreeCAD import Vector
import Part
from FreeCAD import Matrix
import MeshPart

from dataclasses import dataclass

@dataclass
class Moments:
    Volume: float
    cm: Vector
    M: Matrix


def MonteCarloMoments(shape, nsim: int) -> Moments:
    moments = Moments(0, Vector(0, 0, 0), Matrix())
    b = shape.BoundBox
    print(f"xmin={b.XMin} xmax={b.XMax}")
    pmin = Vector(b.XMin, b.YMin, b.ZMin)
    pmax = Vector(b.XMax, b.YMax, b.ZMax)

    cm = Vector(0, 0, 0)
    Ixx = 0.
    Ixy = 0.
    Ixz = 0.
    Iyy = 0.
    Iyz = 0.
    Izz = 0.

    numIn: int = 0

    start = time.perf_counter()
    for i in range(nsim):
        x = pmin.x + random.random() * (pmax.x - pmin.x)
        y = pmin.y + random.random() * (pmax.y - pmin.y)
        z = pmin.z + random.random() * (pmax.z - pmin.z)
        r = Vector(x, y, z)
        if shape.isInside(r, 0.1, True):
            numIn += 1
            cm += r

            Ixx += z * z + y * y
            Ixy += -x * y
            Ixz += -x * z
            Iyy += x * x + z * z
            Iyz += -y * z
            Izz += x * x + y * y
    end = time.perf_counter()
    print(f" Monte Carlo time = {end-start} seconds")

    # calculate moments about center of mass, using parallel axis theorem.
    cm /= numIn
    Ixx -= numIn * (cm.y * cm.y + cm.z * cm.z)
    Iyy -= numIn * (cm.x * cm.x + cm.z * cm.z)
    Izz -= numIn * (cm.x * cm.x + cm.y * cm.y)
    Ixy += numIn * cm.x * cm.y
    Ixz += numIn * cm.x * cm.z
    Iyz += numIn * cm.y * cm.z

    moments.cm = cm
    moments.M = Matrix(Ixx, Ixy, Ixz, 0,
                       Ixy, Iyy, Iyz, 0,
                       Ixz, Iyz, Izz, 0,
                       0, 0, 0, 1)

    A = moments.M.A
    B = (A[i] / numIn for i in range(len(A)))
    moments.M = Matrix(*B)

    MCVolume = b.XLength * b.YLength * b.ZLength
    print(f"numIn= {numIn} nsim={nsim}")
    MCVolume *= float(numIn) / nsim
    moments.Volume = MCVolume

    return moments


def topShapes(part: App.Part) -> list:
    ''' return a list containing the top shapes contained in part
    a top shape may contain other shapes, but cannot be contained in
    other shapes within Part
    '''
    shapeObjs = []
    for obj in part.OutList:
        if hasattr(obj, 'Shape') and obj.Shape.Volume > 0:
            # print(f"Adding {obj.Label}")
            shapeObjs.append(obj)
    
    # now remove shapes contained in other shapes
    copy = shapeObjs[:]
    for obj in copy:
        for parentObj in obj.InList:
            if parentObj in copy:
                # print(f"removing {obj.Label}")
                shapeObjs.remove(obj)

    return shapeObjs


def partCM(part: App.Part) -> tuple:
    '''
    Note, what we are calling the total moment of inertia is NOT the total moment
    of inertia at all!
    So as not to over-emphasize large bodies compared to small ones in the "total",
    we actually calculate the moment of inertia PER UNIT VOLUME for each object
    The sum of the moments of inertia per unit volume is a useful tool to see
    if the setup in geant matches the setup in FreeCAD, but IT IS NOT a moment of
    inertia, not even, per unit volume. The latter would be sum( volume_i * inertia_density_i)/total volume
    '''
    outerShapes = topShapes(part)
    cm: Vector = Vector(0, 0, 0)
    totVol: float = 0
    II: Matrix = Matrix(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    for obj in part.OutList:
        if hasattr(obj, 'Shape') and obj in outerShapes:
            vol = obj.Shape.Volume
            cm0 = obj.Shape.CenterOfGravity
            cm += vol * cm0
            totVol += vol
            shape = obj.Shape
            if hasattr(shape, 'MatrixOfInertia'):
                II0 = obj.Shape.MatrixOfInertia
                A = II0.A
            else:
                # print(f"{obj.Label} has no MatrixOfInertia. Calculate using MonteCarlo")
                # moments: Moments = MonteCarloMoments(shape, 50000)
                # print(f"{obj.Label} has no MatrixOfInertia. Calculate using MonteCarlo")
                mesh = MeshPart.meshFromShape(Shape=obj.Shape, LinearDeflection=0.1, AngularDeflection=0.524,
                                              Relative=False)
                shape = Part.Shape()
                shape.makeShapeFromMesh(mesh.Topology, 0.05)
                solid = Part.makeSolid(shape)
                II0 = solid.MatrixOfInertia
                A = II0.A
                # print(moments.M)
                # A = moments.M.A
                # B = (A[i] for i in range(len(A)))

            B = (A[i] / vol for i in range(len(A)))
            II0 = Matrix(*B)
            II += II0

            prettyPrint(obj, vol, cm0, II0)

        elif obj.TypeId == 'App::Part':
            partVol, cm0, II0 = partCM(obj)
            cm += partVol * cm0
            totVol += partVol
            II += II0
            prettyPrint(obj, partVol, cm0, II0)

    cm = part.Placement * (cm / totVol)
    return (totVol, cm, II)


def prettyPrint(part, vol, cm0, II0):
    volunit = "mm^3"
    if vol > 1e9:
        vol /= 1e9
        volunit = " m^3"
    elif vol > 1000:
        vol /= 1000
        volunit = "cm^3"

    cmunit = "mm";
    cm = cm0
    if cm0.Length >= 1000.0:
        cm /= 1000
        cmunit = " m"
    elif cm0.Length >= 100.0:
        cm /= 10
        cmunit = "cm"

    inertiaUnit = "mm^2"
    maxInertia = max(II0.A)
    if maxInertia > 1e6:
        B = (II0.A[i] / 1e6 for i in range(len(II0.A)))
        II0 = Matrix(*B)
        inertiaUnit = " m^3"
    elif maxInertia > 1000:
        B = (II0.A[i] / 100 for i in range(len(II0.A)))
        II0 = Matrix(*B)
        inertiaUnit = "cm^3"


    s = f'{part.Label:20s} Volume: {vol:8.2f} {volunit}'
    s += f' --> center of "mass": ({cm.x:7.2f},{cm.y:7.2f},{cm.z:7.2f}) {cmunit}, '
    s += f'moments of inetria: ( ({II0.A11:7.2f},{II0.A12:7.2f},{II0.A13:7.2f}), '
    s += f'({II0.A21:7.2f},{II0.A22:7.2f},{II0.A23:7.2f}), '
    s += f'({II0.A31:7.2f},{II0.A32:7.2f},{II0.A33:7.2f}) ) {inertiaUnit}'
    print(s)


obj = App.Gui.Selection.getSelection()[0]
cm = Vector(0, 0, 0)
totVol = 0
II: Matrix = Matrix(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
for part in obj.OutList:
    if part.TypeId == 'App::Part':
        vol, cm0, II0 = partCM(part)
        prettyPrint(part, vol, cm0, II0)

        cm += vol * cm0
        totVol += vol
        II += II0

cm = cm / totVol
print()
print(f'Center of Geometry =  ({cm.x/10:.2f}, {cm.y/10:.2f}, {cm.z/10:.2f})  cm')
print(f'Total Volume =  {totVol/1000:.2f} cm^3')
s = f'moments of inetria: ({II.A11/100:0.2f}, {II.A12/100:0.2f}, {II.A13/100:0.2f}), '
s += f'({II.A21/100:0.2f}, {II.A22/100:0.2f}, {II.A23/100:0.2f}), '
s += f'({II.A31/100:0.2f}, {II.A32/100:0.2f}, {II.A33/100:0.2f}) cm^2 '
print(s)

