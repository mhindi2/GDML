#
# First select the World volume and then execute this macro
#

import sys
import FreeCAD as App
from FreeCAD import Vector
import Part
from FreeCAD import Matrix

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
    for obj in shapeObjs[:]:
        for parentObj in obj.InList:
            if parentObj in shapeObjs:
                # print(f"removing {obj.Label}")
                shapeObjs.remove(obj)
    
    
    return shapeObjs


def partCM(part: App.Part) -> tuple:
    outerShapes = topShapes(part)
    # print(outerShapes)
    cm: Vector = Vector(0, 0, 0)
    totVol: float = 0
    II: Matrix = Matrix(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    for obj in part.OutList:
        if hasattr(obj, 'Shape') and obj in outerShapes:
            vol = obj.Shape.Volume
            cm += vol * obj.Shape.CenterOfGravity
            totVol += vol
            shape = obj.Shape
            if hasattr(shape, 'MatrixOfInertia'):
                A = obj.Shape.MatrixOfInertia.A
                B = (A[i]/vol for i in range(len(A)))
                II += Matrix(*B)
            else:
                print(f"{obj.Label} has not MatrixOfInertia")
        elif obj.TypeId == 'App::Part':
            partVol, cm0, II0 = partCM(obj)
            cm += partVol * cm0
            totVol += partVol
            II += II0
    
    cm = part.Placement * (cm / totVol)
    return (totVol, cm, II)


def prettyPrint(vol, cm0, II0):
    volunit = "mm^3"
    if vol > 1000000:
        vol /= 1000000
        volunit = " m^3"
    elif vol > 10000:
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
        prettyPrint(vol, cm0, II0)

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

