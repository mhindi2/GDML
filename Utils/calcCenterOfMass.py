#
# First select the World volume and then execute this macro
#

import sys
import FreeCAD as App
from FreeCAD import Vector
import Part


def topShapes(part: App.Part) -> list:
    ''' return a list containing the top shapes contained in part
    a top shape may contain other shapes, but cannot be contained in
    other shapes within Part
    '''
    shapeObjs = []
    for obj in part.OutList:
        if hasattr(obj, 'Shape'):
            print(f"Adding {obj.Label}")
            shapeObjs.append(obj)
    
    # now remove shapes contained in other shapes
    for obj in shapeObjs[:]:
        for parentObj in obj.InList:
            if parentObj in shapeObjs:
                print(f"removing {obj.Label}")
                shapeObjs.remove(obj)
    
    
    return shapeObjs

def partCM(part: App.Part) -> tuple:
    outerShapes = topShapes(part)
    print(outerShapes)
    cm: Vector = Vector(0, 0, 0)
    totVol: float = 0
    for obj in part.OutList:
        if hasattr(obj, 'Shape') and obj in outerShapes:
            vol = obj.Shape.Volume
            cm += vol * obj.Shape.CenterOfGravity
            totVol += vol
            print(f'{obj.Label} {totVol} {cm}')
        elif obj.TypeId == 'App::Part':
            partVol, partCM = partCM(obj)
            cm += partVol * partCM
            totVol += partVol
    
    cm = part.Placement * (cm / totVol)
    return (totVol, cm)

obj = App.Gui.Selection.getSelection()[0]
cm = Vector(0, 0, 0)
totVol = 0
for part in obj.OutList:
    if part.TypeId == 'App::Part':
        vol, partCM = partCM(part) 
        print(f'{obj.Label} vol: {vol} cm: {cm}')
        cm += vol * partCM
        totVol += vol

cm = cm / totVol
print(f'Center of Geometry =  {cm}  mm')
print(f'Total Volume =  {totVol} mm^3')
print('Moments of Intertia:')

