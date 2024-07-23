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
    for obj in part.OutList:
        if hasattr(obj, 'Shape') and obj in outerShapes:
            vol = obj.Shape.Volume
            cm += vol * obj.Shape.CenterOfGravity
            totVol += vol
        elif obj.TypeId == 'App::Part':
            partVol, cm0 = partCM(obj)
            cm += partVol * cm0
            totVol += partVol
    
    cm = part.Placement * (cm / totVol)
    return (totVol, cm)


obj = App.Gui.Selection.getSelection()[0]
cm = Vector(0, 0, 0)
totVol = 0
for part in obj.OutList:
    if part.TypeId == 'App::Part':
        vol, cm0 = partCM(part)
        cm0_in_cm = cm0/10
        print(f'{part.Label} volume: {vol/1000:.2f} cm^3  '
              f'--> center of "mass" : ({cm0_in_cm.x:.2f}, {cm0_in_cm.y:.2f},  {cm0_in_cm.z:.2f}) cm')
        cm += vol * cm0
        totVol += vol

cm = cm / totVol
print(f'Center of Geometry =  ({cm.x/10:.2f}, {cm.y/10:.2f}, {cm.z/10:.2f})  cm')
print(f'Total Volume =  {totVol/1000} cm^3')
print('Moments of Intertia:')

