import sys
from lxml import etree 

volList = []
solidList = []
positionList = []
rotationList = []

def processVol(vol) :
   global volList, solidList
   print(vol)
   #print(vol.attrib)
   vname = vol.attrib.get('name')
   print('volume : ' + vname)
   if vname not in volList : volList.append(vname)
   solid = vol.find('solidref')
   #print(solid.attrib)
   sname = solid.attrib.get('ref')
   print('solid : ' + sname)
   if sname not in solidList : solidList.append(sname)
   material = vol.find('materialref')
   if material is not None :
      #print('material : '+str(material.attrib))
      print('material : ' + material.attrib.get('ref'))
   for pv in vol.findall('physvol') :
      volref = pv.find('volumeref')
      print('physvol : '+volref.attrib.get('ref'))
      posref = pv.find('positionref')
      if posref is not None :
         posname = posref.attrib.get('ref')
         print('Position ref : '+posname)
         if posname not in positionList : positionList.append(posname)
      rotref = pv.find('rotationref')
      if rotref is not None :
         rotname = rotref.attrib.get('ref')
         print('Rotation ref : '+rotname)
         if rotname not in rotationList : rotationList.append(rotname)

if len(sys.argv)<4:
  print ("Usage: sys.argv[0] <Volume> <in_file> <Out_file")
  sys.exit(1)

volume = sys.argv[1]
iname = sys.argv[2]
oname = sys.argv[3]

print('\nExtracting Volume : '+volume+' from : '+iname+' to '+oname)
tree = etree.parse(iname)
root = tree.getroot()
structure = tree.find('structure')
#print(etree.fromstring(structure))
# Following works
#vol = structure.find('volume[@name="World"]')
vol = structure.find(f"volume[@name='{volume}']")
if vol is not None :
   processVol(vol)
else :
   print('Volume : '+volume+' Not found')

print('Vol List')
print(volList)
print('Solid List')
print(solidList)
print('Position List')
print(positionList)
print('Rotation List')
print(rotationList)
