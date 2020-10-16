import sys
from lxml import etree 

volList = []
solidList = []
positionList = []
rotationList = []

def processSolid(sname) :
   # Passed olid name volume, booleans
   global solidList, solids
   solid = oldSolids.find(f"*[@name='{sname}']")
   print('Adding : '+sname)
   print(solid.attrib)
   if sname not in solidList : solidList.append(sname)
   if solid.tag == 'subtraction' or solid.tag == 'union' or \
      solid.tag == 'intersection' :
      print('Found Boolean')
      first = solid.find('first')
      print('first : '+str(first.attrib))
      fname = first.attrib.get('ref')
      processSolid(fname)
      second = solid.find('second')
      print('second : '+str(second.attrib))
      sname = second.attrib.get('ref')
      processSolid(sname)

def processVol(vol) :
   global volList, solidList, oldSolids
   print(vol)
   print(vol.attrib)
   # Need to process physvols first
   for pv in vol.findall('physvol') :
      volref = pv.find('volumeref')
      pname = volref.attrib.get('ref')
      print('physvol : '+pname)
      pvol = structure.find(f"volume[@name='{pname}']")
      if pvol is not None :
         processVol(pvol)
      posref = pv.find('positionref')
      if posref is not None :
         posname = posref.attrib.get('ref')
         print('Position ref : '+posname)
         if posname not in positionList : 
            positionList.append(posname)
      rotref = pv.find('rotationref')
      if rotref is not None :
         rotname = rotref.attrib.get('ref')
         print('Rotation ref : '+rotname)
         if rotname not in rotationList : rotationList.append(rotname)

   vname = vol.attrib.get('name')
   print('volume : ' + vname)
   if vname not in volList : volList.append(vname)
   solid = vol.find('solidref')
   sname = solid.attrib.get('ref')
   processSolid(sname)
   material = vol.find('materialref')
   if material is not None :
      #print('material : '+str(material.attrib))
      print('material : ' + material.attrib.get('ref'))

def processAssembly(assem) :
    name = assem.attrib.get('name')
    print('Process Assembly ; '+name)


if len(sys.argv)<5:
  print ("Usage: sys.argv[0] <parms> <Volume> <in_file> <Out_file> <materials>")
  print("/n For parms the following are or'ed together")
  print("     1 - gdml otherwise xml")
  print("     2 - minmum materials")
  print("     4 - separate materials file") 
  sys.exit(1)

parms  = int(sys.argv[1])
if ( parms and 1 ) == 1 :
   gdml = True 
   xml = etree.Element('gdml')
else :
   xml = etree.Element('xml')

volume = sys.argv[2]
iname = sys.argv[3]
oname = sys.argv[4]

print('\nExtracting Volume : '+volume+' from : '+iname+' to '+oname)
tree = etree.parse(iname)
root = tree.getroot()
structure = tree.find('structure')
oldSolids = tree.find('solids')
#print(etree.fromstring(structure))
# Following works
#vol = structure.find('volume[@name="World"]')
# Test if Volume
vol = structure.find(f"volume[@name='{volume}']")
if vol is not None :
   processVol(vol)
else : 
   # Test if Assembly
   vol = structure.find(f"assembly[@name='{volume}']")
   if vol is not None :
      processAssembly(vol)
   else :
      print(volume+' :  Not found as Volume or Assembly')
      exit(0)
newDefine = etree.SubElement(xml,'define')
mats = tree.find('materials')
xml.append(mats)
newSolids = etree.SubElement(xml,'solids')
newVols   = etree.SubElement(xml,'structure')
oldDefine = tree.find('define') 
oldSolids = tree.find('solids')
oldVols   = tree.find('structure')
for posName in positionList :
    p = oldDefine.find(f"position[@name='{posName}']")
    newDefine.append(p)
for rotName in rotationList :
    p = oldDefine.find(f"rotation[@name='{rotName}']")
    newDefine.append(p)
for solidName in solidList :
    print('Solid : '+solidName)
    s = oldSolids.find(f"*[@name='{solidName}']")
    #print(s.attrib)
    newSolids.append(s)
for volName in volList :
    v = oldVols.find(f"volume[@name='{volName}']")
    newVols.append(v)

print('Vol List')
print(volList)
print('Solid List')
print(solidList)
print('Position List')
print(positionList)
print('Rotation List')
print(rotationList)
if gdml == True :
   set = etree.SubElement(xml,'setup', {'name':'Default', 'version':'1.0'})
   etree.SubElement(set,'world', { 'ref' : volList[-1]})
etree.ElementTree(xml).write(oname)
