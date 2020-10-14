import sys
from lxml import etree 

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
   print(vol)
   print(vol.attrib)
   print(vol.attrib.get('name'))
   solid = vol.find('solidref')
   print(solid.attrib)
   print(solid.attrib.get('ref'))
   material = vol.find('materialref')
   if material is not None :
      print('material : '+str(material.attrib))
      print(material.attrib.get('ref'))
   for pv in vol.findall('physvol') :
      volref = pv.find('volumeref')
      print(volref.attrib)
      posref = pv.find('positionref')
      if posref is not None :
         print(posref.attrib.get('ref'))
      rotref = pv.find('rotationref')
      if rotref is not None :
         print(rotref.attrib.get('ref'))
else :
   print('Volume : '+volume+' Not found')
#solidList = []
#for s in tree.xpath('//solids/*') :
#    #print(s.tag)
#    if s.tag not in solidList :
#        solidList.append(s.tag)
#l = 'List of Solids in : '+iname
#print(l)
#print('=' * len(l))
#for i in solidList :
#    print(i)        
