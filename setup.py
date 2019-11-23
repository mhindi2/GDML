from setuptools import setup
#from pygears import __version__
__version__ = '0.1'

import os
if 'posix' in os.name:
    setup(name='FreeCad_GDML_Workbench',
      version=str(__version__),
      packages=['freecad','lxml']
      maintainer="keithsloan52",
      maintainer_email="keith@sloan-home.co.uk",
      url="https://github.com/KeithSloan/FreeCAD_GDML_Workbench",
      description="GDML Workbench for FreeCAD",
      install_requires=['lxml'],
      include_package_data=True
      include_package_data=True)

else:
    setup(name='FreeCad_GDML_Workbench',
      version=str(__version__),
      packages=['freecad','lxml']
      maintainer="keithsloan52",
      maintainer_email="keith@sloan-home.co.uk",
      url="https://github.com/KeithSloan/FreeCAD_GDML_Workbench",
      description="GDML Workbench for FreeCAD",
      install_requires=['python3-lxml'],
      include_package_data=True
      include_package_data=True)

