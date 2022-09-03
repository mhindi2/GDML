#########################################################
#
#   Assembly dictionary
#
#########################################################
# Geant4 generates its own internal names for the
# physical volumes that are included in the assembly.
# The names are of the form
#
# av_WWW_impr_XXX_YYY_ZZZ
#
# where:
# WWW - assembly volume instance number. I gather this is a number that starts
#       at 1 when the first assembly structure is encountered in the gdml file
#       and then is incremented by one every time a new <assembly structure
#       is encounterd
# XXX - assembly volume imprint number. I am not sure what this is. I assume that
#       a given assembly can be placed inside several logical volume mothers
#       and that the imprint number keeps track of how many times the one and same
#       assembly has been placed. For now I am going to take this as 1. I.e, no more
#       than one placement for each assembly
# YYY - the name of the placed logical volume. geant seems to append an '_pv'
#       to out generated logical volume names, which for current exportGDML is
#        'V-' + copied solid name
# ZZZ - the logical volume index inside the assembly volume. This is a counter
#       that increases sequentially with the number of solids in the assembly
# to keep track of what a bordersurface is replaced by we use the following
# dictionaries. Each has as key the name of the App::Part  that appears
# in the PV1, PV2 properties of the bordersurface

from .exportGDML import isAssembly, assemblyHeads


class AssemblyHelper:
    instCount = 0

    def __init__(self, assemblyVol):
        AssemblyHelper.instCount += 1
        self.assemblyVol = assemblyVol
        self.xxx = 1
        self.www = AssemblyHelper.instCount

    def getSubVols(self, vol):
        volsList = []
        for obj in assemblyHeads(vol):
            if isAssembly(obj):
                volsList += self.getSubVols(obj)
            else:
                volsList.append(obj.Name)

        return volsList

