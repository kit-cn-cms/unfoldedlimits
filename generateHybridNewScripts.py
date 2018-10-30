import ROOT
import sys
import os
import csv
from array import array

ROOT.gROOT.SetBatch(True)

mode = sys.argv[1]
if mode == "CMS":
    print "creating HybridNew scripts for CMS run"
elif mode == "ATLAS":
    print "creating HybridNew scripts for ATLAS run"
elif mode == "combined":
    print "creating HybridNew scripts for combined run"
else:
    sys.exit("Only 'CMS', 'ATLAS' or 'combined' valid options!")

script = """#!/bin/bash 

cd /nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/CMSSW_8_1_0/src
eval `scramv1 runtime -sh` 
cd /nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits
mkdir -p workdir
cd workdir\n"""

script += "mkdir -p " + mode + "\n"
script += "cd " + mode + "\n"

axialsamples = [
    # "Axial_MonoJ_NLO_Mphi-500_Mchi-200_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-300_Mchi-100_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-1500_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-1750_Mchi-10_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-1750_Mchi-100_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-300_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-2250_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-2000_Mchi-200_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-500_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-2250_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-100_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-2000_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-1750_Mchi-150_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-2000_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-2000_Mchi-100_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-1500_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-500_Mchi-150_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-1750_Mchi-500_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-750_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-1000_Mchi-600_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-750_Mchi-600_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-300_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-1750_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-100_Mchi-40_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-1000_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-1750_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-2000_Mchi-10_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-300_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-500_Mchi-500_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-1000_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-2250_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-750_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-500_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-2000_Mchi-150_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-100_Mchi-200_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-1750_Mchi-200_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-2000_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-2250_Mchi-150_gSM-0p25_gDM-1p0_13TeV-madgraph",
    # "Axial_MonoJ_NLO_Mphi-1500_Mchi-600_gSM-0p25_gDM-1p0_13TeV-madgraph"
]

vectorsamples = [
    "Vector_MonoJ_NLO_Mphi-1000_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-1000_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-1000_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-1000_Mchi-600_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-100_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-100_Mchi-200_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-100_Mchi-40_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-1500_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-1500_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-1500_Mchi-600_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-1750_Mchi-150_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-1750_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-1750_Mchi-200_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-1750_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-1750_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-1750_Mchi-500_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2000_Mchi-100_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2000_Mchi-10_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2000_Mchi-150_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2000_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2000_Mchi-200_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2000_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2000_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2000_Mchi-500_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2250_Mchi-100_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2250_Mchi-10_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2250_Mchi-150_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2250_Mchi-200_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2250_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2250_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-2250_Mchi-500_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-300_Mchi-100_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-300_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-300_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-300_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-500_Mchi-150_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-500_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-500_Mchi-200_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-500_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Vector_MonoJ_NLO_Mphi-500_Mchi-500_gSM-0p25_gDM-1p0_13TeV-madgraph",
]

pseudosamples = [
    "Pseudo_MonoJ_NLO_Mphi-1000_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-1000_Mchi-350_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-100_Mchi-10_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-100_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-10_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-10_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-200_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-200_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-20_Mchi-0_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-20_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-300_Mchi-100_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-300_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-300_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-350_Mchi-100_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-350_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-350_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-400_Mchi-100_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-400_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-500_Mchi-100_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-500_Mchi-200_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-500_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-50_Mchi-10_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-50_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Pseudo_MonoJ_NLO_Mphi-750_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph"
]

scalarsamples = [
    "Scalar_MonoJ_NLO_Mphi-100_Mchi-40_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-10_Mchi-10_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-10_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-10_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-200_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-200_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-20_Mchi-0_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-20_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-300_Mchi-100_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-300_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-350_Mchi-100_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-350_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-350_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-400_Mchi-100_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-400_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-400_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-500_Mchi-100_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-500_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-500_Mchi-200_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-500_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-50_Mchi-10_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-50_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph",
    "Scalar_MonoJ_NLO_Mphi-750_Mchi-350_gSM-1p0_gDM-1p0_13TeV-madgraph",
]

listnames = ["axial", "vector", "pseudo", "scalar"]
# lists = [axialsamples, vectorsamples, pseudosamples, scalarsamples]
lists = [axialsamples]

samplescript = ""
for samples, listname in zip(lists, listnames):
    for i, sample in enumerate(samples):
        samplescript += "echo '" + sample + "'\n"
        samplescript += "mkdir -p '" + sample + "'\n"
        samplescript += "cd '" + sample + "'\n"
        samplescript += "cd realData \n"
        if mode == "CMS":
            samplescript += """root -l -b -q '/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/prepareUnfoldedWorkspace_v3.C(\"""" + \
                sample + """\")' 2>&1 | tee ws.log\n"""
        elif mode == "ATLAS":
            samplescript += """root -l -b -q '/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/prepareUnfoldedWorkspace_v3ATLAS.C(\"""" + \
                sample + """\")' 2>&1 | tee ws.log\n"""
        elif mode == "combined":
            samplescript += """root -l -b -q '/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/prepareCombinedWorkSpace_v2Test.C(\"""" + \
                sample + """\")' 2>&1 | tee ws.log\n"""
        n = 0      
        HybridScript = ""
        for point in [x * 0.2 for x in range(1, 50)]:
            # print point, n
            HybridScript += """combine -M HybridNew -v 2 -s -1 --LHCmode LHC-limits --singlePoint """ + str(point) + """ --saveToys --saveHybridResult --iterations 10 -T 100 --clsAcc 0 ws.root -D data\n """
            if n>0:
                if n % 10 == 0 or n == len([x * 0.2 for x in range(1, 50)])-1:
                    HybridScript += "cd /nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits  \n"
                    directory = "/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/workdir/"+mode+'/scripts/'
                    filename = directory+"Job_HybridNew_"+listname+"_"+str(i)+"_"+str(n)+'.sh'
                    if not os.path.exists(directory):
                        os.makedirs(directory)
                    f = open(filename, 'w')
                    f.write(script+samplescript+HybridScript)
                    f.close()
                    HybridScript = ""
                    print filename                
            n+=1
        samplescript += "cd .. \n"

print "created scripts in ", directory


# // combineTool.py -M Impacts -d workspaces/workspaceUnfolded.root -D data -m 125 --doInitialFit --robustFit 1 --rMin -2 --rMax 2 -t -1 --expectSignal=0 -n asimov_bonly
# // combineTool.py -M Impacts -d workspaces/workspaceUnfolded.root -D data -m 125 --robustFit 1 --doFits --parallel 5 --rMin -2 --rMax 2 -t -1 --expectSignal=0 -n asimov_bonly
# // combineTool.py -M Impacts -d workspaces/workspaceUnfolded.root -m 125 -o impacts.json --rMin -2 --rMax 2 -t -1 --expectSignal=0 -n asimov_bonly
# // plotImpacts.py -i impacts.json -o impacts
