import ROOT
import sys
import os
import csv
from array import array

ROOT.gROOT.SetBatch(True)

mode = sys.argv[1]
if mode == "CMS":
    print "creating scripts for CMS run"
elif mode == "ATLAS":
    print "creating scripts for ATLAS run"
elif mode == "combined":
    print "creating scripts for combined run"
else:
    sys.exit("Only 'CMS', 'ATLAS' or 'combined' valid options!")

opt = sys.argv[2:]
doImpacts = False
doPulls = False
print opt
if "doImpacts" in opt:
    doImpacts = True
    print "creating script with impacts production included"
elif "doPulls" in opt:
    doPulls = True
    print "creating script with pulls production included"
else:
    print "not recognizing option ", opt, " ... possible options are 'doPulls' or 'doImpacts' -> creating standard scripts"

script = """#!/bin/bash 

cd /nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/CMSSW_8_1_0/src
eval `scramv1 runtime -sh` 
cd /nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits
mkdir -p workdir
cd workdir\n"""

script += "mkdir -p " + mode + "\n"
script += "cd " + mode + "\n"

axialsamples = [
    "Axial_MonoJ_NLO_Mphi-500_Mchi-200_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-300_Mchi-100_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-1500_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-1750_Mchi-10_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-1750_Mchi-100_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-300_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-2250_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-2000_Mchi-200_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-500_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-2250_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-100_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-2000_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-1750_Mchi-150_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-2000_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-2000_Mchi-100_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-1500_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-500_Mchi-150_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-1750_Mchi-500_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-750_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-1000_Mchi-600_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-750_Mchi-600_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-300_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-1750_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-100_Mchi-40_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-1000_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-1750_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-2000_Mchi-10_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-300_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-500_Mchi-500_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-1000_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-2250_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-750_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-500_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-2000_Mchi-150_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-100_Mchi-200_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-1750_Mchi-200_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-2000_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-2250_Mchi-150_gSM-0p25_gDM-1p0_13TeV-madgraph",
    "Axial_MonoJ_NLO_Mphi-1500_Mchi-600_gSM-0p25_gDM-1p0_13TeV-madgraph"
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
lists = [axialsamples, vectorsamples, pseudosamples, scalarsamples]

samplescript = ""
for samples, listname in zip(lists, listnames):
    for i, sample in enumerate(samples):
        samplescript += "echo '" + sample + "'\n"
        samplescript += "mkdir -p '" + sample + "'\n"
        samplescript += "cd '" + sample + "'\n"
        if mode == "CMS":
            samplescript += """root -l -b -q '/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/prepareUnfoldedWorkspace_v3.C(\"""" + \
                sample + """\")' 2>&1 | tee ws.log\n"""
        elif mode == "ATLAS":
            samplescript += """root -l -b -q '/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/prepareUnfoldedWorkspace_v3ATLAS.C(\"""" + \
                sample + """\")' 2>&1 | tee ws.log\n"""
        elif mode == "combined":
            samplescript += """root -l -b -q '/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/prepareCombinedWorkSpace_v2Test.C(\"""" + \
                sample + """\")' 2>&1 | tee ws.log\n"""

        samplescript += """combine -M FitDiagnostics -t -1 --forceRecreateNLL --expectSignal=0 --saveWithUncertainties --saveNLL --rMin=-10 --rMax=10 --minos all --cminDefaultMinimizerTolerance 0.001 -D data ws.root -n """ + \
            sample + """asimov_bonly 2>&1 | tee BestFit_bonly.log\n"""
        samplescript += """combine -M FitDiagnostics -t -1 --forceRecreateNLL --expectSignal=1 --saveWithUncertainties --saveNLL --rMin=-10 --rMax=10 --minos all --cminDefaultMinimizerTolerance 0.001 -D data ws.root -n """ + \
            sample + """asimov_bs 2>&1 | tee BestFit_bs.log\n"""
        samplescript += """combine -M FitDiagnostics --forceRecreateNLL --saveWithUncertainties --saveNLL --rMin=-10 --rMax=10 --minos all --cminDefaultMinimizerTolerance 0.001 -D data ws.root -n """ + \
            sample + """obs 2>&1 | tee BestFit_bs.log\n"""
        # samplescript += """combine -M AsymptoticLimits --minosAlgo=stepping --rMin=-10 --rMax=10 --minos all --cminDefaultMinimizerTolerance 0.001 --run="both" ws.root -v 1 -D data -n """ + \
        samplescript += """combine -M AsymptoticLimits --minosAlgo=stepping --rMin=-10 --rMax=10 --cminDefaultMinimizerTolerance 0.001 --run="both" --noFitAsimov ws.root -v 1 -D data -n """ + \
            sample + """ | tee limit.log\n"""

        if (doImpacts):
            samplescript += """combineTool.py -M Impacts -d ws.root -D data -m 125 --doInitialFit  -t -1 --forceRecreateNLL --expectSignal=0 --saveNLL --rMin=-10 --rMax=10 --cminDefaultMinimizerTolerance 0.001 | tee initialFit_asimov_bonly.log\n"""
            samplescript += """combineTool.py -M Impacts -d ws.root -D data -m 125 --doFits --parallel 2  -t -1 --forceRecreateNLL --expectSignal=0 --saveNLL --rMin=-10 --rMax=10  --cminDefaultMinimizerTolerance 0.001 | tee doFits_asimov_bonly.log\n"""
            samplescript += """combineTool.py -M Impacts -d ws.root -D data -m 125 -o impacts_asimov_bonly.json --rMin -2 --rMax 2 -t -1 --expectSignal=0 | tee createJson_asimov_bonly.log\n"""
            samplescript += """plotImpacts.py -i impacts_asimov_bonly.json -o impacts_b_only  | tee plotImpact_asimov_bonly.log\n"""

            samplescript += """combineTool.py -M Impacts -d ws.root -D data -m 125 --doInitialFit  -t -1 --forceRecreateNLL --expectSignal=1  --saveNLL --rMin=-10 --rMax=10  --cminDefaultMinimizerTolerance 0.001 | tee initialFit_asimov_bs.log\n"""
            samplescript += """combineTool.py -M Impacts -d ws.root -D data -m 125 --doFits --parallel 2  -t -1 --forceRecreateNLL --expectSignal=1  --saveNLL --rMin=-10 --rMax=10 --cminDefaultMinimizerTolerance 0.001 | tee doFits_asimov_bs.log\n"""
            samplescript += """combineTool.py -M Impacts -d ws.root -D data -m 125 -o impacts_asimov_bs.json --rMin -2 --rMax 2 -t -1 --expectSignal=1 | tee createJson_asimov_bs.log\n"""
            samplescript += """plotImpacts.py -i impacts_asimov_bs.json -o impacts_bs | tee plotImpact_asimov_bs.log\n"""

        if doPulls:
            samplescript += """python /nfs/dust/cms/user/swieland/Darkmatter/pyroot-plotscripts/tools/diffNuisances.py fitDiagnostics""" + sample + \
                """asimov_bonly.root -g pulls_bonly -a 2>&1 | tee makePullsasimov_bonly.log\n"""
            samplescript += """python /nfs/dust/cms/user/swieland/Darkmatter/pyroot-plotscripts/tools/diffNuisances.py fitDiagnostics""" + sample + \
                """asimov_bs.root -g pulls_bs -a 2>&1 | tee makePullsasimov_bs.log\n"""                
        samplescript += "cd .. \n"

        # if i > 0:
            # if i % 10 == 0 or i == len(samples)-1:
        samplescript += "cd /nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits  \n"
        # samplescript += "set +x"
        directory = "/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/workdir/"+mode+'/scripts/'
        if (doImpacts):
            filename = directory+"Impacts_Job_" + \
                listname+"_"+str(i)+'.sh'
        elif (doPulls):
            filename = directory+"pulls_Job_"+listname+"_"+str(i)+'.sh'
        else:
            filename = directory+"Job_"+listname+"_"+str(i)+'.sh'
        if not os.path.exists(directory):
            os.makedirs(directory)
        f = open(filename, 'w')
        f.write(script+samplescript)
        f.close()
        samplescript = ""
        print filename


print "created scripts in ", directory


# // combineTool.py -M Impacts -d workspaces/workspaceUnfolded.root -D data -m 125 --doInitialFit --robustFit 1 --rMin -2 --rMax 2 -t -1 --expectSignal=0 -n asimov_bonly
# // combineTool.py -M Impacts -d workspaces/workspaceUnfolded.root -D data -m 125 --robustFit 1 --doFits --parallel 5 --rMin -2 --rMax 2 -t -1 --expectSignal=0 -n asimov_bonly
# // combineTool.py -M Impacts -d workspaces/workspaceUnfolded.root -m 125 -o impacts.json --rMin -2 --rMax 2 -t -1 --expectSignal=0 -n asimov_bonly
# // plotImpacts.py -i impacts.json -o impacts
