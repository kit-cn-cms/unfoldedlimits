from array import array
from optparse import OptionParser
from optparse import OptionGroup
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import sys
import os
import stat
import subprocess
import time
import shutil
import imp
from ast import literal_eval

ROOT.gROOT.SetBatch(True)
ROOT.gDirectory.cd('PyROOT:/')

directory = os.path.dirname(os.path.abspath(sys.argv[0]))
basefolder = os.path.abspath(os.path.join(directory, "base"))

if not basefolder in sys.path:
    sys.path.append(basefolder)

pathToCMSSWsetup="/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/setupCMSSW_8_1_0.txt"



usage = "usage: %prog [options]"
parser = OptionParser(usage = usage)

parser.add_option("--mode",
dest="mode",
help="Calculate Limits for which experiment",
metavar = "arg"
)

parser.add_option("--generateOnly",
dest="generateOnly",
action="store_true",
help="only generate Script, don't automatically submit them",
metavar = "generateOnly",
default = False
)

parser.add_option("--doPulls",
dest="doPulls",
action="store_true",
help="include production of PullPlots",
metavar = "doPulls",
default = False
)

parser.add_option("--doImpacts",
dest="doImpacts",
action="store_true",
help="include production of ImpactPlots",
metavar = "doImpacts",
default = False
)

parser.add_option("--useMassPointsList",
dest="useMassPointsList",
action="store_true",
help="get mass points from file specified with '--MassPointsList'",
metavar = "useMassPointsList",
default = False
)


parser.add_option("--MassPointsList",
dest="MassPointsList",
action="store_true",
help="file to provide mass points",
metavar = "MassPointsList",
default = "/nfs/dust/cms/user/swieland/Darkmatter/DM_Unfolding/axial_masspoints.txt"
)

parser.add_option("--SamplePath",
dest="SamplePath",
action="store_true",
help="Path of Samples, default=/nfs/dust/cms/user/swieland/Darkmatter/ntuples_signal_madgraph_tagging",
metavar = "SamplePath",
default = "/nfs/dust/cms/user/swieland/Darkmatter/ntuples_signal_madgraph_tagging"
)

parser.add_option("--specificSample",
dest="specificSample",
help="Name of specific Sample to Run",
metavar = "specificSample",
default = ""
)

parser.add_option("--specificCoupling",
dest="specificCoupling",
help="Name of specific Coupling to Run",
metavar = "specificCoupling",
default = ""
)

(options, args) = parser.parse_args()

if options.mode == None:
    parser.error("type has to be specified!")

mode = options.mode
doPulls = options.doPulls
doImpacts = options.doImpacts
useMassPointsList = options.useMassPointsList
MassPointsList = options.MassPointsList
generateOnly = options.generateOnly
SamplePath = options.SamplePath
specificSample = options.specificSample
specificCoupling = options.specificCoupling

axialsamples = []
vectorsamples = []
pseudosamples = []
scalarsamples = []

listnames = ["Axial", "Vector", "Pseudo", "Scalar"]
lists = [axialsamples, vectorsamples, pseudosamples, scalarsamples]

specificCouplingSampleList = []
if specificCoupling != "":
    listnames = [specificCoupling]
    lists = [specificCouplingSampleList]
    if useMassPointsList == True:
        print "using mass point list"
        l=[]
        with open(MassPointsList) as f:
            for line in f:
                lists[0].append(line.strip('\n'))
            # l = [list(literal_eval(line)) for line in f]
        # print l
    else:
        print "not using mass point list"

if specificSample=="" and not useMassPointsList:
    for x in os.listdir(SamplePath):
        for i,listname in enumerate(listnames):
            if x.startswith(listname):
                lists[i].append(x)
else:
    lists[0].append(specificSample)
# print axialsamples



# commonfitoptions = " -D data  --cminDefaultMinimizerTolerance 0.001 --cminDefaultMinimizerStrategy 0 --rMin=-10 --rMax=10 --saveToys -d ws.root"
# commonImpactfitoptions = " -D data  --cminDefaultMinimizerTolerance 0.001 --cminDefaultMinimizerStrategy 0 --rMin=-10 --rMax=10 -d ws.root"
# # commonFitDiagnosticsOptions=" --robustFit=1 --forceRecreateNLL --saveWithUncertainties --saveNLL "
# commonFitDiagnosticsOptions = " --minos all --forceRecreateNLL --saveWithUncertainties --saveNLL "

commonfitoptions = " -D data  --cminDefaultMinimizerTolerance 0.1 --cminDefaultMinimizerStrategy 0 --saveToys -d ws.root"
commonImpactfitoptions = " -D data  --cminDefaultMinimizerTolerance 0.1 --cminDefaultMinimizerStrategy 0 -d ws.root"
commonFitDiagnosticsOptions=" --robustFit=1 --forceRecreateNLL --saveWithUncertainties --saveNLL"
# commonFitDiagnosticsOptions = " --minos all --forceRecreateNLL --saveWithUncertainties --saveNLL "

AsymptoticLimitsOption = """ --run="both" --noFitAsimov """
# AsymptoticLimitsOption = """ --run="both" """

range_as_bonly = " --rMin=-0.2 --rMax=0.2 "
range_as_bs = " --rMin=0.5 --rMax=1.5 "
range_obs = " --rMin=-3 --rMax=3 "
range_limit = " --rMin=0 --rMax=1 "

# range_10 = " --rMin=-10 --rMax=10 "
# range_as_bonly = range_10
# range_as_bs = range_10
# range_obs = range_10

commands=[]
script = """#!/bin/bash 

cd /nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/CMSSW_8_1_0/src
eval `scramv1 runtime -sh` 
cd /nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits
mkdir -p workdir
cd workdir\n"""

script += "mkdir -p " + mode + "\n"
script += "cd " + mode + "\n"


samplescript = ""
for samples, listname in zip(lists, listnames):
    for i, sample in enumerate(samples):
        for k,string in enumerate(["cd realData\n","cd MCData\n"]):
            wsString = ""
            if k == 0:
                samplescript += "echo '" + sample + "'\n"
                samplescript += "mkdir -p '" + sample + "'\n"
                samplescript += "cd '" + sample + "'\n"
                samplescript += "mkdir -p realData\n"
                samplescript += "mkdir -p MCData\n"
                wsString = ", true"                
            if k == 1: wsString = ", false"

 
            samplescript += string

            if mode == "CMS":
                samplescript += """root -l -b -q '/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/prepareUnfoldedWorkspace_v3.C(\"""" + \
                    sample+"\"" + wsString + """)' 2>&1 | tee ws.log\n"""
            elif mode == "ATLAS":
                samplescript += """root -l -b -q '/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/prepareUnfoldedWorkspace_v3ATLAS.C(\"""" + \
                    sample + """)' 2>&1 | tee ws.log\n"""
            elif mode == "combined":
                samplescript += """root -l -b -q '/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/prepareUnfoldedWorkspace_v3.C(\"""" + \
                    sample+"\"" + wsString + """, true)' 2>&1 | tee ws.log\n"""

            samplescript += """combine -M FitDiagnostics -t -1 """+ range_as_bonly + commonfitoptions + commonFitDiagnosticsOptions + """ --expectSignal=0   -n """ + \
                sample + """asimov_bonly 2>&1 | tee BestFit_bonly.log\n"""
            samplescript += """combine -M FitDiagnostics -t -1  """+ range_as_bs + commonfitoptions + commonFitDiagnosticsOptions + """ --expectSignal=1 -n """ + \
                sample + """asimov_bs 2>&1 | tee BestFit_bs.log\n"""
            samplescript += """combine -M FitDiagnostics """+ range_obs + commonfitoptions + commonFitDiagnosticsOptions + """ -n """ + \
                sample + """obs 2>&1 | tee BestFit.log\n"""
            samplescript += """combine -M AsymptoticLimits """+ range_limit  + commonfitoptions + AsymptoticLimitsOption + """ -n """ + \
                sample + """ | tee limit.log\n"""

            if (doImpacts):
                samplescript += """combineTool.py -M Impacts -m 125 --doInitialFit  """ + range_as_bonly + """ -t -1 --expectSignal=0 --saveNLL  """ + commonImpactfitoptions + """ | tee initialFit_asimov_bonly.log\n"""
                samplescript += """combineTool.py -M Impacts -m 125 --doFits --parallel 2 """ + range_as_bonly + """  -t -1 --expectSignal=0 --saveNLL  """ + commonImpactfitoptions + """ | tee doFits_asimov_bonly.log\n"""
                samplescript += """combineTool.py -M Impacts -m 125 -d ws.root -o impacts_asimov_bonly.json -t -1 --expectSignal=0 | tee createJson_asimov_bonly.log\n"""
                samplescript += """plotImpacts.py -i impacts_asimov_bonly.json -o impacts_b_only"""+sample+""" --per-page 50 | tee plotImpact_asimov_bonly.log\n"""

                samplescript += """combineTool.py -M Impacts -m 125 --doInitialFit """ + range_as_bs + """ -t -1 --expectSignal=1  --saveNLL  """ + commonImpactfitoptions + """ | tee initialFit_asimov_bs.log\n"""
                samplescript += """combineTool.py -M Impacts -m 125 --doFits --parallel 2 """ + range_as_bs + """ -t -1 --expectSignal=1  --saveNLL  """ + commonImpactfitoptions + """ | tee doFits_asimov_bs.log\n"""
                samplescript += """combineTool.py -M Impacts -m 125 -d ws.root -o impacts_asimov_bs.json -t -1 --expectSignal=1 | tee createJson_asimov_bs.log\n"""
                samplescript += """plotImpacts.py -i impacts_asimov_bs.json -o impacts_bs"""+sample+""" --per-page 50 | tee plotImpact_asimov_bs.log\n"""

                samplescript += """combineTool.py -M Impacts -m 125 --doInitialFit --forceRecreateNLL """ + range_obs + """ --saveNLL  """ + commonImpactfitoptions + """ | tee initialFit_obs.log\n"""
                samplescript += """combineTool.py -M Impacts -m 125 --doFits --parallel 2 """ + range_obs + """  --saveNLL  """ + commonImpactfitoptions + """ | tee doFits_obs.log\n"""
                samplescript += """combineTool.py -M Impacts -m 125 -d ws.root -o impacts_obs.json  | tee createJson_obs.log\n"""
                samplescript += """plotImpacts.py -i impacts_obs.json -o impacts_obs"""+sample+""" --per-page 50 | tee plotImpact_obs.log\n"""

            if doPulls:
                samplescript += """python /nfs/dust/cms/user/swieland/Darkmatter/pyroot-plotscripts/tools/diffNuisances.py fitDiagnostics""" + sample + \
                    """asimov_bonly.root -g pulls_bonly"""+sample+""" -a 2>&1 | tee makePullsasimov_bonly.log\n"""
                samplescript += """python /nfs/dust/cms/user/swieland/Darkmatter/pyroot-plotscripts/tools/diffNuisances.py fitDiagnostics""" + sample + \
                    """asimov_bs.root -g pulls_bs"""+sample+""" -a 2>&1 | tee makePullsasimov_bs.log\n""" 
                samplescript += """python /nfs/dust/cms/user/swieland/Darkmatter/pyroot-plotscripts/tools/diffNuisances.py fitDiagnostics""" + sample + \
                    """obs.root -g pulls_obs"""+sample+""" -a 2>&1 | tee makePulls_obs.log\n"""               
            samplescript += "cd .. \n"
        samplescript += "cd .. \n"

            # if i > 0:
                # if i % 10 == 0 or i == len(samples)-1:
        samplescript += "cd /nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits  \n"
        # samplescript += "set +x"
        directory = "/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/workdir/"+mode+'/scripts/'
        if (doImpacts):
            filename = directory+"Impacts_"+sample+'.sh'
        elif (doPulls):
            filename = directory+"pulls_"+sample+'.sh'
        else:
            filename = directory+sample+'.sh'
        if not os.path.exists(directory):
            os.makedirs(directory)
        f = open(filename, 'w')
        f.write(script+samplescript)
        f.close()
        samplescript = ""
        print filename
        commands.append(filename)

if not generateOnly:
    submitstring = "python submit.py -a workdir/LimitArrayJob.sh " 
    for command in commands:
        submitstring += command
        submitstring += " "

    os.system(submitstring)  