# python PlotLimits.py Axial expected CMS workdir/CMS/Axial_MonoJ_NLO_Mphi-*/higgsCombineAxial_MonoJ_NLO_Mphi-*.AsymptoticLimits.mH120.root
# python PlotLimits.py Axial expected ATLAS workdir/ATLAS/Axial_MonoJ_NLO_Mphi-*/higgsCombineAxial_MonoJ_NLO_Mphi-*.AsymptoticLimits.mH120.root
# python PlotLimits.py Axial expected combined workdir/combined/Axial_MonoJ_NLO_Mphi-*/higgsCombineAxial_MonoJ_NLO_Mphi-*.AsymptoticLimits.mH120.root

# python PlotLimits.py --coupling Axial --type=observed --datatype=both --dontExcludeFailedFits --LimitsPath=workdir/combined -o combined

import ROOT
import csv
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
import glob

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPalette(1);

def GetLimit(file_,type="expected"):
    rootfile = ROOT.TFile.Open(file_)
    tree = rootfile.Get("limit")
    if type=="expected":
        tree.GetEntry(2)
    elif type=="down1":
        tree.GetEntry(1)
    elif type=="down2":
        tree.GetEntry(0)
    elif type=="up1":
        tree.GetEntry(3)
    elif type=="up2":
        tree.GetEntry(4)
    elif type=="observed":
        tree.GetEntry(5)
    else:
        print "which kind of limit?"
        exit()
    limit = tree.limit
    rootfile.Close()
    return limit

def GetMasses(file_, datatype_):
    string_ = file_
    print string_
    index_l = string_.find("Mphi")
    index_m = string_.find("Mchi")
    index_r = string_.find("gSM")
    mphi = string_[index_l:index_m-1]
    mchi = string_[index_m:index_r-1] 

    if datatype == "none":
        mphi = mphi.replace("Mphi_","")
        mchi = mchi.replace("Mchi_","")
    else:
        mphi = mphi.replace("Mphi-","")
        mchi = mchi.replace("Mchi-","")


    return [float(mphi),float(mchi)]

def GetXS(file_reader,coupling,mphi,mchi,datatype):
    xs = None
    for row in file_reader:
        #print row['sample']
        if not coupling in str(row['sample']):
            continue
        if datatype == "none":            
            if row['sample'].find('Mphi_'+str(int(mphi))+"_")!=-1 and row['sample'].find('Mchi_'+str(int(mchi))+"_")!=-1:
                print "---------------------------------------------------------------------------------------"
                if str(row['xs'])!="":
                    xs = float(row['xs'])
        else:
            if row['sample'].find('Mphi-'+str(int(mphi))+"_")!=-1 and row['sample'].find('Mchi-'+str(int(mchi))+"_")!=-1:
                print "---------------------------------------------------------------------------------------"
                if str(row['xs'])!="":
                    xs = float(row['xs'])            
    return xs
 
def h2graph(in_hist): # could be TH2F, etc.
  temp_can=ROOT.TCanvas()        # This avoids disturbing canvases already in use.
  temp_can.cd()            # Yes, we do need a canvas for this to work, and
  in_hist.Draw("contlist") # yes, we do need to call Update() on it after
  temp_can.Update()        # Draw(). No, ROOT is not very friendly.

  plah = ROOT.gROOT.GetListOfSpecials().FindObject("contours")

  # agr = malloc(plah.GetSize() * sizeof(TGraph *))
  agr=[]
  for i in range(plah.GetSize()):
    list_ = plah.At(i)
    # agr[i] = list->First()
    agr.append(list_.First())
    agr[i].SetLineColor(in_hist.GetLineColor()) # totally optional
  return agr

def set_color_env():
    NRGBs = 5
    # NCont = 255
    NCont = 1

    stops = [ 0.00, 0.34, 0.61, 0.84, 1.00 ]
    red   = [ 0.00, 0.00, 0.87, 1.00, 0.51 ]
    green = [ 0.00, 0.81, 1.00, 0.20, 0.00 ]
    blue  = [ 0.51, 1.00, 0.12, 0.00, 0.00 ]
    stopsArray = array('d', stops)
    redArray = array('d', red)
    greenArray = array('d', green)
    blueArray = array('d', blue)
    ROOT.TColor.CreateGradientColorTable(NRGBs, stopsArray, redArray, greenArray, blueArray, NCont)
    ROOT.gStyle.SetNumberContours(NCont)

# set_color_env()


usage = "usage: %prog [options]"
parser = OptionParser(usage = usage)

parser.add_option("--LimitsPath",
dest="LimitsPath",
help="get limits from this directory, default=workdir/CMS/",
metavar = "LimitsPath",
default = "workdir/CMS/"
)

parser.add_option("--coupling",
dest="coupling",
help="Plot Limits for which kind of coupling?",
metavar = "coupling",
default = "Axial"
)

parser.add_option("--type",
dest="limit_type",
help="expected or observed limit?",
metavar = "limit_type",
default = "expected"
)

parser.add_option("--datatype",
dest="datatype",
help="realData/MCData/both/none (meaning: exp. from MCData, obs. from realData",
metavar = "datatype",
default = "realData"
)

parser.add_option("--dontExcludeFailedFits",
action="store_true",
dest="dontExcludeFailedFits",
default=False,
help="Don't exclude Samples, where Fit to Data failed"
)

parser.add_option("-o",
dest="plotname",
help="additional outputString, default=CMS?",
metavar = "plotname",
default = "CMS"
)

(options, args) = parser.parse_args()


LimitsPath = options.LimitsPath
limit_type = options.limit_type
datatype = options.datatype
coupling = options.coupling
plotname = options.plotname
dontExcludeFailedFits = options.dontExcludeFailedFits

print limit_type
# if (limit_type != "observed" or limit_type != "expected"):
    # parser.error("possible limit_types are 'observed' and 'expected'")

datatypes = [datatype]
if datatype == "both":
    datatypes = ["realData", "MCData"]
elif datatype == "none":
    datatypes = ["dummy"]

files=[[],[]]
excluded=[]
# only get files, where BestFit on data didn't fail
for subdir, dirs, files_ in os.walk(LimitsPath):
    # print subdir
    for i,datatyp in enumerate(datatypes):
        if datatyp in subdir or datatype == "none":
            # print "found one"
            for file in files_:
                # print file
                filepath=os.path.join(subdir, file)
                # print "filepath:", filepath
                if datatype == "none" and (coupling in filepath) and "Limits" in filepath:
                    print "LimitsFile found:", filepath
                    if "Vector_MonoJ_NLO_Mphi_1500_Mchi_600_gSM_0p25_gDM_1p0" in filepath:
                        continue
                    files[i].append(filepath)


                elif file.startswith("fitDiagnostics") and not ("PseudoExperiment" in os.path.join(subdir, file)) and (coupling in os.path.join(subdir, file)) and (("obs" in os.path.join(subdir, file))):
                # if file.startswith("fitDiagnostics") and (coupling in file) and ("obs" in os.path.join(subdir, file)):
                    # print os.path.join(subdir, file)
                    # print file
                    # filepath=os.path.join(subdir, file)
                    # print filepath
                    rootfile = ROOT.TFile(filepath,"OPEN")

                    fit_s = rootfile.Get("fit_s")
                    # print fit_s
                    if not isinstance(fit_s,ROOT.RooFitResult):
                        print rootfile.GetName(),"does not contain RooFitResult"
                        # noFitResult.append(os.path.basename(os.path.abspath(os.path.join(filepath, os.pardir))))
                        # excluded.append(os.path.basename(os.path.abspath(os.path.join(filepath, os.pardir))))
                        excluded.append(os.path.dirname(os.path.dirname(filepath)))
                        print "continued"
                        continue

                    if fit_s.status() == 0 and fit_s.covQual() == 3:
                        #print "loading values"
                        var = fit_s.floatParsFinal().find("r")
                        val = var.getVal()
                        error = var.getError()
                        # print os.path.abspath(os.path.join(filepath, os.pardir))+"/*AsymptoticLimits*.root"
                        if not dontExcludeFailedFits:
                            files[i].append(glob.glob(os.path.abspath(os.path.join(filepath, os.pardir))+"/*AsymptoticLimits*.root")[0])
                    else:
                        print "something went wrong in the fit, file", rootfile.GetName()
                        print "\tfit status =", fit_s.status()
                        print "\tcovQual() =", fit_s.covQual()
                        # excluded.append(os.path.basename(os.path.abspath(os.path.join(filepath, os.pardir))))
                        excluded.append(os.path.dirname(os.path.dirname(filepath)))

                    if dontExcludeFailedFits:
                        files[i].append(glob.glob(os.path.abspath(os.path.join(filepath, os.pardir))+"/*AsymptoticLimits*.root")[0])

                    fit_s.Delete()


# print files
reader = None
csvfile = open('Madgraph_Signal_XS.csv',"r")
reader = csv.DictReader(csvfile)

ROOT.gStyle.SetPadRightMargin(0.15)
w=1200
h=800
canvas = ROOT.TCanvas()
canvas.cd()
canvas.SetCanvasSize(w,h)
#canvas.cd(1)
canvas.SetLogz()

types=["expected", "up1", "down1", "up2", "down2"]
graphlist = [ROOT.TGraph2D(),ROOT.TGraph2D(),ROOT.TGraph2D(),ROOT.TGraph2D(),ROOT.TGraph2D()]
graphlistMCdata = [ROOT.TGraph2D(),ROOT.TGraph2D(),ROOT.TGraph2D(),ROOT.TGraph2D(),ROOT.TGraph2D()]

if limit_type == "observed":
    graphlist.append(ROOT.TGraph2D())
    graphlistMCdata.append(ROOT.TGraph2D())
    types.append("observed")

points=ROOT.TGraph()
x=ROOT.Double(0)
y=ROOT.Double(0)
# real Data
for i in range(len(files[0])): 
    csvfile.seek(0)
    masses = GetMasses(files[0][i], datatype)
    # xs = GetXS(reader,coupling,masses[0],masses[1],datatype)
    # if xs==None:
        # continue
    print "mphi,mchi ",masses
    # print "xs ",xs
    points.SetPoint(i,masses[0], masses[1])
    # points.Print()

    n = points.GetPoint(i,x,y)
    print x,y
    print n
    # if n == 30:
    #     raw_input()

    for j,type_ in enumerate(types):
        limit = GetLimit(files[0][i],type_)
        print type_, ":" 
        # ratio = limit/xs
        print "limit ",limit
        # print "mu/xs ",ratio
        graphlist[j].SetPoint(i,masses[0],masses[1],limit)
        # graphlist[j].SetPoint(i,masses[0],masses[1],ratio)
        if datatype != "none" or coupling == "Vector":
            graphlist[j].SetMarginBinsContent(500)


#MC data
if datatype == "both":
    for i in range(len(files[1])):
        csvfile.seek(0)
        masses = GetMasses(files[1][i], datatype)
        # xs = GetXS(reader,coupling,masses[0],masses[1],datatype)
        # if xs==None:
            # continue
        print "MCData mphi,mchi ",masses
        # print "xs ",xs
        points.SetPoint(i,masses[0], masses[1])

        for j,type_ in enumerate(types):
            limit = GetLimit(files[1][i],type_)
            print type_, ":" 
            # ratio = limit/xs
            print "limit ",limit
            # print "mu/xs ",ratio
            graphlistMCdata[j].SetPoint(i,masses[0],masses[1],limit)
            # graphlistMCdata[j].SetPoint(i,masses[0],masses[1],ratio)

            graphlistMCdata[j].SetMarginBinsContent(500)


graph=graphlist[0]
# print graph
if datatype == "both":
    graph=graphlistMCdata[0] #get median expected from MCdata
# graph=graphlist[len(graphlist)-1]
# graph=graphlist[4]
graph.SetNpx(500)
graph.SetNpy(500)
graph.Print()
graph.GetHistogram().GetXaxis().SetTitle("m_{med} [GeV]")
graph.GetHistogram().GetXaxis().SetTitleSize(0.05)
graph.GetHistogram().GetXaxis().SetTitleOffset(0.8)

graph.GetHistogram().GetYaxis().SetTitle("m_{DM} [GeV]")
graph.GetHistogram().GetYaxis().SetTitleSize(0.05)
graph.GetHistogram().GetYaxis().SetTitleOffset(0.8)

#graph.GetHistogram().SetMaximum(100.)


graph.GetHistogram().SetMinimum(0.01)
graph.GetHistogram().GetZaxis().SetTitle("#sigma_{95% CL}^{expected}/#sigma_{th}")
graph.GetHistogram().GetZaxis().SetTitleSize(0.05)
graph.GetHistogram().GetZaxis().SetTitleOffset(0.8)
graph.SetTitle("Exclusion limits at 95% CL")
# graph.SetMargin(0.)
# graph.SetMarginBinsContent(0.1);
graph.SetMarkerStyle(20)
graph.GetHistogram().Smooth()
graph.GetHistogram().SetContour(1000)
# graph.GetHistogram().Draw("tri1colz")
# graph.GetHistogram().Draw("colz")
# canvas.SetTheta(90)
# canvas.SetPhi(0)
# graph.Draw("tri1colz")
# raw_input()
graph.Draw("colz")

graph.GetHistogram().SetAxisRange(0.02, 50.,"Z");

# graph.Draw("colz2 tri1")

cms = ROOT.TLatex(0.12, 0.91, 'CMS private work in progress'  )
cms.SetNDC()
cms.SetTextSize(0.04)
ROOT.gStyle.SetPalette(1)
lumi = ROOT.TLatex(0.7, 0.91, '35.9 fb^{-1} (13 TeV)'  )
lumi.SetNDC()
lumi.SetTextSize(0.03)


contours = array('d',[1.0])
contgraphs=[]

if datatype=="both" and datatype != "none":
    for k, gr in enumerate(graphlistMCdata):
        graph_copy = gr.Clone()
        graph_copy.GetHistogram().SetContour(1,contours)
        x = h2graph(graph_copy.GetHistogram())
        contgraphs.append(x[0].Clone())
        contgraphs[k].SetName(types[k])  

        # graph_copy.Draw("cont3 same") 
        # print x      
        # x[0].Draw()
        # x[0].Print()
        # graph.Draw("colz same")
        # print contgraphs
        # graph_copy.GetHistogram().SetLineStyle(10)
        # graph_copy.Draw("cont LIST same")
        # ROOT.gPad.Update()
    graph_copy = graphlist[len(graphlist)-1].Clone()
    graph_copy.GetHistogram().SetContour(1,contours)
    # graph_copy.Draw("cont3 same") 

    x = h2graph(graph_copy.GetHistogram())
    contgraphs.append(x[0].Clone()) 
    contgraphs[len(contgraphs)-1].SetName("observed_realData")    
else:    
    for k, gr in enumerate(graphlist):
        graph_copy = gr.Clone()
        graph_copy.GetHistogram().SetContour(1,contours)
        x = h2graph(graph_copy.GetHistogram())
        # print x
        contgraphs.append(x[0].Clone())  
      
        # x[0].Draw()
        # x[0].Print()
        # graph.Draw("colz same")
        # print contgraphs
        # graph_copy.GetHistogram().SetLineStyle(10)
        # graph_copy.Draw("cont LIST same")
        # ROOT.gPad.Update()

print contgraphs

if datatype != "none" or coupling == "Vector":
    for gr in graphlist:
        gr.SetMarginBinsContent(0)
    for gr in graphlistMCdata:
        gr.SetMarginBinsContent(0)
#     print "prior:", gr.GetN()
#     for i in range(gr.GetN()):

#     # for i in range(50, 139):
#         gr.RemovePoint(i)
#     print "after:", gr.GetN()

#     gr.Print()
#     gr.Draw("AC")
#     canvas.Update()
# canvas.Update()

cms.Draw("same")
lumi.Draw("same")
canvas.cd()
cms.Draw("same")
lumi.Draw("same")


contgraphs[0].Draw("Csame")
contgraphs[1].Draw("C")
contgraphs[1].SetLineStyle(5)
contgraphs[2].Draw("C")
contgraphs[2].SetLineStyle(5)
# contgraphs[3].Draw("C")
# contgraphs[3].SetLineStyle(6)
# contgraphs[4].Draw("C")
# contgraphs[4].SetLineStyle(6)


if limit_type == "observed":
    # for i in range(10, contgraphs[len(contgraphs)-1].GetN()+1):
    # # for i in range(50, 139):
    #     contgraphs[len(contgraphs)-1].RemovePoint(i)
    contgraphs[len(contgraphs)-1].Draw("C")
    # contgraphs[len(contgraphs)-1].Print()
    contgraphs[len(contgraphs)-1].SetLineColor(ROOT.kRed)
points.Draw("*")

legend = ROOT.TLegend(0.1,0.7,0.48,0.9);
legend.AddEntry(contgraphs[0],"Median expected","l");
legend.AddEntry(contgraphs[1],"Median expected #pm 1 #sigma","l");
legend.AddEntry(contgraphs[len(contgraphs)-1],"Observed","l");
legend.Draw();
canvas.Print(coupling+"_"+limit_type+"_hist_"+datatype+plotname+".pdf")

file = ROOT.TFile("Limitcontours"+coupling+"_"+limit_type+"_hist_"+datatype+plotname+".root","RECREATE")

for i,graph in enumerate(graphlistMCdata):
    graph.SetName("MCgraph_"+types[i])
    graph.Write()

for i,graph in enumerate(graphlist):
    graph.SetName("graph_"+types[i])
    graph.Write()

for cont in contgraphs:
    cont.Write() 

points.SetName("points")
points.Write()
points.Print()

canvas.Write()

file.Close()
print "excluded following samples:\n" , excluded



