import ROOT
import sys

if len(sys.argv)>1:
	wspath = sys.argv[1]	
else:
	wspath ="workdir/CMS/Axial_MonoJ_NLO_Mphi-1000_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph/ws.root"

toypath = "workdir/CMS/Axial_MonoJ_NLO_Mphi-1000_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph/higgsCombineAxial_MonoJ_NLO_Mphi-1000_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraphasimov_bs.FitDiagnostics.mH120.123456.root"


file=ROOT.TFile(wspath)
# file.ls()
toyfile=ROOT.TFile(toypath)
# toyfile.cd("toys")
# toyfile.ls()

toy = toyfile.Get("toys/toy_asimov")
toy.Print("v")

w = file.Get("w")
# w.Print()



ModelConfig = w.obj("ModelConfig")

ModelConfig.Print()

pdf = ModelConfig.GetPdf()
pdf.Print()

obs = ModelConfig.GetObservables()
obs.Print()

# data = ROOT.RooDataSet("data", "", obs);
# data.add(obs);
# nll = pdf.createNLL(data)


# bsToy = ROOT.RooDataSet("data", "", toy);
# data.add(toy);
nll = pdf.createNLL(toy)


ROOT.RooMinuit(nll).migrad() ;

# r=ModelConfig.GetParametersOfInterest()
r=w.var("r")
r.Print("v")

frame1 = r.frame(-2,2,10) ;
nll.plotOn(frame1) ;

frame1.Draw()
# // Create likelihood function
# RooAbsReal* nll = w::model.createNLL(*data,NumCPU(2)) ;

# // Minimize likelihood
# RooMinuit(*nll).migrad() ;

# // Plot likelihood scan in parameter frac
# RooPlot* frame1 = w::frac.frame(Bins(10),Range(0.01,0.95)) ;
# nll->plotOn(frame1,ShiftToZero()) ;

# // Plot the profile likelihood in frac
# RooAbsReal* pll_frac = nll->createProfile(w::frac) ;
# pll_frac->plotOn(frame1,LineColor(kRed)) ;


raw_input("Press a key...")