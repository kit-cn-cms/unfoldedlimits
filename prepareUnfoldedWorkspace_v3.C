// StandardHypoTestInvDemo("workdir/CMS/Axial_MonoJ_NLO_Mphi-1000_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph/ws.root","w","ModelConfig","ModelConfig_bonly","data",3, 3, 1, 50, 0, 4, 100, 0, 0)
// combine -M FitDiagnostics workspaces/workspaceUnfolded.root --forceRecreateNLL -D data
// combine -M FitDiagnostics workspaces/workspaceUnfolded.root --forceRecreateNLL -D data -t -1 --expectSignal=0,1 --rMin -2 --rMax 2

// combineTool.py -M Impacts -d workspaces/workspaceUnfolded.root -D data -m 125 --doInitialFit --robustFit 1 --rMin -2 --rMax 2 -t -1 --expectSignal=0 -n asimov_bonly
// combineTool.py -M Impacts -d workspaces/workspaceUnfolded.root -D data -m 125 --robustFit 1 --doFits --parallel 5 --rMin -2 --rMax 2 -t -1 --expectSignal=0 -n asimov_bonly
// combineTool.py -M Impacts -d workspaces/workspaceUnfolded.root -m 125 -o impacts.json --rMin -2 --rMax 2 -t -1 --expectSignal=0 -n asimov_bonly
// plotImpacts.py -i impacts.json -o impacts

// combine -M AsymptoticLimits workspaces/workspaceUnfolded.root -D data -t -1 --expectSignal=1 --rMin -2 --rMax 2
#include <vector>
#include <map>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <cstdlib>  //getenv
#include <string>   //stod

#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TString.h"

#ifndef __CINT__
// you need to include this for compiled macro.
// But for CINT, it needs to be in this ifndef/endif condition
#include "RooGlobalFunc.h"
#endif

#include "RooStats/RooStatsUtils.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/SamplingDistPlot.h"

#include "HiggsAnalysis/CombinedLimit/interface/SimpleGaussianConstraint.h"


using namespace RooFit;
using namespace RooStats;
using namespace std;

TFile* ReadTFile(TString path, TString errorMessage = "Could not read file")
{
  TFile* tf = new TFile(path, "READ");
  if (!tf || tf->IsZombie())
  {
    cerr << errorMessage.Data()  << endl;
    exit(EXIT_FAILURE);
  }

  return tf;
}

TH1D* ReadTH1F(TFile* tf, TString saveName, TString newHistName, TString errorMessage = "Could not read histogram")
{
  TH1D* h1 = 0;

  TObject* obj = 0;
  obj = tf->Get(saveName);
  if (!obj || !obj->InheritsFrom(TH1::Class())) {
    cerr << errorMessage.Data()  << endl;
    exit(EXIT_FAILURE);
  }

  h1 = (TH1D*)obj->Clone(newHistName);
  h1->SetDirectory(0);
  delete obj;
  return h1;
}

struct Sample
{
  TString name;
  TString histName = "";
  std::vector<std::tuple<TString, TH1*>> shapeSysts;
  vector<std::tuple<TString, float>> LNSysts;
  TH1* h_nominal;
  bool signal = false;
} unfolded, unfoldedATLAS, sig, z_nunu_jets, w_lnu_jets, Diboson, z_ll_jets, ttbar, qcd, singletop, gamma_jets;


//unfolded
std::vector<std::tuple<TString, TH1*>>unfShapeSysts = {
  std::make_tuple("CMS_scale_j", nullptr),
  std::make_tuple("CMS_res_j", nullptr),

  // std::make_tuple("DataStat", nullptr),
  // std::make_tuple("fakeStat", nullptr),
  // std::make_tuple("fakeScale", nullptr),
  // std::make_tuple("MCStat", nullptr),

  std::make_tuple("CMS_btag_lf", nullptr),
  std::make_tuple("CMS_btag_hf", nullptr),
  std::make_tuple("CMS_btag_hfstats1", nullptr),
  std::make_tuple("CMS_btag_lfstats1", nullptr),
  std::make_tuple("CMS_btag_hfstats2", nullptr),
  std::make_tuple("CMS_btag_lfstats2", nullptr),
  std::make_tuple("CMS_btag_cferr1", nullptr),
  std::make_tuple("CMS_btag_cferr2", nullptr),

  std::make_tuple("Weight_PU", nullptr),
};

vector<std::tuple<TString, float>> unfLNSysts{ // not working in current context
  // std::make_tuple("Lumi", 1.025),
};

//unfolded ATLAS
std::vector<std::tuple<TString, TH1*>>unfShapeSystsATLAS = {
  std::make_tuple("ATLASCMS_scale_j", nullptr),
  std::make_tuple("ATLASCMS_res_j", nullptr),

  std::make_tuple("ATLASCMS_btag_lf", nullptr),
  std::make_tuple("ATLASCMS_btag_hf", nullptr),
  std::make_tuple("ATLASCMS_btag_hfstats1", nullptr),
  std::make_tuple("ATLASCMS_btag_lfstats1", nullptr),
  std::make_tuple("ATLASCMS_btag_hfstats2", nullptr),
  std::make_tuple("ATLASCMS_btag_lfstats2", nullptr),
  std::make_tuple("ATLASCMS_btag_cferr1", nullptr),
  std::make_tuple("ATLASCMS_btag_cferr2", nullptr),

  std::make_tuple("ATLASWeight_PU", nullptr),
};

vector<std::tuple<TString, float>> unfLNSystsATLAS{
  // std::make_tuple("Lumi", 1.025),
};

//signal
std::vector<std::tuple<TString, TH1*>> sigShapeSysts = {
  std::make_tuple("Weight_scale_variation_muR", nullptr),
  std::make_tuple("Weight_scale_variation_muF", nullptr),

};

vector<std::tuple<TString, float>> sigLNSysts{
  std::make_tuple("XS_signal", 1.05),
};

// z_nunu_jets
std::vector<std::tuple<TString, TH1*>> z_nunu_ShapeSysts = {
  std::make_tuple("Weight_scale_variation_muR", nullptr),
  std::make_tuple("Weight_scale_variation_muF", nullptr),
  std::make_tuple("Weight_PDF", nullptr),


  std::make_tuple("BosonWeight_QCD1", nullptr),
  std::make_tuple("BosonWeight_QCD2", nullptr),
  std::make_tuple("BosonWeight_QCD3", nullptr),

  std::make_tuple("BosonWeight_EW1", nullptr),
  std::make_tuple("ZvvBosonWeight_EW2", nullptr),
  std::make_tuple("ZvvBosonWeight_EW3", nullptr),
  std::make_tuple("ZvvBosonWeight_Mixed", nullptr),

  // std::make_tuple("BosonWeight_EW2", nullptr),
  // std::make_tuple("BosonWeight_EW3", nullptr),
  // std::make_tuple("BosonWeight_Mixed", nullptr),
  std::make_tuple("BosonWeight_Alpha", nullptr),

};

vector<std::tuple<TString, float>> z_nunu_LNSysts{
  // std::make_tuple("XS_z_nunu_jets", 1.025),
};

// w_lnu_jets
std::vector<std::tuple<TString, TH1*>> w_lnu_ShapeSysts{
  std::make_tuple("Weight_scale_variation_muR", nullptr),
  std::make_tuple("Weight_scale_variation_muF", nullptr),
  std::make_tuple("Weight_PDF", nullptr),

  std::make_tuple("BosonWeight_QCD1", nullptr),
  std::make_tuple("BosonWeight_QCD2", nullptr),
  std::make_tuple("BosonWeight_QCD3", nullptr),

  std::make_tuple("BosonWeight_EW1", nullptr),
  std::make_tuple("WlnuBosonWeight_EW2", nullptr),
  std::make_tuple("WlnuBosonWeight_EW3", nullptr),
  std::make_tuple("WlnuBosonWeight_Mixed", nullptr),

  // std::make_tuple("BosonWeight_EW2", nullptr),
  // std::make_tuple("BosonWeight_EW3", nullptr),
  // std::make_tuple("BosonWeight_Mixed", nullptr),
  std::make_tuple("BosonWeight_Alpha", nullptr),
};

vector<std::tuple<TString, float>> w_lnu_LNSysts{
  // std::make_tuple("XS_w_lnu_jets", 1.025),
};

// Diboson
std::vector<std::tuple<TString, TH1*>> Diboson_ShapeSysts{
  // std::make_tuple("Weight_scale_variation_muR", nullptr),
  // std::make_tuple("Weight_scale_variation_muF", nullptr),
  // std::make_tuple("Weight_PDF", nullptr),
};

vector<std::tuple<TString, float>> Diboson_LNSysts{
  std::make_tuple("XS_Diboson", 1.2),
};

// z_ll_jets
std::vector<std::tuple<TString, TH1*>> z_ll_jets_ShapeSysts{
  std::make_tuple("Weight_scale_variation_muR", nullptr),
  std::make_tuple("Weight_scale_variation_muF", nullptr),
  std::make_tuple("Weight_PDF", nullptr),

  std::make_tuple("BosonWeight_QCD1", nullptr),
  std::make_tuple("BosonWeight_QCD2", nullptr),
  std::make_tuple("BosonWeight_QCD3", nullptr),

  std::make_tuple("BosonWeight_EW1", nullptr),
  std::make_tuple("ZllBosonWeight_EW2", nullptr),
  std::make_tuple("ZllBosonWeight_EW3", nullptr),
  std::make_tuple("ZllBosonWeight_Mixed", nullptr),

  // std::make_tuple("BosonWeight_EW2", nullptr),
  // std::make_tuple("BosonWeight_EW3", nullptr),
  // std::make_tuple("BosonWeight_Mixed", nullptr),
  std::make_tuple("BosonWeight_Alpha", nullptr),

};

vector<std::tuple<TString, float>> z_ll_jets_LNSysts{
  // std::make_tuple("XS_z_ll_jets", 1.025),
};

// ttbar
std::vector<std::tuple<TString, TH1*>> ttbar_ShapeSysts{
  std::make_tuple("Weight_scale_variation_muR", nullptr),
  std::make_tuple("Weight_scale_variation_muF", nullptr),
  // std::make_tuple("Weight_PDF", nullptr),
};

vector<std::tuple<TString, float>> ttbar_LNSysts{
  std::make_tuple("XS_ttbar", 1.10),
  // std::make_tuple("XS_Top", 1.10),
};

// singletop
std::vector<std::tuple<TString, TH1*>> singletop_ShapeSysts{
  std::make_tuple("Weight_scale_variation_muR", nullptr),
  std::make_tuple("Weight_scale_variation_muF", nullptr),
  // std::make_tuple("Weight_PDF", nullptr),
};

vector<std::tuple<TString, float>> singletop_LNSysts{
  std::make_tuple("XS_singleTop", 1.10),
  // std::make_tuple("XS_Top", 1.10),

};

// qcd
std::vector<std::tuple<TString, TH1*>> qcd_ShapeSysts{
  std::make_tuple("Weight_scale_variation_muR", nullptr),
  std::make_tuple("Weight_scale_variation_muF", nullptr),
  // std::make_tuple("Weight_PDF", nullptr),
};

vector<std::tuple<TString, float>> qcd_LNSysts{
  std::make_tuple("XS_qcd", 2.0),
};

// gamma_jets
std::vector<std::tuple<TString, TH1*>> gamma_jets_ShapeSysts{
  std::make_tuple("Weight_scale_variation_muR", nullptr),
  std::make_tuple("Weight_scale_variation_muF", nullptr),
  // std::make_tuple("Weight_PDF", nullptr),
};

vector<std::tuple<TString, float>> gamma_jets_LNSysts{
  std::make_tuple("XS_gamma_jets", 1.20),
};

template<class C, class T>
auto contains(const C& v, const T& x)
-> decltype(end(v), true)
{
  return end(v) != std::find(begin(v), end(v), x);
}

// Vector_MonoJ_NLO_Mphi-300_Mchi-100_gSM-0p25_gDM-1p0_13TeV-madgraph
// Axial_MonoJ_NLO_Mphi-1000_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph
// AxialMediator100DM100
void prepareUnfoldedWorkspace_v3(const char* signalname = "Axial_MonoJ_NLO_Mphi-1000_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph", bool realdata = true, bool doCombination = false)
{
  // doCombination = true;
  bool doEigenStuff = false;
  fprintf(stdout, "Starting execution\n");
  TFile*  file = nullptr;
  TFile*  ATLASfile = nullptr;
  if (realdata) {
    // file  = ReadTFile("/nfs/dust/cms/user/swieland/Darkmatter/DM_Unfolding/rootfiles/scaledData.root", "Could not read input file");
    file  = ReadTFile("/nfs/dust/cms/user/swieland/Darkmatter/DM_Unfolding/rootfiles/data.root", "Could not read input file");
    ATLASfile  = ReadTFile("/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/ATLAS_CMSCopy.root", "Could not read input file");
  }
  else {
    file  = ReadTFile("/nfs/dust/cms/user/swieland/Darkmatter/DM_Unfolding/rootfiles/MCdata.root", "Could not read input file");
    ATLASfile  = ReadTFile("/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/ATLAS_CMSCopyMCData.root", "Could not read input file");
  }
  TFile*  Signalfile  = ReadTFile("/nfs/dust/cms/user/swieland/Darkmatter/DM_Unfolding/rootfiles/signals.root", "Could not read input file");
//   TFile*  Signalfile  = ReadTFile("/nfs/dust/cms/user/swieland/Darkmatter/DM_Unfolding/rootfiles/AxialSignals.root", "Could not read input file");
  TFile* MCDatafile  = ReadTFile("/nfs/dust/cms/user/swieland/Darkmatter/DM_Unfolding/rootfiles/MCdata.root", "Could not read input file");

  file  = ReadTFile("/nfs/dust/cms/user/swieland/Darkmatter/DM_Unfolding/rootfiles/data.root", "Could not read input file");


  bool isSplitSample = false;
  // TFile*  file  = ReadTFile("Split.root", "Could not read input file");
  // Get Covariance Matrix
  TH2D* h2_cov = nullptr;
  // file->GetObject("ErrorMatrixTotal", h2_cov);
  // file->GetObject("ErrorMatrix_DataMCStat", h2_cov);
  // MCDatafile->GetObject("ErrorMatrix_DataMCStat", h2_cov);
  // file->GetObject("ErrorMatrix_DataMCFakeStat", h2_cov);
  if (realdata) file->GetObject("ErrorMatrix_DataMCStat", h2_cov);
  else MCDatafile->GetObject("ErrorMatrix_DataMCStat", h2_cov);

  if (h2_cov == nullptr) {
    cerr << "Failed to retrieve Monojet covariance matrix" << endl;
    return;
  }
  TH2* h2_covATLAS = nullptr;
  if (doCombination) h2_covATLAS =  (TH2*) ATLASfile->Get("ErrorMatrix_DataMCStat");

  int Nbins = h2_cov->GetNbinsX();
  // int Nbins = 9;
  // int Nbins = 5;
  int NbinMin = 1;

  // int NbinMin = 6;
  int NBinsCov = Nbins - NbinMin + 1;
  TMatrixDSym cov(NBinsCov);

  std::vector<Sample> Samples;
  std::vector<TString> includeSystematics; // keep track of all used systematics
  // Declare Samples
  unfolded.name = "unfolded";
  unfolded.histName = "unfolded_Gen_Hadr_Recoil_Pt";
  // add stat unc, for each bin
  // for (int iBin = NbinMin; iBin <= Nbins; iBin++) {
  //   unfShapeSysts.push_back(std::make_tuple(TString::Format("Input_Bin%i", iBin), nullptr));
  // }
  unfolded.shapeSysts = unfShapeSysts;
  unfolded.LNSysts = unfLNSysts;
  Samples.push_back(unfolded);
  for (auto &shapesys : unfShapeSysts) includeSystematics.push_back(get<0>(shapesys));

  for (auto &lnsys : unfLNSysts) includeSystematics.push_back(get<0>(lnsys));

  sig.name = "signal";
  sig.histName = signalname;
  sig.shapeSysts = sigShapeSysts;
  sig.LNSysts = sigLNSysts;
  sig.signal = true;
  Samples.push_back(sig);
  for (auto &shapesys : sigShapeSysts) includeSystematics.push_back(get<0>(shapesys));
  for (auto &lnsys : sigLNSysts) includeSystematics.push_back(get<0>(lnsys));

  if (doCombination) {
    unfoldedATLAS.name = "unfoldedATLAS";
    unfoldedATLAS.histName = "unfoldednotATLAS_Gen_Hadr_Recoil_Pt";
    unfoldedATLAS.shapeSysts = unfShapeSystsATLAS;
    unfoldedATLAS.LNSysts = unfLNSystsATLAS;
    Samples.push_back(unfoldedATLAS);
    for (auto &shapesys : unfShapeSystsATLAS) includeSystematics.push_back(get<0>(shapesys));
    for (auto &lnsys : unfLNSystsATLAS) includeSystematics.push_back(get<0>(lnsys));
  }

  z_nunu_jets.name = "z_nunu_jets";
  z_nunu_jets.histName = "z_nunu_jets_Gen_Hadr_Recoil_Pt";
  z_nunu_jets.shapeSysts = z_nunu_ShapeSysts;
  z_nunu_jets.LNSysts = z_nunu_LNSysts;
  Samples.push_back(z_nunu_jets);
  for (auto &shapesys : z_nunu_ShapeSysts) includeSystematics.push_back(get<0>(shapesys));
  for (auto &lnsys : z_nunu_LNSysts) includeSystematics.push_back(get<0>(lnsys));

  w_lnu_jets.name = "w_lnu_jets";
  w_lnu_jets.histName = "w_lnu_jets_Gen_Hadr_Recoil_Pt";
  w_lnu_jets.shapeSysts = w_lnu_ShapeSysts;
  w_lnu_jets.LNSysts = w_lnu_LNSysts;
  Samples.push_back(w_lnu_jets);
  for (auto &shapesys : w_lnu_ShapeSysts) includeSystematics.push_back(get<0>(shapesys));
  for (auto &lnsys : w_lnu_LNSysts) includeSystematics.push_back(get<0>(lnsys));

  Diboson.name = "diboson";
  Diboson.histName = "diboson_Gen_Hadr_Recoil_Pt";
  Diboson.shapeSysts = Diboson_ShapeSysts;
  Diboson.LNSysts = Diboson_LNSysts;
  Samples.push_back(Diboson);
  for (auto &shapesys : Diboson_ShapeSysts) includeSystematics.push_back(get<0>(shapesys));
  for (auto &lnsys : Diboson_LNSysts) includeSystematics.push_back(get<0>(lnsys));

  z_ll_jets.name = "z_ll_jets";
  z_ll_jets.histName = "z_ll_jets_Gen_Hadr_Recoil_Pt";
  z_ll_jets.shapeSysts = z_ll_jets_ShapeSysts;
  z_ll_jets.LNSysts = z_ll_jets_LNSysts;
  Samples.push_back(z_ll_jets);
  for (auto &shapesys : z_ll_jets_ShapeSysts) includeSystematics.push_back(get<0>(shapesys));
  for (auto &lnsys : z_ll_jets_LNSysts) includeSystematics.push_back(get<0>(lnsys));

  ttbar.name = "ttbar";
  ttbar.histName = "ttbar_Gen_Hadr_Recoil_Pt";
  ttbar.shapeSysts = ttbar_ShapeSysts;
  ttbar.LNSysts = ttbar_LNSysts;
  Samples.push_back(ttbar);
  for (auto &shapesys : ttbar_ShapeSysts) includeSystematics.push_back(get<0>(shapesys));
  for (auto &lnsys : ttbar_LNSysts) includeSystematics.push_back(get<0>(lnsys));

  singletop.name = "singletop";
  singletop.histName = "singletop_Gen_Hadr_Recoil_Pt";
  singletop.shapeSysts = singletop_ShapeSysts;
  singletop.LNSysts = singletop_LNSysts;
  Samples.push_back(singletop);
  for (auto &shapesys : singletop_ShapeSysts) includeSystematics.push_back(get<0>(shapesys));
  for (auto &lnsys : singletop_LNSysts) includeSystematics.push_back(get<0>(lnsys));

  qcd.name = "qcd";
  qcd.histName = "qcd_Gen_Hadr_Recoil_Pt";
  qcd.shapeSysts = qcd_ShapeSysts;
  qcd.LNSysts = qcd_LNSysts;
  Samples.push_back(qcd);
  for (auto &shapesys : qcd_ShapeSysts) includeSystematics.push_back(get<0>(shapesys));
  for (auto &lnsys : qcd_LNSysts) includeSystematics.push_back(get<0>(lnsys));

  gamma_jets.name = "gamma_jets";
  gamma_jets.histName = "gamma_jets_Gen_Hadr_Recoil_Pt";
  gamma_jets.shapeSysts = gamma_jets_ShapeSysts;
  gamma_jets.LNSysts = gamma_jets_LNSysts;
  Samples.push_back(gamma_jets);
  for (auto &shapesys : gamma_jets_ShapeSysts) includeSystematics.push_back(get<0>(shapesys));
  for (auto &lnsys : gamma_jets_LNSysts) includeSystematics.push_back(get<0>(lnsys));

  TString fileName = nullptr;
  if (!doCombination) {
    if (realdata) fileName = "/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/workdir/CMS/" + sig.histName + "/realData/ws.root";
    else fileName = "/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/workdir/CMS/" + sig.histName + "/MCData/ws.root";
  }
  else {
    if (realdata) fileName = "/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/workdir/combined/" + sig.histName + "/realData/ws.root";
    else fileName = "/nfs/dust/cms/user/swieland/Darkmatter/CombineStuff/unfoldedlimits/workdir/combined/" + sig.histName + "/MCData/ws.root";
  }
  RooWorkspace w("w");
  w.factory("r[0,-10,10]");
  w.var("r")->setConstant(false);

  RooRealVar* poi = (RooRealVar*) w.var("r");
  cout << poi->getVal() << endl;
  RooArgSet  globalObservables;
  RooArgSet  nuisances;
  RooArgSet  observables;
  RooArgSet  observables_eigen;
  RooArgSet  muVec;


  //TOCHECK: Mean = nominal nuicance?
  //define all nuisance parameters
  //shape Systs -> Gaussian
  for (auto Sample : Samples) {
    cout << "looking at Sample " << Sample.name << endl;
    for (auto& shapetuple : Sample.shapeSysts) {
      TString shape = get<0>(shapetuple);
      if (contains(includeSystematics, shape)) {
        printf("adding %s\n", shape.Data());
        w.factory("NUI_" + shape + "_In[0,-7,7]");
        w.factory("NUI_" + shape + "[0,-7,7]");
        // two different instances of NUI:
        // First will be set constant and added as global observable
        // Second one will be allowed to float and is used to vary n_obs
        // w.factory(TString::Format("SimpleGaussianConstraint:%s_Pdf(%s,%s,1)", ("NUI_" + shape).Data(), ("NUI_" + shape ).Data(), ("NUI_" + shape  + "_In").Data()));
        w.factory(TString::Format("Gaussian:%s_Pdf(%s,%s,1)", ("NUI_" + shape).Data(), ("NUI_" + shape ).Data(), ("NUI_" + shape  + "_In").Data()));
        w.var("NUI_" + shape + "_In")->setConstant(true);
        w.var("NUI_" + shape)->setError(1);
        w.var("NUI_" + shape)->setAsymError(-1, 1);
        globalObservables.add( *w.var(("NUI_" + shape + "_In")) );
        nuisances.add( *w.var(("NUI_" + shape)) );
      }
    }
    //lognormal Systs
    for (auto& lnsys : Sample.LNSysts) {
      TString name = std::get<0>(lnsys);
      float kappa = std::get<1>(lnsys);
      if (contains(includeSystematics, name)) {
        printf("adding %s\n", name.Data());
        w.factory("NUI_" + name + "_In[0,-7,7]");
        w.factory("NUI_" + name + "[0,-7,7]");
        w.factory(name + TString::Format("_kappa[%f]", kappa));
        // two different instances of NUI:
        // First will be set constant and added as global observable
        // Second one will be allowed to float and is used to vary n_obs
        w.factory( TString::Format("expr::alpha_%s('pow(%s,%s)',%s, %s)",
                                   ("NUI_" + name).Data(), (name + "_kappa").Data(), ("NUI_" + name).Data(), (name + "_kappa").Data(), ("NUI_" + name).Data()));
        // w.factory(TString::Format("SimpleGaussianConstraint:%s_Pdf(%s,%s,1)", ("NUI_" + name).Data(), ("NUI_" + name ).Data(), ("NUI_" + name + "_In").Data()));
        // w.factory(TString::Format("SimpleGaussianConstraint:%s_Pdf(%s,%s,1)", ("NUI_" + name).Data(), ("NUI_" + name ).Data(), ("NUI_" + name + "_In").Data()));
        w.factory(TString::Format("Gaussian:%s_Pdf(%s,%s,1)", ("NUI_" + name).Data(), ("NUI_" + name ).Data(), ("NUI_" + name + "_In").Data()));
        w.var("NUI_" + name + "_In")->setConstant(true);
        w.var(name + "_kappa")->setConstant(true);
        w.var("NUI_" + name)->setError(1);
        w.var("NUI_" + name)->setAsymError(-1, 1);
        globalObservables.add(*w.var(("NUI_" + name + "_In")));
        nuisances.add(*w.var(("NUI_" + name)));
      }
    }
  }

  // multiply all nuicance PDFs
  TString string_nui = "PROD:nuisancePdf(";
  for (TString& name : includeSystematics) {
    string_nui += TString::Format("%s_Pdf,", ("NUI_" + name).Data());
  }
  string_nui += ")";
  w.factory(string_nui);
  cout << "Created Nuicance PDF" << endl;



  // bool isSplitSample = true;

  for (int i = 0; i < NBinsCov; i++) {
    for (int j = i; j < NBinsCov; j++) {
      if (j != i) {
        if (doCombination) {
          cov(i, j) = h2_cov->GetBinContent(NbinMin + i  , NbinMin + j) + h2_covATLAS->GetBinContent(NbinMin + i  , NbinMin + j);
          cov(j, i) = h2_cov->GetBinContent(NbinMin + i  , NbinMin + j) + h2_covATLAS->GetBinContent(NbinMin + i  , NbinMin + j);
        }
        else {
          cov(i, j) = h2_cov->GetBinContent(NbinMin + i  , NbinMin + j) ;
          cov(j, i) = h2_cov->GetBinContent(NbinMin + i  , NbinMin + j) ;
          // cov(i, j) = 0;
          // cov(j, i) = 0;
        }
      }
      else
      {
        if (doCombination) cov(i, i) = h2_cov->GetBinContent(NbinMin + i  , NbinMin + j) + h2_covATLAS->GetBinContent(NbinMin + i  , NbinMin + j);
        else cov(i, i) = h2_cov->GetBinContent(NbinMin + i  , NbinMin + j);
      }
    }
  }

  // read Histos and calculate sigma
  for (auto& Sample : Samples) {
    TString name = Sample.name;
    TString histname = Sample.histName;
    cout << "getting " << name << endl;
    // get nominal template
    TH1*  h_nom = nullptr;
    if (!Sample.signal && !Sample.name.Contains("ATLAS")) file->GetObject(histname , h_nom);
    else if (!Sample.signal && Sample.name.Contains("ATLAS")) ATLASfile->GetObject(histname, h_nom);
    else if (Sample.signal) Signalfile->GetObject(histname, h_nom);

    //Scale MC Samples with factor 2 if doing combination
    if (!name.Contains("unfolded") && doCombination) h_nom->Scale(2.);
    Sample.h_nominal = h_nom;
    Sample.h_nominal->Print();

    // get shape sys Templates
    TString systName_1d, systName_1u;
    for (auto& shapetuple : Sample.shapeSysts) {
      TString sys = get<0>(shapetuple) ;
      cout << sys << endl;
      if (contains(includeSystematics, sys)) {
        systName_1d = sys + "Up";
        systName_1u = sys + "Down";
        TH1* h_sysd = nullptr;
        if (!Sample.signal && !Sample.name.Contains("ATLAS")) h_sysd = (TH1*) file->Get(histname + "_" + systName_1d);
        else if (!Sample.signal && Sample.name.Contains("ATLAS")) h_sysd = (TH1D*) ATLASfile->Get(histname + "_" + systName_1d);
        else if (Sample.signal) h_sysd = (TH1*) Signalfile->Get(histname + "_" + systName_1d);
        //check, whether variation exists
        if (h_sysd == nullptr) {
          cout << "couldn't open " << histname + "_"  + systName_1d << endl;
          continue;
        }

        TH1* h_sysu = nullptr;
        if (!Sample.signal && !Sample.name.Contains("ATLAS")) h_sysu = (TH1*) file->Get(histname + "_" + systName_1u);
        else if (!Sample.signal && Sample.name.Contains("ATLAS")) h_sysu = (TH1D*) ATLASfile->Get(histname + "_" + systName_1u);
        else if (Sample.signal) h_sysu = (TH1*) Signalfile->Get(histname + "_" + systName_1u);
        //check, whether variation exists
        if (h_sysu == nullptr) {
          cout << "couldn't open " << histname + "_"  + systName_1u << endl;
          continue;
        }

        //Scale MC Samples with factor 2 if doing combination
        if (!name.Contains("unfolded") && doCombination) h_sysd->Scale(2.);

        TH1* diff_1d = (TH1*) h_sysd->Clone();
        diff_1d->Add(h_nom, -1);

        //Scale MC Samples with factor 2 if doing combination
        if (!name.Contains("unfolded") && doCombination) h_sysu->Scale(2.);

        TH1* diff_1u = (TH1*) h_sysu->Clone();
        diff_1u->Add(h_nom, -1);

        for (Int_t i = NbinMin; i <= Nbins; i++) {
          // double sigma = 0.5 * (fabs(diff_1d->GetBinContent(i)) + fabs(diff_1u->GetBinContent(i)));
          double sigma = std::max({fabs(diff_1d->GetBinContent(i)), fabs(diff_1u->GetBinContent(i))});
          // cout << "down: " << diff_1d->GetBinContent(i) << " up: " << diff_1u->GetBinContent(i) << " = " << sigma << endl;
          diff_1d->SetBinContent(i, sigma);
        }
        diff_1d->SetName("sigma_" + histname + "_" + sys);
        get<1>(shapetuple) = diff_1d;
      }
    }
  }

  //calculate total expected according to nuisance parameters
  //for example: n_obs' = n_obs*alpha_XS*alpha_Lumi*.. + sigma_syst1*nui_1 + sigma_syst2*nui2 + ..
  for (auto Sample : Samples) {
    Sample.h_nominal->Print();
    cout << "integral: " << Sample.h_nominal->Integral() << endl;
    TH1* nom = Sample.h_nominal;
    for (Int_t iBin = NbinMin ; iBin <= Nbins; iBin++) {
      double n_nom = nom->GetBinContent(iBin);
      // if (Sample.signal) n_nom /= 2.;
      cout << n_nom << " for Sample " << Sample.name << endl;
      double n_nom_sigma = nom->GetBinError(iBin);
      cout << n_nom_sigma << endl;
      if (doCombination and Sample.name == "unfolded") {
        n_nom += Samples[2].h_nominal->GetBinContent(iBin);
      }

      if (!Sample.signal && !Sample.name.Contains("unfolded")) { // add stat error of bkgs to cov matrix
        double staterror = n_nom_sigma * n_nom_sigma;
        cout << "adding " << staterror << " to diagonal of cov for " << Sample.name << " in Bin " << iBin << endl;
        cov(iBin - NbinMin, iBin - NbinMin) += staterror;
        cout << "diagonal at " << iBin - NbinMin << " is: " << cov(iBin - NbinMin, iBin - NbinMin) << endl;
      }
      RooRealVar* realVar_n_nom = (RooRealVar*)w.factory(TString::Format("n_nom_%s_Bin%d[%.3e]", (Sample.name).Data(), iBin, n_nom));
      realVar_n_nom->setVal(n_nom);
      realVar_n_nom->setConstant(true);
      // cout << realVar_n_nom->getVal() << " vs " << n_nom << endl;
      if (Sample.name == "unfolded") observables.add(*realVar_n_nom); //TODO dont hardcode
      TString n_nom_LNNui_string =  TString::Format("prod:n_nom_LNNui_%s_Bin%d(n_nom_%s_Bin%d,", (Sample.name).Data(), iBin, (Sample.name).Data(), iBin);
      // LN nuicances
      for (auto& lnsys : Sample.LNSysts) { // nominal x alpha
        n_nom_LNNui_string += TString::Format("alpha_%s,", ("NUI_" + get<0>(lnsys)).Data());
      }
      // multiply with r if signal
      if (Sample.signal == true) {
        n_nom_LNNui_string += "r,";
      }
      //remove trailing comma
      n_nom_LNNui_string = n_nom_LNNui_string.Strip(TString::kTrailing, ',');
      n_nom_LNNui_string += ")";
      // cout << n_nom_LNNui_string << endl;
      w.factory(n_nom_LNNui_string);

      // shape nuicances
      // if (!Sample.name.Contains("ATLAS")) {
      TString nexp_string =  TString::Format("sum:n_exp_%s_Bin%d(n_nom_LNNui_%s_Bin%d,", (Sample.name).Data(), iBin, (Sample.name).Data(), iBin);
      for (auto& shapesys : Sample.shapeSysts) {
        TString sysname = get<0>(shapesys);
        TH1* h_sys = get<1>(shapesys);
        // h_sys->Print();
        double n_sigma = h_sys->GetBinContent(iBin);
        RooRealVar* realVar_n_sigma = (RooRealVar*)w.factory(TString::Format("n_sigma_%s_%s_Bin%d[%.3e]",
                                      (Sample.name).Data(), sysname.Data(),  iBin, n_sigma));
        realVar_n_sigma->setVal(n_sigma);
        realVar_n_sigma->setConstant(true);
        // cout << n_sigma << " vs. " << realVar_n_sigma->getVal() << endl;
        // sigma * nuicance
        TString n_sigmuNui_str = TString::Format("prod:n_sigmaNui_%s_%s_Bin%d(n_sigma_%s_%s_Bin%d, NUI_%s",
                                 (Sample.name).Data(), sysname.Data(), iBin, (Sample.name).Data(), sysname.Data(),  iBin, sysname.Data());
        if (Sample.signal == true) {
          n_sigmuNui_str += ",r)";
        }
        else n_sigmuNui_str += ")";
        w.factory(n_sigmuNui_str);
        nexp_string += TString::Format("n_sigmaNui_%s_%s_Bin%d,", (Sample.name).Data(), sysname.Data(), iBin);
      }
      nexp_string = nexp_string.Strip(TString::kTrailing, ',');
      nexp_string += ")";
      // cout << nexp_string << endl;
      RooRealVar* realVar_n_exp = (RooRealVar*)w.factory(nexp_string);
      // }
    }
  }




  //bkg+signal+unfShapeSysts //TODO LNSYS for unfolded?!?!?!?
  for (Int_t iBin = NbinMin ; iBin <= Nbins; iBin++) {
    TString str_bkg_signal_unfShapeSys_sum = TString::Format("sum:bkg_signal_unfShape_sum_Bin%d(", iBin);
    for (auto Sample : Samples) { //bkg + sig
      if (!Sample.name.Contains("unfolded")) str_bkg_signal_unfShapeSys_sum += TString::Format("n_exp_%s_Bin%d,", (Sample.name).Data(), iBin );
    }
    for (auto& shapesys : unfolded.shapeSysts) { //unfShapeSysts
      TString sysname = get<0>(shapesys);
      str_bkg_signal_unfShapeSys_sum += TString::Format("n_sigmaNui_unfolded_%s_Bin%d,", sysname.Data(), iBin);
    }
    if (doCombination) {
      for (auto& shapesys : unfoldedATLAS.shapeSysts) { //unfShapeSysts
        TString sysname = get<0>(shapesys);
        str_bkg_signal_unfShapeSys_sum += TString::Format("n_sigmaNui_unfoldedATLAS_%s_Bin%d,", sysname.Data(), iBin);
      }
    }

    //remove trailing comma
    str_bkg_signal_unfShapeSys_sum = str_bkg_signal_unfShapeSys_sum.Strip(TString::kTrailing, ',');
    str_bkg_signal_unfShapeSys_sum += ")";
    // cout << str_bkg_signal_unfShapeSys_sum << endl;
    RooRealVar* str_bkg_signal_unfShapeSys_sum_Bin_realvar = (RooRealVar*)w.factory(str_bkg_signal_unfShapeSys_sum);
    muVec.add(*str_bkg_signal_unfShapeSys_sum_Bin_realvar);
  }
  // muVec.Print("v");


  cout << "### COV ###" << endl;
  cov.Print();

  TMatrixDSym covi(cov);
  covi.Invert();
  cout << "### COV inv ###" << endl;
  covi.Print();


  TDecompSVD svd(covi);
  svd.Decompose();
  TMatrixD u = svd.GetU();
  TMatrixD ut(TMatrixD::kTransposed, u);
  cout << "### UT ###" << endl;
  ut.Print();
  // TMatrixD v = svd.GetV();
  // TMatrixD vt(TMatrixD::kTransposed,v);
  TVectorD vec_s = svd.GetSig();

  TMatrixD s(NBinsCov, NBinsCov);
  TMatrixDDiag d(s);
  d = vec_s;
  cout << "### S ###" << endl;
  s.Print();
  vec_s.Print();

  if (doEigenStuff) {
    // do U^T * obs and U^T * exp to get them in eigenbasis
    int iBin0 = 0; // need second index starting at zero for matrix acces
    for ( Int_t iBin = NbinMin; iBin <= Nbins; iBin++ ) {
      double n_obs_eigen = 0.;
      double n_bs_eigen = 0.;

      TString str_bkg_signal_unfShapeSys_sum_eigen = TString::Format("expr:bkg_signal_unfShapeSys_sum_eigen_Bin%d('", iBin);

      int jBin0 =  0; // need second index starting at zero for matrix acces
      for ( Int_t jBin = NbinMin; jBin <= Nbins; jBin++ ) {
        double n_obs_j = w.var(TString::Format("n_nom_unfolded_Bin%d", jBin))->getVal();
        double n_bs_j = w.function(TString::Format("bkg_signal_unfShape_sum_Bin%d", jBin))->getVal();
        n_obs_eigen += n_obs_j * ut(iBin0, jBin0);
        n_bs_eigen += n_bs_j * ut(iBin0, jBin0) ;

        str_bkg_signal_unfShapeSys_sum_eigen += TString::Format("bkg_signal_unfShape_sum_Bin%d*%.3e+", jBin, ut(iBin0, jBin0));
        cout << n_obs_j << " * " << ut(iBin0, jBin0) << endl;
        jBin0++;
      }
      cout << "--> " << n_obs_eigen << endl;
      cout << "--> " << n_bs_eigen << endl;

      //Remove trailing plus
      str_bkg_signal_unfShapeSys_sum_eigen = str_bkg_signal_unfShapeSys_sum_eigen.Strip(TString::kTrailing, '+');
      str_bkg_signal_unfShapeSys_sum_eigen += "',";

      for ( Int_t jBin = NbinMin; jBin <= Nbins; jBin++ ) {
        str_bkg_signal_unfShapeSys_sum_eigen += TString::Format("bkg_signal_unfShape_sum_Bin%d,", jBin);
      }
      //Remove trailing comma
      str_bkg_signal_unfShapeSys_sum_eigen = str_bkg_signal_unfShapeSys_sum_eigen.Strip(TString::kTrailing, ',');
      str_bkg_signal_unfShapeSys_sum_eigen += ")";
      cout << "str_bkg_signal_unfShapeSys_sum_eigen Bin " << iBin << ":" << endl;
      cout << str_bkg_signal_unfShapeSys_sum_eigen << endl;
      w.factory(str_bkg_signal_unfShapeSys_sum_eigen);

      double n_min, n_max;
      if (n_bs_eigen > 0.) {
        n_min = std::min(0., -10 * n_bs_eigen);
        n_max = 10 * n_bs_eigen;
      }
      else {
        n_max = std::max(0., - 10 * n_bs_eigen);
        n_min = 10 * n_bs_eigen;
      }

      // RooRealVar* realVar_n_obs_eigen = (RooRealVar*)w.factory(TString::Format("n_obs_eigen_%d[%.3e,%.3e,%.3e]", iBin, n_obs_eigen, n_min, n_max));
      RooRealVar* realVar_n_obs_eigen = (RooRealVar*)w.factory(TString::Format("n_obs_eigen_%d[%.3e]", iBin, n_obs_eigen));
      realVar_n_obs_eigen->setVal(n_obs_eigen);
      realVar_n_obs_eigen->setConstant(true);
      observables_eigen.add(*realVar_n_obs_eigen);
      // cout << realVar_n_obs_eigen->getVal() << " vs " << n_obs_eigen << endl;

      w.factory(TString::Format("Gaussian:pdf_Bin_%d(n_obs_eigen_%d,bkg_signal_unfShapeSys_sum_eigen_Bin%d,%.3e)", iBin, iBin, iBin, 1. / sqrt(s(iBin0, iBin0))));
      w.pdf(TString::Format("pdf_Bin_%d", iBin))->Print("v");
      iBin0++;
    }

    //multiply each BinPDF with each Nuicance PDF
    TString string_model = "PROD:modelCMS(";
    for ( Int_t iBin = NbinMin; iBin <= Nbins; iBin++ ) {
      string_model += TString::Format("pdf_Bin_%d,", iBin);
    }
    string_model += "nuisancePdf";
    string_model = string_model.Strip(TString::kTrailing, ',');
    string_model += ")";
    // cout << string_model.Data() << endl;

    w.factory(string_model);
  }

  RooMultiVarGaussian mvg("mvg", "mvg", observables, muVec, cov);
  mvg.Print();
  w.import(mvg);
  w.factory("PROD:model2(mvg,nuisancePdf)");

  if (doEigenStuff) {
    // standard gaussians in eigenbasis

    RooDataSet data("data", "", observables_eigen);
    data.add(observables_eigen);
    w.import(data);
    w.import(observables_eigen);

    ModelConfig model_b("ModelConfig_bonly", &w);
    model_b.SetPdf                 ( *w.pdf("modelCMS") );
    model_b.SetParametersOfInterest( RooArgSet( *w.var("r") ) );
    model_b.SetObservables         ( observables_eigen );
    model_b.SetGlobalObservables   ( globalObservables );
    model_b.SetNuisanceParameters  ( nuisances );

    poi->setVal(0.);
    model_b.SetSnapshot(*poi);

    w.import(model_b);

    ModelConfig model_bs("ModelConfig", &w);

    model_bs.SetPdf                 ( *w.pdf("modelCMS") );
    model_bs.SetParametersOfInterest( RooArgSet( *w.var("r") ) );
    model_bs.SetObservables         ( observables_eigen );
    model_bs.SetGlobalObservables   ( globalObservables );
    model_bs.SetNuisanceParameters  ( nuisances );

    poi->setVal(1.);
    model_bs.SetSnapshot(*poi);
    w.import(model_bs);

    const auto* list_nuis = model_bs.GetNuisanceParameters();
    if (list_nuis) {
      cout << "Found " << list_nuis->getSize() << " nuisance parameters:" << endl;
      list_nuis->writeToStream(std::cout, false);

    }

    const auto* list_global = model_bs.GetGlobalObservables();
    if (list_global) {
      cout << "Found " << list_global->getSize() << " global observables:" << endl;
      list_global->writeToStream(std::cout, false);

    }

    const auto* list_obs = model_bs.GetObservables();
    if (list_obs) {
      cout << "Found " << list_obs->getSize() << " observables:" << endl;
      list_obs->writeToStream(std::cout, false);

    }
  }
  else {
    //MultiVariateGaussian
    RooDataSet data("data", "", observables);
    data.add(observables);
    w.import(data);

    ModelConfig model_b("ModelConfig_bonly", &w);
    model_b.SetPdf                 ( *w.pdf("model2") );
    model_b.SetParametersOfInterest( RooArgSet( *w.var("r") ) );
    model_b.SetObservables         ( observables );
    model_b.SetGlobalObservables   ( globalObservables );
    model_b.SetNuisanceParameters  ( nuisances );

    poi->setVal(0.);
    model_b.SetSnapshot(*poi);

    w.import(model_b);

    ModelConfig model_bs("ModelConfig", &w);
    model_bs.SetPdf                 ( *w.pdf("model2") );
    model_bs.SetParametersOfInterest( RooArgSet( *w.var("r") ) );
    model_bs.SetObservables         ( observables );
    model_bs.SetGlobalObservables   ( globalObservables );
    model_bs.SetNuisanceParameters  ( nuisances );

    poi->setVal(1.);
    model_bs.SetSnapshot(*poi);
    w.import(model_bs);

    const auto* list_nuis = model_bs.GetNuisanceParameters();
    if (list_nuis) {
      cout << "Found " << list_nuis->getSize() << " nuisance parameters:" << endl;
      list_nuis->writeToStream(std::cout, false);

    }

    const auto* list_global = model_bs.GetGlobalObservables();
    if (list_global) {
      cout << "Found " << list_global->getSize() << " global observables:" << endl;
      list_global->writeToStream(std::cout, false);

    }

    const auto* list_obs = model_bs.GetObservables();
    if (list_obs) {
      cout << "Found " << list_obs->getSize() << " observables:" << endl;
      list_obs->writeToStream(std::cout, false);

    }
  }
  // w.Print();
  // printf("Done printing the current workspace\n");




  // write workspace in the file (recreate file if already existing)
  w.writeToFile(fileName, true);
  // w.Print();




}



