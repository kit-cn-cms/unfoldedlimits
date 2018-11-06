using namespace std;

void tweakCMS() {

	TFile*  file  = new TFile("/nfs/dust/cms/user/swieland/Darkmatter/DM_Unfolding/rootfiles/data_normedmuRmuF.root", "OPEN");
	// TFile*  file  = new TFile("data_CMSsmeared.root","OPEN");

	TRandom3 rand;
	// TFile *file = new TFile("hsimple.root");
	TFile* file_tweaked = new TFile("ATLAS_CMSCopy.root", "RECREATE");
	double bins_met[13] = {200., 240., 300., 380., 470., 570., 670., 790., 910., 1040., 1180., 1330., 2000. };
	TObject *obj;
	TKey *key;
	TIter next(file->GetListOfKeys());
	while ((key = (TKey *) next())) {
		obj = file->Get(key->GetName()); // copy object to memory
		// do something with obj
		printf(" found object:%s\n", key->GetName());
		// obj->Print();
		TH1 *h = (TH1*)key->ReadObj();
		TH1* hnew = (TH1*)h->Clone();
		// TH1* hnew = new TH1D("hist", "hist", 12, bins_met );
		int iBin0 = 1;
		for ( Int_t iBin = 1; iBin <= h->GetNbinsX(); iBin++) {
			// cout << h->GetBinContent(iBin) << endl;
			double content = h->GetBinContent(iBin);
			if ( content != 0 && iBin0 < 12) {
				// double g = rand.Poisson(content);
				// cout << "old: " << content << " new: " << g << endl;
				hnew->SetBinContent(iBin0, content);
				iBin0++;
			}
		}
		TString histoname = TString(key->GetName());
		histoname.ReplaceAll("unfolded", "unfoldednotATLAS");
		histoname.ReplaceAll("Gen_Hadr_Recoil_Pt_", "Gen_Hadr_Recoil_Pt_ATLAS_");
		hnew->SetName(histoname);
		hnew->Write();
		cout << hnew->GetNbinsX() << endl;
	}

	// file->ls();
	file->Close();


	file  = new TFile("/nfs/dust/cms/user/swieland/Darkmatter/DM_Unfolding/rootfiles/MCdata.root", "OPEN");

	file_tweaked = new TFile("ATLAS_CMSCopyMCData.root", "RECREATE");
	TIter next2(file->GetListOfKeys());
	while ((key = (TKey *) next2())) {
		obj = file->Get(key->GetName()); // copy object to memory
		// do something with obj
		printf(" found object:%s\n", key->GetName());
		// obj->Print();
		TH1 *h = (TH1*)key->ReadObj();
		TH1* hnew = (TH1*)h->Clone();
		// TH1* hnew = new TH1D("hist", "hist", 12, bins_met );
		int iBin0 = 1;
		for ( Int_t iBin = 1; iBin <= h->GetNbinsX(); iBin++) {
			// cout << h->GetBinContent(iBin) << endl;
			double content = h->GetBinContent(iBin);
			if ( content != 0 && iBin0 < 12) {
				// double g = rand.Poisson(content);
				// cout << "old: " << content << " new: " << g << endl;
				hnew->SetBinContent(iBin0, content);
				iBin0++;
			}
		}
		TString histoname = TString(key->GetName());
		histoname.ReplaceAll("unfolded", "unfoldednotATLAS");
		histoname.ReplaceAll("Gen_Hadr_Recoil_Pt_", "Gen_Hadr_Recoil_Pt_ATLAS_");
		hnew->SetName(histoname);
		hnew->Write();
		cout << hnew->GetNbinsX() << endl;
	}

	// file->ls();
	file->Close();


}
