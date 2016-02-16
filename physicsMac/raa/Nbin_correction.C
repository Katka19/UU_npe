void Nbin_correction(){

 //TFile *inFile = TFile::Open("out_all1000Again.root");
 //TFile *inFile = TFile::Open("out_allMoreStat.root");
 //TFile *inFile = TFile::Open("out_combined.root");
 TFile *inFile = TFile::Open("out_allHists.root");
 //TFile *inFile = TFile::Open("out_MoreStat_1000Again.root");
 
 TH1F *multCentral = (TH1F*)inFile->Get("mult_cent")->Clone("multCentral");
 TH1F *multMB = (TH1F*)inFile->Get("mult_mb")->Clone("multMB");
 //TH1F *multAll = (TH1F*)inFile->Get("mult")->Clone("multAll");
 TH1F *multCut = (TH1F*)inFile->Get("mult_cent_cut")->Clone("multCut");
 TH1F *centrality = (TH1F*)inFile->Get("centrality")->Clone("centrality");

 cout<<"number of entries mult_cent = "<<multCentral->GetEntries()<<endl;
 cout<<"number of entries = "<<centrality->GetEntries()<<endl;
 cout<<"underflow bin centrality = "<<centrality->GetBinContent(0)<<endl;
 cout<<"70-80% = "<<centrality->GetBinContent(1)<<endl;
 cout<<"60-70% = "<<centrality->GetBinContent(2)<<endl;
 cout<<"50-60% = "<<centrality->GetBinContent(3)<<endl;
 cout<<"40-50% = "<<centrality->GetBinContent(4)<<endl;
 cout<<"30-40% = "<<centrality->GetBinContent(5)<<endl;
 cout<<"20-30% = "<<centrality->GetBinContent(6)<<endl;
 cout<<"10-20% = "<<centrality->GetBinContent(7)<<endl;
 cout<<"5-10% = "<<centrality->GetBinContent(8)<<endl;
 cout<<"0-5% = "<<centrality->GetBinContent(9)<<endl;

 TCanvas *c = new TCanvas("c","c",600,450);
/*
 centrality->SetStats(0);
 centrality->SetTitle(0);
 centrality->GetXaxis()->SetTitle("centrality bin");
 centrality->GetYaxis()->SetTitle("Counts");
 centrality->Draw("hist");
 
 gStyle->SetLegendBorderSize(0);
 TLegend *legCent = new TLegend(0.8, 0.15, 0.9, 0.2);
 legCent->AddEntry((TObject*)0, "0-5%","");
 legCent->SetFillColor(0);
 legCent->Draw("same");

 TLegend *legCent3 = new TLegend(0.72, 0.15, 0.8, 0.2);
 legCent3->AddEntry((TObject*)0, "5-10%","");
 legCent3->SetFillColor(0);
 legCent3->Draw("same");

 TLegend *legCent2 = new TLegend(0.62, 0.15, 0.72, 0.2);
 legCent2->AddEntry((TObject*)0, "10-20%","");
 legCent2->SetFillColor(0);
 legCent2->Draw("same");

 c->SetLogy();
 c->SaveAs("centrality.pdf");*/

 multCentral->SetLineColor(kBlue);
 multCentral->SetStats(0);
 multCentral->SetTitle(0);
 multCentral->GetXaxis()->SetTitle("refMult");
 multCentral->SetMaximum(5000);
 //multAll->SetLineColor(kBlack);
 //multAll->SetStats(0);
 //multAll->SetTitle(0);
 multMB->SetLineColor(kBlack);
 multMB->GetXaxis()->SetTitle("refMult");
 multMB->SetStats(0);
 multMB->SetTitle(0);
 multCut->SetLineColor(kRed+1);
 multCut->SetStats(0);
 multCut->SetTitle(0);
 multCut->SetFillStyle(3005);
 //multCut->SetFillColor(kGreen+3);

 //double bin1 = multMB->GetXaxis()->FindBin(610);
 //double bin2 = multMB->GetXaxis()->FindBin(700);
 double mbInt = multMB->Integral(600,700);
 cout<<"integral v MB = "<<mbInt<<endl;
 double centInt = multCentral->Integral(600,700);
 cout<<"integral central 5 = "<<centInt<<endl;

 multMB->Scale(1/mbInt);
 multCentral->Scale(1/centInt);
 multCut->Scale(1/centInt);

 multCentral->SetMaximum(1);
 multCentral->GetXaxis()->SetRangeUser(500,700);
 multMB->GetXaxis()->SetRangeUser(500,700);
 multCut->GetXaxis()->SetRangeUser(500,700);
 //multAll->Draw("hist");
 multCentral->Draw("EP");
 multMB->Draw("EPsame");
 multCut->Draw("EPsame");

 gStyle->SetLegendBorderSize(-1);
 TLegend *leg = new TLegend(0.5, 0.7, 0.85, 0.85);
 leg->AddEntry(multCut, "StRefMultCorr cut, central 0-5%","lf");
 leg->AddEntry(multCentral, "central 0-5%","l");
 leg->AddEntry(multMB, "min bias","l");
 leg->SetFillColor(0);
 leg->Draw("same");

 c->SetLogy();
 c->SaveAs("../multWeight/refMult_scaled_allHists.pdf");

 //..calculate the correction factor
 /*c->cd();
 multMB->GetXaxis()->SetRangeUser(500,700);
 multCentral->GetXaxis()->SetRangeUser(500,700);
 multCut->GetXaxis()->SetRangeUser(500,700);

 multMB->Draw("");
 multCut->Draw("same");

 c->SaveAs("refMult_zoom.pdf");*/

 //TH1F *hWeight = new TH1F("hWeight","hWeight",200,500,700);
 TH1F *hWeight = new TH1F("hWeight","hWeight",800,0,800);
 hWeight->Sumw2();

 //double weight[200] = {0}; 

 double binStart = mult_mb->GetXaxis()->FindBin(500);
 cout<<"binStart = "<<binStart<<endl;
 double binEnd = mult_mb->GetXaxis()->FindBin(700);
 cout<<"binEnd = "<<binEnd<<endl;
 
 for(int iw = 0; iw<800; iw++){

	double w_mb = multMB->GetBinContent(iw+1);
	double mb_err = multMB->GetBinError(iw+1);
	if(iw > 460){
		double w_cent = multCentral->GetBinContent(iw+1);
		double cent_err = multCentral->GetBinError(iw+1);
	}
	else{
		double w_cent = 0;
		double cent_err = 0;
	}

	//weight[iw-binStart] = (w_mb/w_cent);
	if(w_cent == 0 || w_mb == 0){
		double weight = 0;
		double err = 0;
	}
	else{
		//double weight = w_mb/w_cent;
		double weight = w_mb/w_cent;
		double err = weight*TMath::Sqrt((mb_err/w_mb)*(mb_err/w_mb) + (cent_err/w_cent)*(cent_err/w_cent));
	}
	hWeight->SetBinContent(iw+1, weight);
	hWeight->SetBinError(iw+1, err);

 }

 TFile *result = new TFile("weight_factor_allHists.root","recreate");
 hWeight->Write();
 result->Close();
 delete result;

 c->cd();
 c->SetLogy(0);
 hWeight->GetXaxis()->SetRangeUser(500,700);
 hWeight->Draw("EP");
 //c->SaveAs("weight_factor_1000.pdf");
 
}
