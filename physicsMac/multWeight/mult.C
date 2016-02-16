void mult(){

	//TFile *inFile = TFile::Open("MultCorrAll.root");
	//TFile *inFile = TFile::Open("MultCorrNoList2.root");
	TFile *inFile = TFile::Open("../../../MultCorrWeight2.root");	

	TH1F *mult = (TH1F*)inFile->Get("mult_runid")->Clone("mult");
	TH1F *event = (TH1F*)inFile->Get("event_runid")->Clone("event");
	TH1F *mult_cent = (TH1F*)inFile->Get("mult_runid_cent")->Clone("mult_cent");
	TH1F *event_cent = (TH1F*)inFile->Get("event_runid_cent")->Clone("event_cent");
	TH1F *mult_cut = (TH1F*)inFile->Get("mult_runid_cent_cut")->Clone("mult_cut");
  TH1F *event_cut = (TH1F*)inFile->Get("event_runid_cent_cut")->Clone("event_cut");
/*
	TH1F *mult_zoom = (TH1F*)inFile->Get("mult_runid_zoom")->Clone("mult_zoom");
	TH1F *mult_zoom_corr = (TH1F*)inFile->Get("mult_runid_zoom_corr")->Clone("mult_zoom_corr");
	TH1F *event_zoom = (TH1F*)inFile->Get("event_runid_zoom")->Clone("event_zoom");
	TH1F *mult_centcut_zoom = (TH1F*)inFile->Get("mult_runid_cent_cut_zoom")->Clone("mult_centcut_zoom");
	TH1F *mult_centcut_zoom_corr = (TH1F*)inFile->Get("mult_runid_cent_cut_zoom_corr")->Clone("mult_centcut_zoom_corr");
	TH1F *event_centcut_zoom = (TH1F*)inFile->Get("event_runid_cent_cut_zoom")->Clone("event_centcut_zoom");
*/
	TH1F *multCorr = (TH1F*)inFile->Get("multCorr")->Clone("multCorr");
	TH1F *multCorr_cent = (TH1F*)inFile->Get("multCorr_cent")->Clone("multCorr_cent");
	TH1F *multCorr_cent_cut = (TH1F*)inFile->Get("multCorr_cent_cut")->Clone("multCorr_cent_cut");

cout<<"entries = "<<mult_cut->GetEntries()<<endl;

	/*TH1F *weight = (TH1F*)inWeight->Get("hWeight")->Clone("weight");
	TH1F *weight2 = (TH1F*)inWeight2->Get("hWeight")->Clone("weight2");
	TH1F *weight3 = (TH1F*)inWeight3->Get("hWeight")->Clone("weight3");
*/
	mult->Divide(event);
	mult_cent->Divide(event_cent);
	mult_cut->Divide(event_cut);
/*
	mult_zoom->Divide(event_zoom);
	mult_zoom_corr->Divide(event_zoom);
	mult_centcut_zoom->Divide(event_centcut_zoom);
	mult_centcut_zoom_corr->Divide(event_centcut_zoom);
*/
cout<<"entries = "<<mult_cut->GetEntries()<<endl;

	TCanvas *can = new TCanvas("can","can",800,500);
	mult->SetTitle(0);
	mult_cut->SetTitle(0);
	mult_cent->SetTitle(0);
	mult->GetYaxis()->SetTitle("refMult/Events");
	mult->GetXaxis()->SetTitle("runID");
	//mult->GetXaxis()->SetRangeUser(13128950,13129100);
	//mult_cent->GetXaxis()->SetRangeUser(13128950,13129100);
	//mult_cut->GetXaxis()->SetRangeUser(13128950,13129100);
	//mult->GetXaxis()->SetRangeUser(13128000,13128100);
	//mult_cent->GetXaxis()->SetRangeUser(13128000,13128100);
	//mult_cut->GetXaxis()->SetRangeUser(13128000,13128100);
	mult->GetXaxis()->SetNdivisions(505);
	mult_cut->GetXaxis()->SetNdivisions(505);
	mult_cent->GetXaxis()->SetNdivisions(505);
	mult_cut->SetLineColor(kRed);
	mult_cent->SetLineColor(kBlue);
	mult_cut->SetStats(0);
	mult_cut->Draw();
	mult_cent->Draw("same");
	mult->Draw("same");

	TLegend *leg = new TLegend(0.6, 0.4, 0.8, 0.55);
	leg->AddEntry(mult,"MB");
	leg->AddEntry(mult_cent,"central5");
	leg->AddEntry(mult_cut,"central5+StRefMultCorr");
	leg->Draw("same");

	//can->SaveAs("multEvent_vs_RunID_allMuDst.pdf");
/*
	mult_zoom->GetXaxis()->SetTitle("refMult/Events");
	mult_zoom->GetYaxis()->SetTitle("runID");
	mult_zoom->GetXaxis()->SetNdivisions(505);
	mult_centcut_zoom->GetXaxis()->SetNdivisions(505);
	mult_zoom_corr->GetXaxis()->SetNdivisions(505);
	mult_centcut_zoom_corr->GetXaxis()->SetNdivisions(505);
	mult_zoom->SetLineColor(kBlack);
	mult_centcut_zoom->SetLineColor(kRed);
	mult_zoom_corr->SetLineColor(kBlue);
	mult_centcut_zoom_corr->SetLineColor(kMagenta);
	//mult_centcut_zoom->Draw();
	//mult_zoom->Draw("same");
	//mult_zoom_corr->Draw("same");
	//mult_centcut_zoom_corr->Draw("same");
	
	//can->SaveAs("multEvent_vs_RunID_allMuDst_zoom_noList2.pdf");
*/
	multCorr->Scale(1./multCorr->Integral(600,800));
	multCorr_cent->Scale(1./multCorr_cent->Integral(600,800));
	multCorr_cent_cut->Scale(1./multCorr_cent_cut->Integral(600,800));

	multCorr->SetStats(0);
	multCorr->SetTitle(0);
	multCorr->GetXaxis()->SetTitle("refMult");
	multCorr->SetLineColor(kBlack);
	multCorr_cent->SetLineColor(kBlue+1);
	multCorr_cent_cut->SetLineColor(kRed+2);
	multCorr_cent_cut->SetFillStyle(3005);
	multCorr_cent_cut->SetFillColor(kRed+2);
	multCorr->SetMinimum(1e-5);
	multCorr->Draw();	
	multCorr_cent->Draw("same");
	multCorr_cent_cut->Draw("same");
	multCorr_cent_cut->Draw("hist same");

	TLegend *leg2 = new TLegend(0.13, 0.2, 0.4, 0.4);
	leg2->AddEntry(multCorr, "MB");
	leg2->AddEntry(multCorr_cent, "central5");
	leg2->AddEntry(multCorr_cent_cut, "central5+StRegMultCorr");
	leg2->Draw("same");

	can->SetLogy();
	//can->SaveAs("refMult_allMuDst_NoList3.pdf");

	//..calculate weight

	TH1F *weight = new TH1F("weight","weight; MB/central5", 800, 0, 800);
	
	int nbin = multCorr->GetNbinsX();
	cout<<"nbin = "<<nbin<<endl;
	for(int i=0; i<nbin; i++)
	{

		double w = multCorr->GetBinContent(i+1)/multCorr_cent_cut->GetBinContent(i+1);
		if (multCorr_cent_cut->GetBinContent(i+1) == 0 || multCorr->GetBinContent(i+1) == 0) w = 0;
		weight->SetBinContent(i+1, w);

		double a = multCorr->GetBinError(i+1)/multCorr->GetBinContent(i+1);
		double b = multCorr_cent_cut->GetBinError(i+1)/multCorr_cent_cut->GetBinContent(i+1);
		double error = w*TMath::Sqrt(a*a + b*b);
		if (multCorr_cent_cut->GetBinContent(i+1) == 0 || multCorr->GetBinContent(i+1) == 0) error = 0;
		//if(i > 500) cout<<"weight error = "<<error<<endl;
		weight->SetBinError(i+1, error);

	}

	TCanvas *can2 = new TCanvas("can2","",600,450);

	weight->GetXaxis()->SetRangeUser(520, 800);
	weight->SetMaximum(3);
	weight->Draw("P");

	TLine *line = new TLine(520, 1, 800, 1);
	line->SetLineStyle(2);
	line->Draw("same");

	//can2->SaveAs("Weight_allMuDstNoList3.pdf");

/*	TFile *f = new TFile("weight_allMuDst_noList2.root","recreate");
	weight->Write();
	f->Close();

	delete f;
*/
	/*weight->SetLineColor(kBlack);
	weight2->SetLineColor(kRed+2);
	weight3->SetLineColor(kBlue);

	TLegend *l = new TLegend(0.4, 0.7, 0.6, 0.85);
	l->AddEntry(weight,"dataset 2");
	l->AddEntry(weight2,"dataset 3");
	l->AddEntry(weight3,"dataset 1");

	weight->GetXaxis()->SetRangeUser(420,750);
	weight->SetStats(0);
	weight->Draw();
	weight2->Draw("same");
	weight3->Draw("same");
	l->Draw();
	
	can->SaveAs("weights.pdf");
*/

	
}
