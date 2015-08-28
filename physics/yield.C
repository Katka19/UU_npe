void customize_histo(TH1F *h, int marker, int color);
void customize_graph(TGraphErrors *h, int marker, int color);
//------------------------------------------------
void yield(){

gStyle->SetLegendBorderSize(0);
gStyle->SetFillColor(-1);

TFile *DanielFile = TFile::Open("../Daniel_AuAu/Spectra_pho_npe_ratio/Data/npe_spec.root");
TFile *vzFile = TFile::Open("physics_central5.root");
TFile *ptFile = TFile::Open("pT_spectra.root");
//TFile *reffFile = TFile::Open("../Reconstruction_efficiency/eff_pho_vs_pt.root"); //..using Au+Au embedding
//TFile *single_reffFile = TFile::Open("../AuAu_embed/reff_single.root"); //..using Au+Au embedding
TFile *efficiencyFile = TFile::Open("efficiencies_with_errors.root");
TFile *emcFile = TFile::Open("../UU_embedding/single_rec/EMC_eff_embedding.root");
TFile *singleFile = TFile::Open("../UU_embedding/single_rec/single_rec_efficiency.root");
TFile *phoFile = TFile::Open("../UU_embedding/photo_rec/photonic_rec_efficiency.root");

TH1F *Daniel_npe_phe = (TH1F*)DanielFile->Get("R_var0");
TH1F *inc_ele = (TH1F*)ptFile->Get("electron_pt")->Clone("inc_ele");
TH1F *pho_ele = (TH1F*)ptFile->Get("photonic_pt")->Clone("pho_ele");
//TF1 *eff_pho = (TF1*)reffFile->Get("fit_0_5")->Clone("eff_pho");
//TF1 *eff_single = (TF1*)single_reffFile->Get("fun4")->Clone("eff_single");
//TH1F *reff_single = (TH1F*)single_reffFile->Get("efficiency")->Clone("reff_single");
TH1F *purity = (TH1F*)efficiencyFile->Get("PurityEff")->Clone("purity");
TH1F *vz = (TH1F*)vzFile->Get("Vz")->Clone("vz");
//TH1F *emc_eff = (TH1F*)efficiencyFile->Get("EMC_eff")->Clone("emc_eff");
TH1F *emc_eff = (TH1F*)emcFile->Get("EMC_combined")->Clone("emc_eff");
TH1F *nSigma_eff = (TH1F*)efficiencyFile->Get("nsigma_eff")->Clone("nSigma_eff");
TH1F *single_eff = (TH1F*)singleFile->Get("finalStat")->Clone("single_eff");
//TH1F *single_eff = (TH1F*)efficiencyFile->Get("Single_eff")->Clone("single_eff");
//TH1F *photo_eff = (TH1F*)efficiencyFile->Get("Photo_eff")->Clone("photo_eff");
TH1F *photo_eff = (TH1F*)phoFile->Get("final")->Clone("photo_eff");
TH1F *JPsi_yield = (TH1F*)efficiencyFile->Get("final_jpsi")->Clone("JPsi_yield");

double Range_array[15] = {1.2, 1.4, 1.6, 1.8, 2, 2.4, 2.8, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6};
int n_bin = 14;
//double Shift_array[14] = {1.29342, 1.49393, 1.69437, 1.89476, 2.18104, 2.5831, 2.98477, 3.38613, 3.78727, 4.18824, 4.58907, 4.98979, 5.39043, 5.79099};
double Shift_array[14] = {1.3, 1.5, 1.7, 1.9, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0, 5.4, 5.8};
double Correction[14] = {1.01455, 1.0162, 1.01554, 1.01664, 1.01542, 1.02414, 1.00108, 1.03622, 1.00381, 1.04126, 1.004, 0.901204, 0.973212, 0.915594};

TH1F *inc_reb = inc_ele->Rebin(n_bin, "inc_reb", Range_array);
TH1F *pho_reb = pho_ele->Rebin(n_bin, "pho_reb", Range_array);

double inc[15];
double inc_error[15];
double pho[15];
double pho_error[15];
double pur[15];
double pur_error[15];
double effPho[15];
double effPho_error[15];
double npe[15];
double npe_error[15];

double yield[15];
double yield_subtr[15];
double yield_jpsi[15];
double jpsi_contribution[15];
double yield_jpsi_error[15];
double yield_error_stat[15];
double error_syst[15];
double yield_error_syst[15];
double pT[15];
double delta_pT[15];
double effNsigma[15];
double nsigma_err[15];
double effEmc[15];
double emc_err[15];
double effEmcHighpt = 0.555;
double effSingle[15];
double single_err[15];
double deltay = 1.4;

double num[15];
double denom[15];
double ratio[15];
double error1[15];
double ratio_error[15];

double rel_stat_err[15];
double rel_syst_err[15];
double syst1[15];
double syst2[15];
double syst3[15];
double syst4[15];
double syst5[15];
double syst6[15];
double syst_comb[15];

int bin1 = vz->GetXaxis()->FindBin(-30);
int bin2 = vz->GetXaxis()->FindBin(30);
double N_events = vz->Integral(bin1,bin2);

TH1F *n_npe = new TH1F("n_npe","n_npe;#font[12]{p_{T}} (GeV/#font[12]{c}); N_{npe}",n_bin, Range_array);
TH1F *yield_npe = new TH1F("yield_npe","yield_npe;#font[12]{p_{T}} (GeV/#font[12]{c}); 1/2#pip_{T}*d^{2}N/dp_{T}d#eta",n_bin, Range_array); //..for Raa
TH1F *yield_npe_syst = new TH1F("yield_npe_syst","yield_npe_syst;#font[12]{p_{T}} (GeV/#font[12]{c}); 1/2#pip_{T}*d^{2}n/dp_{T}d#eta",n_bin,Range_array); //..for Raa
TH1F *npe_over_phe = new TH1F("npe_over_phe","npe_over_phe;#font[12]{p_{T}} (GeV/#font[12]{c}); NPE/PHE",n_bin,Range_array);
TH1F *pho_phoEff = new TH1F("pho_phoEff","pho_phoEff;#font[12]{p_{T}} (GeV/#font[12]{c}); Counts",n_bin,Range_array);
TH1F *relative_error = new TH1F("relative_error","relative_error;#font[12]{p_{T}} (GeV/#font[12]{c});relative error [%]",n_bin,Range_array); //..statistical
TH1F *relative_error1 = new TH1F("relative_error1","relative_error1;#font[12]{p_{T}} (GeV/#font[12]{c});relative error [%]",n_bin,Range_array); //..EMC eff
TH1F *relative_error2 = new TH1F("relative_error2","relative_error2;#font[12]{p_{T}} (GeV/#font[12]{c});relative error [%]",n_bin,Range_array); //..nsigma eff
TH1F *relative_error3 = new TH1F("relative_error3","relative_error3;#font[12]{p_{T}} (GeV/#font[12]{c});relative error [%]",n_bin,Range_array); //..rec eff
TH1F *relative_error4 = new TH1F("relative_error4","relative_error4;#font[12]{p_{T}} (GeV/#font[12]{c});relative error [%]",n_bin,Range_array); //..purity
TH1F *relative_error5 = new TH1F("relative_error5","relative_error5;#font[12]{p_{T}} (GeV/#font[12]{c});relative error [%]",n_bin,Range_array); //..photonic
TH1F *relative_error_syst = new TH1F("relative_error_syst","relative_error_syst;#font[12]{p_{T}} (GeV/#font[12]{c});relative error [%]",n_bin,Range_array);

double bin_center[15] = {0};
double error[15] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};

for(int i=0; i<n_bin; i++){
	//..calculating N_npe
	inc[i] = inc_reb->GetBinContent(i+1);
	inc_error[i] = inc_reb->GetBinError(i+1);
	pho[i] = pho_reb->GetBinContent(i+1);
	pho_error[i] = pho_reb->GetBinError(i+1);
	pur[i] = purity->GetBinContent(i+1);
	pur_error[i] = purity->GetBinError(i+1);
	//effPho[i] = eff_pho->Eval((Range_array[i]+Range_array[i+1])/2); //..from Au+Au embedding
	effPho[i] = photo_eff->GetBinContent(i+1);
	effPho_error[i] = photo_eff->GetBinError(i+1);
	npe[i] = inc[i]*pur[i] - pho[i]/effPho[i];
	npe_error[i] = TMath::Sqrt((inc_error[i]*pur[i])*(inc_error[i]*pur[i])+(pho_error[i]/effPho[i])*(pho_error[i]/effPho[i])); //..statistical error of NPE

	n_npe->SetBinContent(i+1, npe[i]);
	n_npe->SetBinError(i+1, npe_error[i]);

	//..calculating NPE/PHE ratio
	double pho_over_eff = (pho[i]/effPho[i]);
	ratio[i] = npe[i]/(pho_over_eff);
	ratio_error[i] = ratio[i]*TMath::Sqrt((npe_error[i]/npe[i])*(npe_error[i]/npe[i]) + (pho_error[i]/pho[i])*(pho_error[i]/pho[i]));
	npe_over_phe->SetBinContent(i+1,ratio[i]);
	npe_over_phe->SetBinError(i+1,ratio_error[i]);

	pho_phoEff->SetBinContent(i+1, pho_over_eff);

	//..calculting yield of NPE
	pT[i] = (Range_array[i]+Range_array[i+1])/2;
	delta_pT[i] = (Range_array[i+1] - Range_array[i]);
	effEmc[i] = emc_eff->GetBinContent(i+1);
	emc_err[i] = emc_eff->GetBinError(i+1);
	effNsigma[i] = nSigma_eff->GetBinContent(i+1);
	nsigma_err[i] = nSigma_eff->GetBinError(i+1);
	effSingle[i] = single_eff->GetBinContent(i+1);
        single_err[i] = single_eff->GetBinError(i+1);

	yield[i] = (((1/TMath::TwoPi())*(1/N_events)*(1/pT[i])*(1/delta_pT[i])*(1/deltay)*(npe[i]/(effNsigma[i]*effEmc[i]*effSingle[i])))/2)*(Correction[i]);
	//yield[i] = ((1/TMath::TwoPi())*(1/N_events)*(1/pT[i])*(1/delta_pT[i])*(1/deltay)*(npe[i]/(effNsigma[i]*effEmc[i]*effSingle[i])))/2;
	//...subtract J/Psi contribution
	yield_subtr[i] = yield[i] - (JPsi_yield->GetBinContent(i+1));

	yield_jpsi[i] = JPsi_yield->GetBinContent(i+1);
	yield_jpsi_error[i] = JPsi_yield->GetBinError(i+1);

	jpsi_contribution[i] = (yield_jpsi[i]/yield[i])*100;
	//if(i>=10){yield[i] = ((1/TMath::TwoPi())*(1/N_events)*(1/pT[i])*(1/delta_pT[i])*(1/deltay)*(npe[i]/(effNsigma[i]*effEmcHighpt*effSingle[i])))/2;}

	yield_error_stat[i] = (((1/TMath::TwoPi())*(1/N_events)*(1/pT[i])*(1/delta_pT[i])*(1/deltay)*(npe_error[i]/(effNsigma[i]*effEmc[i]*effSingle[i])))/2)*(Correction[i]);
	//yield_error_stat[i] = ((1/TMath::TwoPi())*(1/N_events)*(1/pT[i])*(1/delta_pT[i])*(1/deltay)*(npe_error[i]/(effNsigma[i]*effEmc[i]*effSingle[i])))/2;	

	//error_syst[i] = (1/(effNsigma[i]*effEmc[i]*effSingle[i]))*TMath::Sqrt((emc_err[i]/effEmc[i])*(emc_err[i]/effEmc[i]) + (nsigma_err[i]/effNsigma[i])*(nsigma_err[i]/effNsigma[i]) + (single_err[i]/effSingle[i])*(single_err[i]/effSingle[i]));
	//yield_error_syst[i] = (yield[i])*error_syst[i];

	syst1[i] = (((1/TMath::TwoPi())*(1/N_events)*(1/pT[i])*(1/delta_pT[i])*(1/deltay)*(npe[i]/(effNsigma[i]*effEmc[i]*effEmc[i]*effSingle[i])))/2)*emc_err[i]; //..EMC
	syst2[i] = (((1/TMath::TwoPi())*(1/N_events)*(1/pT[i])*(1/delta_pT[i])*(1/deltay)*(npe[i]/(effNsigma[i]*effNsigma[i]*effEmc[i]*effSingle[i])))/2)*nsigma_err[i]; //..nsigma
	syst3[i] = (((1/TMath::TwoPi())*(1/N_events)*(1/pT[i])*(1/delta_pT[i])*(1/deltay)*(npe[i]/(effNsigma[i]*effEmc[i]*effSingle[i]*effSingle[i])))/2)*single_err[i]; //..single rec
	syst4[i] = ((1/TMath::TwoPi())*(1/N_events)*(1/pT[i])*(1/delta_pT[i])*(1/deltay)*((inc[i]*pur_error[i])/(effNsigma[i]*effEmc[i]*effSingle[i])))/2; //..purity
	syst5[i] = ((1/TMath::TwoPi())*(1/N_events)*(1/pT[i])*(1/delta_pT[i])*(1/deltay)*(((pho[i]/(effPho[i]*effPho[i]))*effPho_error[i])/(effNsigma[i]*effEmc[i]*effSingle[i])))/2; //..photonic rec eff
	syst6[i] = JPsi_yield->GetBinError(i+1); //..J/psi
	syst_comb[i] = TMath::Sqrt(syst1[i]*syst1[i] + syst2[i]*syst2[i] + syst3[i]*syst3[i] + syst4[i]*syst4[i] + syst5[i]*syst5[i] + syst6[i]*syst6[i])*Correction[i];
	//syst_comb[i] = TMath::Sqrt(syst1[i]*syst1[i] + syst2[i]*syst2[i] + syst3[i]*syst3[i] + syst4[i]*syst4[i] + syst5[i]*syst5[i] + syst6[i]*syst6[i]);

	yield_npe->SetBinContent(i+1, yield_subtr[i]);
	yield_npe->SetBinError(i+1, yield_error_stat[i]);
	yield_npe_syst->SetBinContent(i+1, yield_subtr[i]);
	//yield_npe_syst->SetBinError(i+1, yield_error_syst[i]);
	yield_npe_syst->SetBinError(i+1, syst_comb[i]);
	
	bin_center[i] = (Range_array[i]+Range_array[i+1])/2;

	rel_stat_err[i] = (yield_error_stat[i]/yield[i])*100; //..relative statistical error
	relative_error->SetBinContent(i+1, rel_stat_err[i]);
	//relative_error_syst->SetBinContent(i+1, (yield_error_syst[i]/yield[i])*100);
	rel_syst_err[i] = (syst_comb[i]/yield[i])*100; //..relative systematical error
	relative_error_syst->SetBinContent(i+1, rel_syst_err[i]);

	relative_error1->SetBinContent(i+1, (syst1[i]/yield[i])*100); //..EMC
	relative_error2->SetBinContent(i+1, (syst2[i]/yield[i])*100); //..nsigma
	relative_error3->SetBinContent(i+1, (syst3[i]/yield[i])*100); //..single rec
	relative_error4->SetBinContent(i+1, (syst4[i]/yield[i])*100); //..purity
	relative_error5->SetBinContent(i+1, (syst5[i]/yield[i])*100); //..photonic reconstruction eff

}

//..plot fits of efficiencies with confidence intervals
//TCanvas *syst_can = new TCanvas("syst_can","syst_can",600,450);
//syst_can->Divide(2,2);
TGraphErrors *graph_emc = new TGraphErrors(n_bin, bin_center, effEmc, 0, emc_err);
TGraphErrors *graph_nsigma = new TGraphErrors(n_bin, bin_center, effNsigma, 0, nsigma_err);
TGraphErrors *graph_single = new TGraphErrors(n_bin, bin_center, effSingle, 0, single_err);
TGraphErrors *graph_purity = new TGraphErrors(n_bin, bin_center, pur, 0, pur_error);
TGraphErrors *graph_photo = new TGraphErrors(n_bin, bin_center, effPho, 0, effPho_error);

TCanvas *emc_canvas = new TCanvas("emc_canvas","emc_canvas", 600,450);
customize_graph(graph_emc, 29, kOrange+7);
graph_emc->GetYaxis()->SetTitle("#font[12]{#epsilon_{emc}}");
graph_emc->SetFillColor(kOrange-4);
graph_emc->SetMaximum(1);
graph_emc->SetMinimum(0);
//syst_can->cd(1);
graph_emc->Draw("A3");
graph_emc->Draw("LXsame");
emc_canvas->SaveAs("final_plots/emc_systematics.pdf");

TCanvas *nsigma_canvas = new TCanvas("nsigma_canvas","nsigma_canvas",600,450);
customize_graph(graph_nsigma, 29, kMagenta+2);
graph_nsigma->GetYaxis()->SetTitle("#font[12]{#epsilon_{n#sigma_e}}");
graph_nsigma->SetFillColor(kMagenta-9);
graph_nsigma->SetMaximum(1);
graph_nsigma->SetMinimum(0);
//syst_can->cd(2);
graph_nsigma->Draw("A3");
graph_nsigma->Draw("LXsame");
nsigma_canvas->SaveAs("final_plots/nsigma_systematics.pdf");

TCanvas *single_canvas = new TCanvas("single_canvas","single_canvas",600,450);
customize_graph(graph_single, 29, kGreen+2);
graph_single->GetYaxis()->SetTitle("#font[12]{#epsilon_{rec}}");
graph_single->SetFillColor(kGreen-6);
graph_single->SetMaximum(1);
graph_single->SetMinimum(0);
//syst_can->cd(3);
graph_single->Draw("A3");
graph_single->Draw("LXsame");
single_canvas->SaveAs("final_plots/single_systematics.pdf");

TCanvas *purity_canvas = new TCanvas("purity_canvas","purity_canvas",600,450);
customize_graph(graph_purity, 29, kRed);
graph_purity->GetYaxis()->SetTitle("purity");
graph_purity->SetFillColor(kRed-9);
graph_purity->SetMaximum(1);
graph_purity->SetMinimum(0);
//syst_can->cd(4);
graph_purity->Draw("A3");
graph_purity->Draw("LXsame");
purity_canvas->SaveAs("final_plots/purity_systematics.pdf");

//syst_can->SaveAs("final_plots/systematics_plots_UUemb.pdf");

TCanvas *photo_canvas = new TCanvas("photo_canvas","photo_canvas", 600, 450);
customize_graph(graph_photo, 29, kCyan+3);
graph_photo->GetYaxis()->SetTitle("#font[12]{#epsilon_{pho}}");
graph_photo->SetFillColor(kCyan-3);
graph_photo->SetMaximum(1);
graph_photo->SetMinimum(0);
graph_photo->Draw("A3");
graph_photo->Draw("LXsame");
photo_canvas->SaveAs("final_plots/photonic_systematics.pdf");
//.....................

//..save relative statistical errors of UU yield for raa.C
customize_histo(relative_error, 29, 1);
relative_error->SetMarkerSize(1.5);
relative_error->SetMaximum(100);

TFile *relative = new TFile("UU_relative_stat_error.root","recreate");
relative_error->Write();
relative->Close();
delete relative;
//....................

//..plot relative systematical errors from all contributions
TCanvas *rel_syst_can = new TCanvas("rel_syst_can", "rel_syst_can", 600, 450);
customize_histo(relative_error1, 29, kOrange+7);
customize_histo(relative_error2, 29, kMagenta+2);
customize_histo(relative_error3, 29, kGreen+2);
customize_histo(relative_error4, 29, kRed);
customize_histo(relative_error5, 29, kCyan+2);
customize_histo(relative_error_syst, 29, kBlue+3);

relative_error1->SetMaximum(70);
relative_error1->SetMinimum(-1);

relative_error1->Draw("P");
relative_error2->Draw("Psame");
relative_error3->Draw("Psame");
relative_error4->Draw("Psame");
relative_error5->Draw("Psame");
relative_error_syst->Draw("Psame");

TLegend *lege_rel_syst = new TLegend(0.15, 0.6, 0.7, 0.85);
lege_rel_syst->AddEntry(relative_error1, "systematic from EMC efficiency","p");
lege_rel_syst->AddEntry(relative_error2, "systematic from n#sigma_{e} efficiency","p");
lege_rel_syst->AddEntry(relative_error3, "systematic from single track rec. eff.","p");
lege_rel_syst->AddEntry(relative_error4, "systematic from purity","p");
lege_rel_syst->AddEntry(relative_error5, "systematic from photonic rec. eff.","p");
lege_rel_syst->AddEntry(relative_error_syst, "systematic from all efficiences","p");
lege_rel_syst->SetTextSize(0.04);
lege_rel_syst->Draw("same");

//rel_syst_can->SaveAs("final_plots/relative_systematic_errors_comb.pdf");
rel_syst_can->SaveAs("final_plots/relative_systematic_errors_comb_UUemb_Jpsi.pdf");
//.....................

//..plot ratio of NPE/PHE
TCanvas *npe_over_phe_can = new TCanvas("npe_over_phe_can","npe_over_phe_can",600,450);
customize_histo(npe_over_phe, 29, 1);
npe_over_phe->SetMinimum(0);
npe_over_phe->SetMaximum(4);
npe_over_phe->Draw("");

customize_histo(Daniel_npe_phe, 29, 2);
Daniel_npe_phe->SetMinimum(0);
Daniel_npe_phe->Draw("EPsame");

TLegend *lege = new TLegend(0.5, 0.7, 0.85, 0.85);
lege->AddEntry(npe_over_phe, "U+U centrality 0-5%");
lege->AddEntry(Daniel_npe_phe, "Au+Au centrality 0-5%");
lege->SetTextSize(0.045);
lege->Draw("same");

npe_over_phe_can->SaveAs("final_plots/ratio_npe_phe_UUemb.pdf");
//....................

//..plot yield of NPE after J/psi subtraction, NPE and J/psi yield together, JPsi contribution in %
TGraphErrors *graph_yield_npe = new TGraphErrors(n_bin, Shift_array, yield, 0, yield_error_stat);
//TGraphErrors *graph_yield_npe_syst = new TGraphErrors(n_bin, bin_center, yield, error, yield_error_syst);
TGraphErrors *graph_yield_npe_syst_subtr = new TGraphErrors(n_bin, Shift_array, yield_subtr, error, syst_comb);
TGraphErrors *graph_yield_npe_syst = new TGraphErrors(n_bin, Shift_array, yield, error, syst_comb);
TGraphErrors *graph_yield_npe_subtr = new TGraphErrors(n_bin, Shift_array, yield_subtr, 0, yield_error_stat);
TGraphErrors *graph_yield_jpsi = new TGraphErrors(n_bin, bin_center, yield_jpsi, 0, yield_jpsi_error);
TGraphErrors *graph_jpsi_contribution = new TGraphErrors(n_bin, bin_center, jpsi_contribution, 0, 0);

TCanvas *yield_can = new TCanvas("yield_can","yield_can", 600,450);
yield_can->SetLogy();
customize_graph(graph_yield_npe_subtr, 29, kBlue+3);
customize_graph(graph_yield_npe_syst_subtr, 29, kBlue+3);
graph_yield_npe_syst_subtr->SetFillColor(kBlue-9);
graph_yield_npe_syst_subtr->Draw("AP2");
graph_yield_npe_subtr->SetFillColor(0);
graph_yield_npe_subtr->SetMarkerSize(1.2);
graph_yield_npe_subtr->Draw("PZsame");

TLegend *lege_yield = new TLegend(0.3, 0.8, 0.9, 0.8);
lege_yield->AddEntry((TObject*)0, "U+U centrality 0-5% (ToF+ZDC)","");
lege_yield->SetTextSize(0.045);
lege_yield->SetFillColor(0);
lege_yield->Draw("same");
cout<<"before drawing UU yield"<<endl;
yield_can->SaveAs("final_plots/yield_05_syst_comb_UUemb_JPsi_binshift.pdf");
yield_can->SaveAs("final_plots/yield_05_syst_comb_UUemb_JPsi_binshift.png");

customize_graph(graph_yield_npe, 29, kBlue+3);
customize_graph(graph_yield_npe_syst, 29, kBlue+3);
customize_graph(graph_yield_jpsi, 29, kRed+2);
graph_yield_npe_syst->SetFillColor(kBlue-9);
graph_yield_npe_syst->Draw("AP2");
graph_yield_npe->Draw("PZsame");
graph_yield_jpsi->SetFillColor(0);
graph_yield_jpsi->Draw("PZsame");
lege_yield->Draw("same");

TLegend *lege_yield2 = new TLegend(0.5, 0.55, 0.85, 0.7);
lege_yield2->AddEntry(graph_yield_npe_syst,"yield of NPE","LPF");
lege_yield2->AddEntry(graph_yield_jpsi, "yield of J/#psi","LPF");
lege_yield2->SetTextSize(0.035);
lege_yield2->SetFillColor(0);
lege_yield2->Draw("same");

yield_can->SaveAs("final_plots/yield_npe_plus_JPsi.pdf");

yield_can->SetLogy(0);
customize_graph(graph_jpsi_contribution, 29, kRed+2);
graph_jpsi_contribution->GetYaxis()->SetTitle("J/#psi contribution [%]");
graph_jpsi_contribution->SetMaximum(100);
graph_jpsi_contribution->SetMinimum(0);
graph_jpsi_contribution->Draw("APZ");
//lege_yield->AddEntry((TObject*)0, "J/#psi contribution to NPE yield");
lege_yield->Draw("same");

yield_can->SaveAs("final_plots/JPsi_contribution_to_npe.pdf");
//....................

//..save histograms of yield of NPE for raa.C
TFile *yield_file = new TFile("final_plots/yield_file.root","recreate");
yield_npe->Write();
yield_npe_syst->Write();
yield_file->Close();
delete yield_file;
//....................


//..plot NPE and PHE
TCanvas *npe_pho_can = new TCanvas ("npe_pho_can","npe_pho_can",600,450);
customize_histo(pho_phoEff, 29, 4);
customize_histo(n_npe, 29, 1);
npe_pho_can->SetLogy();

TLegend *legend = new TLegend(0.3, 0.65, 0.85, 0.85);
legend->AddEntry(pho_phoEff, "photonic electrons");
legend->AddEntry(n_npe, "non-photonic electrons");
legend->AddEntry((TObject*)0,"U+U centrality 0-5%","");
legend->SetTextSize(0.045);

pho_phoEff->Draw("EP");
n_npe->Draw("EPsame");
legend->Draw("same");

npe_pho_can->SaveAs("final_plots/npe_pho_UUemb.pdf");
//..................

//..plot inc and pho electrons
double norm = 1;
inc_reb->Scale(norm,"width");
pho_reb->Scale(norm,"width");

TCanvas *cinc = new TCanvas("cinc","cinc",600,450);
inc_reb->SetMarkerStyle(29);
inc_reb->SetMarkerSize(1.2);
inc_reb->SetLineColor(kBlue);
inc_reb->SetMarkerColor(kBlue);
inc_reb->SetMinimum(1);
inc_reb->GetXaxis()->SetTitle("p_{T} (GeV/#font[12]{c})");
inc_reb->GetYaxis()->SetTitle("#bf{#bf{#frac{Counts}{bin width} (#frac{1}{GeV})}}");
pho_reb->SetMarkerStyle(29);
pho_reb->SetMarkerSize(1.2);
pho_reb->SetLineColor(kRed);
pho_reb->SetMarkerColor(kRed);
inc_reb->Draw("EP");
pho_reb->Draw("EPsame");

TLegend *incleg = new TLegend(0.55,0.75, 0.85, 0.85);
incleg->AddEntry(inc_reb,"inclusive e^{+} + e^{-}");
incleg->AddEntry(pho_reb,"photonic e^{+} + e^{-}");
incleg->SetFillColor(0);
incleg->Draw("same");

cinc->SetLogy();
cinc->SaveAs("final_plots/inc_pho_spectra.pdf");
cinc->SaveAs("final_plots/inc_pho_spectra.png");
//..............................
}
//--------------------------------------
//--------------------------------------
//---------------------------------------
void customize_histo(TH1F *h, int marker, int color){

h->GetXaxis()->SetTitleSize(0.045);
h->GetXaxis()->CenterTitle(true);
h->GetXaxis()->SetLabelSize(0.04);
h->GetYaxis()->SetTitleSize(0.05);
h->GetYaxis()->CenterTitle(true);
h->GetYaxis()->SetLabelSize(0.04);
h->SetStats(0);
h->SetTitle(0);
//h->SetMinimum(0);
//h->SetMaximum(2);
h->SetMarkerStyle(marker);
h->SetMarkerColor(color);
h->SetLineColor(color);

}
//---------------------------------------
void customize_graph(TGraphErrors *h, int marker, int color){

h->GetXaxis()->SetTitleSize(0.045);
h->GetXaxis()->SetTitle("p_{T} (GeV/#font[12]{c})");
h->GetXaxis()->CenterTitle(true);
h->GetXaxis()->SetLabelSize(0.04);
h->GetYaxis()->SetTitleSize(0.045);
h->GetYaxis()->SetTitle("1/2#pip_{T}*d^{2}N/dp_{T}d#eta (GeV/#font[12]{c}^{2})");
h->GetYaxis()->CenterTitle(true);
h->GetYaxis()->SetLabelSize(0.04);
h->SetTitle(0);
//h->SetMinimum(0);
//h->SetMaximum(2);
h->SetMarkerStyle(marker);
h->SetMarkerColor(color);
h->SetLineColor(color);

}


