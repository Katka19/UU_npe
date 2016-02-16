void customize_histo(TH1F *h, int marker, int color);
void customize_graph(TGraphErrors *g, int marker, int color);
//--------------------------
void raa(){

gStyle->SetLegendBorderSize(0);
//gStyle->SetFillColor(-1);

TFile *yieldFile = TFile::Open("yield_file_rebinned_pp12.root");
TFile *DanielFile = TFile::Open("../../Daniel_AuAu/Spectra_pho_npe_ratio/Data/nph_yield_Raa_central_trigger_only_0_5pct.root");
//TFile *DanielFile = TFile::Open("../../Daniel_AuAu/Spectra_pho_npe_ratio/Data/nph_yield_Raa.root");
TFile *ppCorrFile = TFile::Open("../raa/final_corr_pp_rebin_pp12.root");
TFile *efficiencyFile = TFile::Open("../efficiencies_with_errors.root");
TFile *Dfile = TFile::Open("../nph_yield_pp_subtract_jpsi_DY.root");

TH1F *uuyield = (TH1F*)yieldFile->Get("yield_npe")->Clone("uuyield");
TH1F *uuyield_syst = (TH1F*)yieldFile->Get("yield_npe_syst")->Clone("uuyield_syst");
TH1F *DanielRaa = (TH1F*)DanielFile->Get("Raa");//..0-5% centrality
TH1F *DanielRaaSys = (TH1F*)DanielFile->Get("Raa_sys");//..0-5% centrality
TH1F *DanielRaaSyspp = (TH1F*)DanielFile->Get("Raa_sys_pp");
TH1F *ppCor = (TH1F*)ppCorrFile->Get("final_correction")->Clone("ppCor");
TH1F *auyield = (TH1F*)DanielFile->Get("nph_rebin")->Clone("auyield");
TH1F *auyield_syst = (TH1F*)DanielFile->Get("nph_sys_rebin")->Clone("auyield_syst");
//TH1F *auyield = (TH1F*)DanielFile->Get("nph_rebin0")->Clone("auyield");
//TH1F *auyield_syst = (TH1F*)DanielFile->Get("nph_sys_rebin0")->Clone("auyield_syst");
//TH1F *pp09 = (TH1F*)efficiencyFile->Get("pp_spectra09")->Clone("pp09");
TH1F *pp09 = (TH1F*)Dfile->Get("nph_pp")->Clone("ppD");

TGraphErrors *ppyield_run12_syst = new TGraphErrors("Run12_cross_section_Raa.dat","%lg %lg %lg");//je to naopak ako tvrdi Xiaozhi
TCanvas *c12 = new TCanvas("c12","c12",600,450);
ppyield_run12_syst->SetMarkerStyle(21);
//ppyield_run12->Draw("AP");
c12->SetLogy();
pp09->Scale(42);
pp09->GetXaxis()->SetRangeUser(3.0,6.0);
pp09->Draw("EP");
ppyield_run12_syst->Draw("Psame");
c12->SaveAs("pp12_vs_pp0508.pdf");

TGraphErrors *ppyield_run12 = new TGraphErrors("Run12_cross_section_Raa_syst.dat","%lg %lg %lg");//je to naopak ako tvrdi Xiaozhi

Double_t Rebin_array[17] = {1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0};
double Center_array[16] = {1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75}; 
int n_rebin = 16;

//double N_bin = 1341; //..from the paper
//double N_bin = 1281;//..from Zhenyu (Hiroshi)
double N_bin = 1156.166;
double N_bin_error = 105;//..from the paper
double N_bin_systError = 125.094;//..difference between correction and Hiroshi's value
double N_binE = TMath::Sqrt((N_bin_error)*(N_bin_error) + (N_bin_systError)*(N_bin_systError));
cout<<"N_binE = "<<N_binE<<endl;

double raa[17] = {0};
double raa_error_uustat[17] = {0};
double raa_error_uusyst[17] = {0};
double raa_error_ppsyst[17] = {0};

double raa_error_UUpp[17] = {0};

double pp_fit[17] = {0};
double pp_corr[17] = {0};
double pp_error[17] = {0};

double uu[17] = {0};
double uu_error_stat[17] = {0};
double uu_error_syst[17] = {0};

double x_error[17] = {0};
double x_error_pp[17] = {0};

double bin_center[17] = {0};
double bin_center_auau[17] = {0};

double rel_err_uu[17] = {0};
double rel_err_pp[17] = {0};
double rel_err_uu_syst[17] = {0};

double pp[17] = {0};
double pp_err[17] = {0};
double pp_err_stat[17] = {0};
double pp_err_syst[17] = {0};

double auau[17] = {0};
double auau_stat[17] = {0};
double auau_syst[17] = {0};
double daniel_raa[17] = {0};
double daniel_raa_stat[17] = {0};
double daniel_raa_syst[17] = {0};
double daniel_raa_ppsyst[17] = {0};
double daniel_error_UUpp[17] = {0};

double pp_stat_corr[17] = {0};
double pp_syst_corr[17] = {0};

TH1F *pp12 = new TH1F("nph_pp","nph_pp",n_rebin, Rebin_array);

TFile *DanielppFile = TFile::Open("../nph_yield_pp_subtract_jpsi_DY.root");
TH1F *Dpp = (TH1F*)DanielppFile->Get("nph_pp")->Clone("Dpp");
cout<<"Dpp[10] = "<<Dpp->GetBinContent(10)<<endl;
double ppp[17] = {0};
double Errppp[17] = {0};

//..Raa from U+U collisions
for(int i=0;i<n_rebin;i++){

	double xxx = 0;
	ppyield_run12->GetPoint(i,xxx,pp[i]);
cout<<"bin pp = "<<xxx<<endl;
	pp[i] = pp[i]/42;
	pp_corr[i] = ppCor->GetBinContent(i+1);
	pp_fit[i] = pp[i]*pp_corr[i];
	pp_err_stat[i] = (ppyield_run12->GetErrorY(i))/42;
	pp_stat_corr[i] = pp_err_stat[i]*pp_corr[i];
	pp_err_syst[i] = (ppyield_run12_syst->GetErrorY(i))/42;
 	pp_syst_corr[i] = pp_err_syst[i]*pp_corr[i];
	pp_err[i] = TMath::Sqrt(pp_err_stat[i]*pp_err_stat[i] + pp_err_syst[i]*pp_err_syst[i]);
	pp_error[i] = pp_err[i]*pp_corr[i];

	pp12->SetBinContent(i+1,pp[i]);
	pp12->SetBinError(i+1,pp_err[i]);

	uu[i] = uuyield->GetBinContent(i+1);
	uu_error_stat[i] = uuyield->GetBinError(i+1);
	uu_error_syst[i] = uuyield_syst->GetBinError(i+1);

//	if(i<9){auau[i] = -5;
//		auau_stat[i] = 0;
//		auau_syst[i] = 0;
//	}else{
		double bin = auyield->GetXaxis()->FindBin(Rebin_array[i]+0.01);
cout<<"binau  = "<<Rebin_array[i]+0.01<<endl;
		double bin_sys = auyield_syst->GetXaxis()->FindBin(Rebin_array[i]+0.01);
		auau[i] = auyield->GetBinContent(bin);
		auau_stat[i] = auyield->GetBinError(bin);
		auau_syst[i] = auyield_syst->GetBinError(bin_sys);
//	}
	raa[i] = (1/N_bin)*(uu[i]/(pp_fit[i]));
	double uustat = (1/N_bin)*(uu_error_stat[i]/pp_fit[i]);
	double ppstat = raa[i]*(pp_stat_corr[i]/pp_fit[i]);
	raa_error_uustat[i] = TMath::Sqrt(uustat*uustat + ppstat*ppstat);
	//raa_error_uustat[i] = (1/N_bin)*(uu_error_stat[i]/pp_fit[i]);
	double uusyst = (1/N_bin)*(uu_error_syst[i]/pp_fit[i]);
	double ppsyst = raa[i]*(pp_syst_corr[i]/pp_fit[i]);
	raa_error_uusyst[i] = TMath::Sqrt(uusyst*uusyst + ppsyst*ppsyst);
	//raa_error_uusyst[i] = (1/N_bin)*(uu_error_syst[i]/pp_fit[i]);
	//raa_error_ppsyst[i] = (1/N_bin)*(uu[i]/((pp_fit[i])*(pp_fit[i])))*(pp_error[i]);
	raa_error_ppsyst[i] = raa[i]*(pp_error[i]/pp_fit[i]);
//cout<<"pp error = "<<raa_error_ppsyst[i]<<endl;

//	raa_error_UUpp[i] = TMath::Sqrt((raa_error_uusyst[i])*(raa_error_uusyst[i]) + (raa_error_ppsyst[i])*(raa_error_ppsyst[i]));
//cout<<"total syst error = "<<raa_error_UUpp[i]<<endl;

//	if(i<9){daniel_raa[i] = -5;
//		daniel_raa_stat[i] = -5;
//		daniel_raa_syst[i] = -5;
//		daniel_error_UUpp[i] = -5;
//	}else{
//		double bin1pp = Dpp->GetXaxis()->FindBin(Rebin_array[i]+0.1);
//		ppp[i] = Dpp->GetBinContent(bin1pp);
//		Errppp[i] = Dpp->GetBinError(bin1pp);
		daniel_raa[i] = (auau[i]/pp[i]);
//		daniel_raa[i] = (auau[i]/ppp[i]);
		double austat = auau_stat[i]/pp[i];
		//double austat = auau_stat[i]/pp[i];
		double ppaustat = daniel_raa[i]*(pp_err_stat[i]/pp[i]);
		daniel_raa_stat[i] = TMath::Sqrt(austat*austat + ppaustat*ppaustat);
		//daniel_raa_stat[i] = (auau_stat[i]/pp[i]);
		double ausyst = auau_syst[i]/pp[i];
		double auppsyst = daniel_raa[i]*(pp_err_syst[i]/pp[i]);
		//daniel_raa_syst[i] = (auau_syst[i]/pp[i]);
		daniel_raa_syst[i] = TMath::Sqrt(ausyst*ausyst + auppsyst*auppsyst);
		daniel_raa_ppsyst[i] = daniel_raa[i]*(pp_err[i]/pp[i]);
//		daniel_error_UUpp[i] = TMath::Sqrt((daniel_raa_syst[i])*(daniel_raa_syst[i]) + (daniel_raa_ppsyst[i])*(daniel_raa_ppsyst[i]));
//	}
cout<<"UU/au = "<<ppsyst/auppsyst<<endl;
	bin_center[i] = (Rebin_array[i]+Rebin_array[i+1])/2;
	if(i>7){bin_center_auau[i] = bin_center[i]+0.1;}
	else{bin_center_auau[i] = bin_center[i];}

	x_error[i] = 0.05;
	x_error_pp[i] = 0.07;

	double nbin_error = (raa[i]/N_bin)*N_binE;
	cout<<"error from nbin = "<<nbin_error/raa[i]<<endl;

}

TFile *file = new TFile("pp12_yield.root","recreate");
pp12->Write();
file->Close();
delete file;

//double bin_center_auau[7] = {2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75};

//..plot final Raa (U+U and Au+Au together)
TCanvas *graph_can = new TCanvas("graph_can","graph_can",600,450);
graph_can->SetLogy();
TGraphErrors *raa_uu_stat = new TGraphErrors(n_rebin, bin_center, raa, 0, raa_error_uustat);
TGraphErrors *raa_uu_syst = new TGraphErrors(n_rebin, bin_center, raa, x_error, raa_error_uusyst);// raa_error_UUpp);
TGraphErrors *raa_daniel = new TGraphErrors(n_rebin, bin_center_auau, daniel_raa, 0, daniel_raa_stat);
TGraphErrors *raa_daniel_sys = new TGraphErrors(n_rebin, bin_center_auau, daniel_raa, x_error, daniel_raa_syst);

TGraphErrors *cronin_upper = new TGraphErrors("../raa/eR-U_cron1_eloss0_S200D+Ball-Y0.0.75.evolve420New","%lg %lg");
TGraphErrors *cronin_lower = new TGraphErrors("../raa/eR-U_cron1_eloss0_S200D+Ball-Y0.1.evolve420New","%lg %lg");
TGraphErrors *eloss_upper = new TGraphErrors("../raa/eR-U_cron1.5_eloss1_S200D+Ball-Y0.0.75.evolve420New","%lg %lg");
TGraphErrors *eloss_lower = new TGraphErrors("../raa/eR-U_cron1.5_eloss1_S200D+Ball-Y0.1.evolve420New","%lg %lg");

cronin_upper->SetLineColor(kGreen+2);
cronin_upper->SetLineWidth(3);
cronin_upper->SetLineStyle(9);
cronin_lower->SetLineColor(kGreen+2);
cronin_lower->SetLineWidth(3);
cronin_lower->SetLineStyle(9);
eloss_upper->SetLineColor(kMagenta+2);
eloss_upper->SetLineWidth(3);
eloss_upper->SetLineStyle(5);
eloss_lower->SetLineColor(kMagenta+2);
eloss_lower->SetLineWidth(3);
eloss_lower->SetLineStyle(5);

customize_graph(raa_uu_stat, 29, kBlue+3);
customize_graph(raa_uu_syst, 29, kBlue+3);
raa_uu_syst->SetFillColor(kBlue-9);
customize_graph(raa_daniel, 29, kRed+3);
customize_graph(raa_daniel_sys, 29, kRed+3);
raa_daniel_sys->SetFillColor(kRed-9);

raa_uu_syst->SetMaximum(15);
raa_uu_syst->SetMinimum(1e-01);
raa_uu_syst->Draw("AP2");
raa_daniel_sys->Draw("P2same");
//raa_uu_syst->Draw("P2same");
raa_daniel->Draw("PZsame");
raa_uu_stat->Draw("PZsame");
cronin_upper->Draw("Lsame");
cronin_lower->Draw("Lsame");
eloss_upper->Draw("Lsame");
eloss_lower->Draw("Lsame");

TLine *line = new TLine(0.9, 1, 6.2, 1);
line->SetLineStyle(2);
line->SetLineColor(1);
line->SetLineWidth(1);

double x_axis[1] = {6.2};
double y_axis[1] = {1};
double x_axis_err[1] = {0.05};
double y_axis_err[1] = {0.14126};//{0.09};//..use relative error, which is the same for all points
TGraphErrors *Nbin_err = new TGraphErrors(1, x_axis, y_axis, x_axis_err, y_axis_err);
Nbin_err->SetFillColor(kGreen+3);
Nbin_err->Draw("2same");

double xau_axis[1] = {6.1};
double yau_axis_err[1] = {0.0262};
TGraphErrors *Nbin_errAu = new TGraphErrors(1, xau_axis, y_axis, x_axis_err, yau_axis_err);
Nbin_errAu->SetFillColor(kGreen-3);
Nbin_errAu->Draw("2same");

TLegend *legend = new TLegend(0.15,0.77,0.73,0.89);//0.77->0.83 when no AuAu data are plotted
legend->AddEntry(raa_uu_syst, "U+U, centrality 0-5% (ToF+ZDC)","LPF");
legend->AddEntry(raa_daniel_sys, "Au+Au, centrality 0-5%","PLF");
legend->SetTextSize(0.04);
legend->SetFillColor(0);

line->Draw("same");
legend->Draw("same");

TLegend *leg_nbin = new TLegend(0.4, 0.48, 0.85, 0.62);
leg_nbin->AddEntry(Nbin_err, "syst. uncertainty from N_{bin}^{UU}","F");
leg_nbin->AddEntry(Nbin_errAu, "syst. uncertainty from N_{bin}^{AuAu}","F");
leg_nbin->SetTextSize(0.04);
leg_nbin->SetFillColor(0);
leg_nbin->Draw("same");

TLegend *leg_syst = new TLegend(0.4, 0.65, 0.85, 0.75);
leg_syst->AddEntry(cronin_lower,"I.Vitev; large Cronin","L");//, CNM e-loss)+(QGP e-loss, dissoc.)","L");
leg_syst->AddEntry(eloss_lower,"I.Vitev; small Cronin","L");//, CNM e-loss)+(QGP e-loss, dissoc.)","L");
leg_syst->SetTextSize(0.04);
leg_syst->SetFillColor(0);
leg_syst->Draw("same");

TText *text = new TText();
text->SetNDC();
text->SetTextFont(1);
text->SetTextColor(kRed);
text->SetTextSize(0.04);
text->SetTextAlign(22);
text->SetTextAngle(0);
text->DrawText(0.25, 0.72, "STAR Preliminary");

graph_can->SaveAs("raa_central05_pp12_log_lowpTAu_preliminary.pdf");
graph_can->SaveAs("raa_central05_pp12_log_lowpTAu_preliminary.png");
//....................

}
//-------------------------------------------
//-------------------------------------------
//--------------------------------------------
void customize_histo(TH1F *h, int marker, int color){

h->GetXaxis()->SetTitleSize(0.045);
h->GetXaxis()->CenterTitle(true);
h->GetXaxis()->SetLabelSize(0.037);
h->GetYaxis()->SetTitleSize(0.045);
h->GetYaxis()->CenterTitle(true);
h->GetYaxis()->SetLabelSize(0.037);
h->SetStats(0);
h->SetTitle(0);
//h->SetMinimum(0);
//h->SetMaximum(2);
h->SetMarkerStyle(marker);
h->SetMarkerColor(color);
h->SetLineColor(color);

}
//---------------------------------------------
void customize_graph(TGraphErrors *g, int marker, int color){

g->GetXaxis()->SetTitleSize(0.045);
g->GetXaxis()->SetTitle("p_{T} (GeV/#font[12]{c})");
g->GetXaxis()->CenterTitle(true);
g->GetXaxis()->SetLabelSize(0.037);
g->GetYaxis()->SetTitle("R_{AA}");
g->GetYaxis()->SetTitleSize(0.045);
g->GetYaxis()->CenterTitle(true);
g->GetYaxis()->SetLabelSize(0.037);
g->SetFillColor(0);
g->SetTitle(0);
g->SetMinimum(0);
g->SetMaximum(2);
g->SetMarkerStyle(marker);
g->SetMarkerSize(1.5);
g->SetMarkerColor(color);
g->SetLineColor(color);

}

