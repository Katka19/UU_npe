#define UUVSAUAUYIELD 0

void customize_histo(TH1F *h, int marker, int color);
void customize_graph(TGraphErrors *g, int marker, int color);
//--------------------------
void raa(){

gStyle->SetLegendBorderSize(0);
//gStyle->SetFillColor(-1);

TFile *yieldFile = TFile::Open("../final_plots/yield_file.root");
TFile *DanielFile = TFile::Open("../../Daniel_AuAu/Raa/nph_yield_Raa_central_trigger_only_0_5pct.root");
//TFile *DanielFile2 = TFile::Open("../../Daniel_AuAu/Raa/nph_yield_Raa.root");
TFile *efficiencyFile = TFile::Open("../efficiencies_with_errors.root");
TFile *relFile = TFile::Open("../UU_relative_stat_error.root");
TFile *ppCorrFile = TFile::Open("final_corr_pp.root");
TFile *ppCorRebFile = TFile::Open("final_corr_pp_rebin.root");
TFile *ppUncorFile = TFile::Open("pp09_uncorr_baseline.root");
TFile *yieldRebFile = TFile::Open("../rebinned_yield/yield_file_rebinned.root");
//TFile *auauFile = TFile::Open("../../for_xiaozhi/Data/nph_yield_Raa_central_trigger_only_0_5pct.root");

TH1F *relUU = (TH1F*)relFile->Get("relative_error")->Clone("relUU");
TH1F *uuyield = (TH1F*)yieldFile->Get("yield_npe")->Clone("uuyield");
TH1F *uuyield_syst = (TH1F*)yieldFile->Get("yield_npe_syst")->Clone("uuyield_syst");
TH1F *uuyieldReb = (TH1F*)yieldRebFile->Get("yield_npe")->Clone("uuyieldReb");
TH1F *uuyieldReb_syst = (TH1F*)yieldRebFile->Get("yield_npe_syst")->Clone("uuyieldReb_syst");
TH1F *DanielRaa = (TH1F*)DanielFile->Get("Raa");//..0-5% centrality
TH1F *DanielRaaSys = (TH1F*)DanielFile->Get("Raa_sys");//..0-5% centrality
TH1F *DanielRaaSyspp = (TH1F*)DanielFile->Get("Raa_sys_pp");
TCanvas *cpp = new TCanvas("cpp","cpp",600,450);
DanielRaaSyspp->SetMarkerStyle(29);
DanielRaaSyspp->SetMarkerColor(kBlack);
DanielRaaSyspp->Draw("EP");
cpp->SaveAs("DanielRaa_pperrors.pdf");
cout<<"bin = "<<DanielRaaSyspp->GetXaxis()->FindBin(2.6)<<endl;
cout<<"yield bin = "<<DanielRaaSyspp->GetBinContent(29)<<endl;
//TH1F *DanielRaa0 = (TH1F*)DanielFile2->Get("Raa0");//..0-10% centrality
//TH1F *ppyield_syst09 = (TH1F*)efficiencyFile->Get("ppyield_systerr09")->Clone("ppyield_syst09");
TH1F *ppyield_syst09 = (TH1F*)efficiencyFile->Get("pp_spectra09")->Clone("ppyield_syst09");
//TH1F *ppyield_syst09 = (TH1F*)efficiencyFile->Get("comb_result")->Clone("ppyield_syst09");//....in reality it is combined 08 and 09 from Daniel, i am just lazy :)
//TH1F *ppyield_syst08 = (TH1F*)efficiencyFile->Get("ppyield_systerr08")->Clone("ppyield_syst08");//...in reality it is run08 only
TH1F *ppyield_syst08 = (TH1F*)efficiencyFile->Get("official_pp")->Clone("ppyield_syst08");
//TH1F *ppyield_syst08 = (TH1F*)efficiencyFile->Get("h_pp")->Clone("ppyield_syst08");
TH1F *ppyield0508 = (TH1F*)efficiencyFile->Get("h_pp")->Clone("ppyield0508");
TF1 *fit_pp08 = (TF1*)efficiencyFile->Get("f")->Clone("fit_pp08");
//TH1F *ppyield_syst09 = (TH1F*)efficiencyFile->Get("ppyield_systerr05")->Clone("ppyield_syst09");//...in reality it is phenix run05
//TH1F *pp_spectra08 = (TH1F*)efficiencyFile->Get("pp_spectra")->Clone("pp_spectra08");
TH1F *ppCor = (TH1F*)ppCorrFile->Get("final_correction")->Clone("ppCor");
TH1F *ppCorReb = (TH1F*)ppCorRebFile->Get("final_correction")->Clone("ppCorReb");
TH1F *auauYield = (TH1F*)DanielFile->Get("nph_rebin")->Clone("auauYield");

Double_t Rebin_array[15] = {1.2,1.4,1.6,1.8,2,2.4,2.8,3.2,3.6,4,4.4,4.8,5.2,5.6,6};
double Center_array[14] = {1.29366, 1.49413, 1.69453, 1.89489, 2.18142, 2.58337, 2.98495, 3.38625, 3.78735, 4.18828, 4.58909, 4.98979, 5.39041, 5.79096}; //..from bin shift correction
int n_rebin = 14;

#if UUVSAUAUYIELD

cout<<"bin = "<<ppyield_syst08->FindBin(2.75)<<endl;//34
//double bin = auauYield->GetNbinsX();
double auX[6] = {0};
double auY[6] = {0};
TH1F *auG = new TH1F("auG","auG",7,2.5,6);
TH1F *auRaa = new TH1F("auRaa","auRaa",7,2.5,6);
double auBin[8] = {2.5,3,3.5,4,4.5,5,5.5,6};
for(int i=29; i<36; i++){
	//auX[i-29] = auauYield->GetXaxis()->GetBinCenter(i);
	//auY[i-29] = auauYield->GetBinContent(i);
	auG->SetBinContent(i-28, auauYield->GetBinContent(i));
cout<<"pT = "<<ppyield_syst08->GetXaxis()->GetBinCenter(i)<<endl;
	//double pp = fit_pp08->Eval((auBin[i-29] + auBin[i-28])/2);
	double pp = ppyield_syst08->GetBinContent(i);
//cout<<"pp pT = "<<(auBin[i-29] + auBin[i-28])/2<<endl;
cout<<"pp = "<<pp<<endl;
	double au = auauYield->GetBinContent(i);
cout<<"au = "<<au<<endl;
	double ratio = au/pp;
cout<<"atio = "<<ratio<<endl;
	double ra = ratio;//..au spectra were already normalized to nbinary
cout<<"au raa = "<<ra<<endl;
	auRaa->SetBinContent(i-28, ra);
}
TCanvas *aucan = new TCanvas("aucan","aucan",600,450);
auRaa->SetMarkerStyle(29);
auRaa->SetMaximum(1.4);
auRaa->SetMinimum(0);
auRaa->Draw("P");
aucan->SaveAs("auRaa.pdf");
//TGraph *auG = new TGraph(6,auX,auY);

TH1F *uuYield = (TH1F*)uuyield->Rebin(7,"uuYield",auBin);
aucan->cd();
uuYield->SetMarkerStyle(29);
uuYield->Draw("P");
aucan->SetLogy();
aucan->SaveAs("rebinned_uuyield.pdf");

TF1 *f = new TF1("f","[0]*TMath::Power(x+[1],-[2])",2.4,6);
f->SetParameter(0,0.001);
f->SetParameter(1,0.3);
f->SetParameter(2,7.5);
double Bbin_Center[9] = {2.6,3,3.4,3.8,4.2,4.6,5,5.4,5.8};
TCanvas *c = new TCanvas("c","c",600,450);
gStyle->SetOptFit(1111);
double error[9] = {0};
TFitResultPtr r = auG->Fit(f,"IS");
r->GetConfidenceIntervals(9,1,0,Bbin_Center,error,0.95);
auG->SetMarkerStyle(29);
auG->Draw("P");
f->Draw("same");
c->SetLogy();
c->SaveAs("auau_yield.pdf");

double au_to_uu[9]={0};


double N_aubin = 1048;
double Bbin[10] = {2.4,2.8,3.2,3.6,4,4.4,4.8,5.2,5.6,6};
//double Bbin_Center[9] = {2.6,3,3.4,3.8,4.2,4.6,5,5.4,5.8};
int nbbin = 9;
TH1F *h = new TH1F("h","h;p_{T};uu/au yield",nbbin,Bbin);
for(int i=0; i<nbbin; i++){
double au = (f->Eval((Bbin[i] + Bbin[i+1])/2));
double erE = error[i];
//double au = auG->GetBinContent(i+1);
double uuu = uuyield->GetBinContent(i+6);
double uuE = uuyield->GetBinError(i+6);
double uuSe = uuyield_syst->GetBinError(i+6);
double uuErr = TMath::Sqrt((uuE)*(uuE) + (uuSe)*(uuSe));
//cout<<" au bin = "<<(Bbin[i] + Bbin[i+1])/2<<endl;
//cout<<" uu bin = "<<Rebin_array[i+5]<<endl;
//cout<<"au = "<<au<<endl;
//cout<<"uu = "<<uu[i+5]<<endl;
au_to_uu[i] = (uuu/au)*(1./1281);
double finalE = au_to_uu[i]*TMath::Sqrt((uuErr/uuu)*(uuErr/uuu) + (erE/au)*(erE/au));
h->SetBinContent(i+1, au_to_uu[i]);
h->SetBinError(i+1,finalE);
}

TCanvas *cc = new TCanvas("cc","cc",600,450);
//TGraphErrors *g = new TGraphErrors(nbbin,Bbin_Center,au_to_uu,0,0);
h->SetTitle("UU yield / AuAu yield");
h->SetStats(0);
//h->SetMaximum(1);
h->SetMinimum(0);
h->SetMarkerStyle(29);
h->Draw("EP");
cc->SaveAs("uu_to_auau.pdf");

#endif



TH1F *relraaUU = new TH1F("relraaUU","relraaUU",n_rebin, Rebin_array);


//double N_bin = 1341; //..from the paper
//double N_bin = 1281;//..from Zhenyu (Hiroshi)
double N_bin = 1156.166;
double N_bin_error = 105;//..from the paper
double N_bin_systError = 125.094;//..difference between correction and Hitoshi's value
double N_binE = TMath::Sqrt((N_bin_error)*(N_bin_error) + (N_bin_systError)*(N_bin_systError));
cout<<"N_binE = "<<N_binE<<endl;

double raa[17] = {0};
double raa08[17] = {0};
double raa_error_uustat[17] = {0};
double raa_error_uustat08[17] = {0};
double raa_error_uusyst[17] = {0};
double raa_error_uusyst08[17] = {0};
double raa_error_ppsyst[17] = {0};
double raa_error_ppsyst08[17] = {0};

double raa_error_UUpp[17] = {0};
double raa_error_UUpp08[17] = {0};

double pp_fit[17] = {0};
double pp08[17] = {0};
double pp_corr[17] = {0};
double pp_error[17] = {0};
double pp_error08[17] = {0};

double uu[17] = {0};
double uu_error_stat[17] = {0};
double uu_error_syst[17] = {0};

double x_error[17] = {0};
double x_error_pp[17] = {0};

double bin_center[15] = {0};

double rel_err_uu[15] = {0};
double rel_err_pp[15] = {0};
double rel_err_uu_syst[15] = {0};

TH1F *ppYield_stat = (TH1F*)efficiencyFile->Get("pp_stat")->Clone("ppYield_stat");
TH1F *ppYield_syst = (TH1F*)efficiencyFile->Get("pp_syst")->Clone("ppYield_syst");

//..Raa from U+U collisions
for(int i=0;i<n_rebin;i++){
	double bin = ppyield_syst09->GetXaxis()->FindBin(Center_array[i]);
	pp_fit[i] = (ppyield_syst09->GetBinContent(bin));
	pp08[i] = (ppyield_syst08->GetBinContent(i+1));
	//pp_corr[i] = ppCor->Eval((Rebin_array[i] + Rebin_array[i+1])/2);
	//pp_corr[i] = ppCor->Eval(Center_array[i]);
	pp_corr[i] = ppCor->GetBinContent(i+1);
	pp_fit[i] = pp_fit[i]*pp_corr[i];
	pp08[i] = pp08[i]*pp_corr[i];
	pp_error[i] = (ppyield_syst09->GetBinError(bin));
	pp_error[i] = pp_error[i]*pp_corr[i];
	//double e1 = ppyield_syst09->GetBinError(i+1)/ppyield_syst09->GetBinContent(i+1);
	//double e2 = ppCor->GetBinError(i+1)/pp_corr[i];
	//pp_error[i] = pp_fit[i]*TMath::Sqrt(e1*e1 + e2*e2);
	pp_error08[i] = (ppyield_syst08->GetBinError(i+1));
	pp_error08[i] = pp_error08[i]*pp_corr[i];
	//double ee1 = ppyield_syst08->GetBinError(i+1)/ppyield_syst08->GetBinContent(i+1);
	//pp_error08[i] = pp08[i]*TMath::Sqrt(ee1*ee1 + e2*e2);
	uu[i] = uuyield->GetBinContent(i+1);
	uu_error_stat[i] = uuyield->GetBinError(i+1);
	uu_error_syst[i] = uuyield_syst->GetBinError(i+1);

	
	raa[i] = (1/N_bin)*(uu[i]/(pp_fit[i]));
cout<<"raa = "<<raa[i]<<endl;
	raa08[i] = (1/N_bin)*(uu[i]/(pp08[i]));
	//*****************
	double binStat = ppYield_stat->GetXaxis()->FindBin(Center_array[i]);
	double ppstat = raa[i]*(ppYield_stat->GetBinError(binStat)/pp_fit[i]);
	double uustat = (1/N_bin)*(uu_error_stat[i]/pp_fit[i]);
	raa_error_uustat[i] = TMath::Sqrt(uustat*uustat + ppstat*ppstat);
	//*****************
	//raa_error_uustat[i] = (1/N_bin)*(uu_error_stat[i]/(pp_fit[i]));
	raa_error_uustat08[i] = (1/N_bin)*(uu_error_stat[i]/(pp08[i]));
	//******************
	double binSyst = ppYield_syst->GetXaxis()->FindBin(Center_array[i]);
	double ppsyst = raa[i]*(ppYield_syst->GetBinError(binSyst)/pp_fit[i]);
	double uusyst = (1/N_bin)*(uu_error_syst[i]/pp_fit[i]);
	raa_error_uusyst[i] = TMath::Sqrt(uusyst*uusyst + ppsyst*ppsyst);
	//******************
	//raa_error_uusyst[i] = (1/N_bin)*(uu_error_syst[i]/(pp_fit[i]));
	raa_error_uusyst08[i] = (1/N_bin)*(uu_error_syst[i]/(pp08[i]));
	//raa_error_ppsyst[i] = (1/N_bin)*(uu[i]/((pp_fit[i])*(pp_fit[i])))*(pp_error[i]);
	raa_error_ppsyst[i] = raa[i]*(pp_error[i]/pp_fit[i]);
	raa_error_ppsyst08[i] = raa08[i]*(pp_error08[i]/pp08[i]);

	raa_error_UUpp[i] = TMath::Sqrt((raa_error_uusyst[i])*(raa_error_uusyst[i]) + (raa_error_ppsyst[i])*(raa_error_ppsyst[i]));
	raa_error_UUpp08[i] = TMath::Sqrt((raa_error_uusyst08[i])*(raa_error_uusyst08[i]) + (raa_error_ppsyst08[i])*(raa_error_ppsyst08[i]));

	bin_center[i] = (Rebin_array[i]+Rebin_array[i+1])/2;

	x_error[i] = 0.05;
	x_error_pp[i] = 0.07;

	rel_err_pp[i] = (raa_error_ppsyst08[i]/raa08[i])*100;
	rel_err_uu[i] = (raa_error_uustat[i]/raa[i])*100;
	relraaUU->SetBinContent(i+1, rel_err_uu[i]);
	rel_err_uu_syst[i] = (raa_error_uusyst[i]/raa[i])*100;

	double nbin_error = (raa[i]/N_bin)*N_binE;
	cout<<"error from nbin = "<<nbin_error/raa[i]<<endl;

	//double au = f->Eval((Rebin_array[i] + Rebin_array[i+1])/2);
	//au_to_uu[i] = uu[i]/au;
}

//..rebinned UU yield to match AuAu binning
double Reb[8] = {2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0};
int nrebin = 7;

double ppReb[7], uuReb[7], ppRebE[7], uuRebStatE[7], uuRebSystE[7], raaReb[7], raaReb_statE[7], raaReb_systE[7], raaReb_ppE[7];
double rebX[7], rebppX[7], corReb[7];
double relReb[7];
double raaReb_errorUUpp[7];

for(int i=0; i<nrebin; i++){

	double bin = ppyield0508->GetXaxis()->FindBin((Reb[i]+Reb[i+1])/2);
cout<<"bin = "<<bin<<endl;
	ppReb[i] = ppyield0508->GetBinContent(bin);	
cout<<"pp yield = "<<ppReb[i]<<endl;
	ppRebE[i] = ppyield0508->GetBinError(bin);
	corReb[i] = ppCorReb->GetBinContent(i+1);
	ppReb[i] = ppReb[i]*corReb[i];
	ppRebE[i] = ppRebE[i]*corReb[i];

	uuReb[i] = uuyieldReb->GetBinContent(i+1);
	uuRebStatE[i] = uuyieldReb->GetBinError(i+1);
	uuRebSystE[i] = uuyieldReb_syst->GetBinError(i+1);

	raaReb[i] = (1/N_bin)*(uuReb[i]/ppReb[i]);
	raaReb_statE[i] = (1/N_bin)*(uuRebStatE[i]/ppReb[i]);
	raaReb_systE[i] = (1/N_bin)*(uuRebSystE[i]/ppReb[i]);
	raaReb_ppE[i] = raaReb[i]*(ppRebE[i]/ppReb[i]);
	
	raaReb_errorUUpp[i] = TMath::Sqrt((raaReb_systE[i])*(raaReb_systE[i]) + (raaReb_ppE[i])*(raaReb_ppE[i]));

	rebX[i] = 0.05;
	rebppX[i] = 0.07;

	double err = TMath::Sqrt(raaReb_systE[i]*raaReb_systE[i] + raaReb_ppE[i]*raaReb_ppE[i]);
	relReb[i] = (err/raaReb[i])*100;
}


//..do not use pp08 point under 2.5 GeV/c
for(int i=0; i<5; i++){
	raa08[i] = -10;
	raa_error_uustat08[i] = 0;
	raa_error_uusyst08[i] = 0;
	raa_error_ppsyst08[i] = 0;
	rel_err_pp[i] = -10;
}



cout<<"!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
cout<<"raa UU = "<<raa[5]<<endl;
cout<<"raa error UU = "<<raa_error_ppsyst[5]<<endl;

double daniel_raa[15] = {0};
double daniel_raa_stat[15] = {0};
double daniel_raa_syst[15] = {0};

double rel_err_auau[15] = {0};
double rel_err_auau_syst[15] = {0};

double bin_center_auau[7] = {2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75};

//..Raa from Au+Au collisions (from Daniel Kikola)
for(int i=0; i<7; i++){
	double bin = DanielRaa->GetXaxis()->FindBin(bin_center_auau[i]);
	daniel_raa[i] = DanielRaa->GetBinContent(bin);
	daniel_raa_stat[i] = DanielRaa->GetBinError(bin);
	daniel_raa_syst[i] = DanielRaaSys->GetBinError(bin);
	double nbinerr = daniel_raa[i]*(27.46719/1048.11384);
cout<<"nbinerr AuAu = "<<nbinerr/daniel_raa[i]<<endl;

	rel_err_auau[i] = (daniel_raa_stat[i]/daniel_raa[i])*100;
	rel_err_auau_syst[i] = (daniel_raa_syst[i]/daniel_raa[i])*100;
}

cout<<"!!!!!!!!!!!!!!!!!!!!"<<endl;
cout<<"raa AuAu = "<<daniel_raa[0]<<endl;
cout<<"raa error AuAu = "<<daniel_raa_syst[0]<<endl;

//..plot final Raa (U+U and Au+Au together)
TCanvas *graph_can = new TCanvas("graph_can","graph_can",600,450);
TGraphErrors *raa_uu_stat = new TGraphErrors(n_rebin, bin_center, raa, 0, raa_error_uustat);
TGraphErrors *raa_uu_syst = new TGraphErrors(n_rebin, bin_center, raa, x_error, raa_error_uusyst);// raa_error_UUpp);
TGraphErrors *raa_pp_syst = new TGraphErrors(n_rebin, bin_center, raa, x_error_pp, raa_error_ppsyst);
TGraphErrors *raa_daniel = new TGraphErrors(7, bin_center_auau, daniel_raa, 0, daniel_raa_stat);
TGraphErrors *raa_daniel_sys = new TGraphErrors(7, bin_center_auau, daniel_raa, x_error, daniel_raa_syst);
//TGraphErrors *raa_uu_stat08 = new TGraphErrors(n_rebin, bin_center, raa08, 0, raa_error_uustat08);
//TGraphErrors *raa_uu_syst08 = new TGraphErrors(n_rebin, bin_center, raa08, x_error, raa_error_uusyst08);
//TGraphErrors *raa_pp_syst08 = new TGraphErrors(n_rebin, bin_center, raa08, x_error_pp, raa_error_ppsyst08);
TGraphErrors *raa_uu_stat08 = new TGraphErrors(nrebin, bin_center_auau, raaReb, 0, raaReb_statE);
TGraphErrors *raa_uu_syst08 = new TGraphErrors(nrebin, bin_center_auau, raaReb, rebX, raaReb_errorUUpp);//raaReb_systE);
TGraphErrors *raa_pp_syst08 = new TGraphErrors(nrebin, bin_center_auau, raaReb, rebppX, raaReb_ppE);

TGraphErrors *cronin_upper = new TGraphErrors("eR-U_cron1_eloss0_S200D+Ball-Y0.0.75.evolve420New","%lg %lg");
TGraphErrors *cronin_lower = new TGraphErrors("eR-U_cron1_eloss0_S200D+Ball-Y0.1.evolve420New","%lg %lg");
TGraphErrors *eloss_upper = new TGraphErrors("eR-U_cron1.5_eloss1_S200D+Ball-Y0.0.75.evolve420New","%lg %lg");
TGraphErrors *eloss_lower = new TGraphErrors("eR-U_cron1.5_eloss1_S200D+Ball-Y0.1.evolve420New","%lg %lg");

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
customize_graph(raa_pp_syst, 29, kBlue+3);
//raa_pp_syst->SetFillColor(kWhite);
raa_pp_syst->SetFillStyle(0);
raa_uu_syst->SetFillColor(kBlue-9);
customize_graph(raa_daniel, 29, kRed+3);
customize_graph(raa_daniel_sys, 29, kRed+3);
raa_daniel_sys->SetFillColor(kRed-9);
//raa_daniel_sys->SetFillStyle(0);
customize_graph(raa_uu_stat08, 29, kGreen+3);
customize_graph(raa_uu_syst08, 29, kGreen+3);
customize_graph(raa_pp_syst08, 29, kGreen+3);
raa_pp_syst08->SetFillStyle(0);
raa_uu_syst08->SetFillColor(kGreen-9);

//raa_uu_syst->SetMaximum();
raa_uu_syst->Draw("AP2");
   //raa_pp_syst08->Draw("P2same");
raa_daniel_sys->Draw("P2same");
raa_uu_syst08->Draw("P2same");
raa_uu_syst->Draw("P2same");
   //raa_pp_syst->Draw("P2same");
//raa_daniel_sys->Draw("P2same");
raa_daniel->Draw("PZsame");
//raa_uu_syst08->Draw("P2same");
//raa_pp_syst08->Draw("P2same");
raa_uu_stat08->Draw("PZsame");
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

TLegend *legend = new TLegend(0.15,0.75,0.73,0.89);
legend->AddEntry(raa_uu_syst, "U+U #sqrt{s_{NN}} = 193 GeV (p+p 2009), centrality 0-5% (ToF+ZDC)","LPF");
legend->AddEntry(raa_uu_syst08, "U+U #sqrt{s_{NN}} = 193 GeV (p+p 2005+2008), centrality 0-5% (ToF+ZDC)","LPF");
legend->AddEntry(raa_daniel_sys, "Au+Au #sqrt{s_{NN}} = 200 GeV (p+p 2005+2008), centrality 0-5%","PLF");
legend->SetTextSize(0.025);
legend->SetFillColor(0);

line->Draw("same");
legend->Draw("same");

TLegend *leg_nbin = new TLegend(0.5, 0.53, 0.85, 0.65);
//leg_syst->AddEntry(raa_daniel_sys,"syst. uncertainty from pp08 baseline","F");
//leg_syst->AddEntry(raa_pp_syst, "syst. err. from pp09 baseline","F");
//leg_syst->AddEntry(raa_pp_syst08, "syst. err. from pp05+pp08 baseline","F");
leg_nbin->AddEntry(Nbin_err, "syst. uncertainty from N_{bin}^{UU}","F");
leg_nbin->AddEntry(Nbin_errAu, "syst. uncertainty from N_{bin}^{AuAu}","F");
leg_nbin->SetTextSize(0.025);
leg_nbin->SetFillColor(0);
leg_nbin->Draw("same");

TLegend *leg_syst = new TLegend(0.3, 0.65, 0.85, 0.75);
leg_syst->AddEntry(cronin_lower,"I.Vitev; (large Cronin, CNM e-loss)+(QGP e-loss, dissoc.)","L");
leg_syst->AddEntry(eloss_lower,"I.Vitev; (small Cronin, CNM e-loss)+(QGP e-loss, dissoc.)","L");
leg_syst->SetTextSize(0.025);
leg_syst->SetFillColor(0);
leg_syst->Draw("same");

//graph_can->SaveAs("../final_plots/raa_central05_run08AND09_Nbincorrected.pdf");
//graph_can->SaveAs("../final_plots/raa_central05_run08AND09_Nbincorrected.png");
//....................

//..again Raa but only run 09
graph_can->cd();
//TCanvas *ccc = new TCanvas("ccc","ccc",600,450);

raa_uu_syst->Draw("AP2");
raa_uu_stat->Draw("PZsame");
//raa_pp_syst->Draw("P2same");
raa_daniel_sys->Draw("P2same");
raa_daniel->Draw("Zsame");
line->Draw("same");
Nbin_err->Draw("2same");
Nbin_errAu->Draw("2same");

TLegend *l1 = new TLegend(0.15,0.75,0.73,0.85);
l1->AddEntry(raa_uu_syst, "U+U #sqrt{s_{NN}} = 193 GeV (p+p 2009), centrality 0-5% (ToF+ZDC)","LPF");
l1->AddEntry(raa_daniel_sys, "Au+Au #sqrt{s_{NN}} = 200 GeV (p+p 2005+2008), centrality 0-5%","PLF");
l1->SetTextSize(0.03);
l1->SetFillColor(0);
l1->Draw("same");

TLegend *l2 = new TLegend(0.5, 0.55, 0.85, 0.63);
//l2->AddEntry(raa_pp_syst, "syst. err. from pp09 baseline","F");
//l2->AddEntry(raa_daniel_sys,"syst. uncertainty from pp08 baseline","F");
l2->AddEntry(Nbin_err, "syst. err. from N_{bin}^{UU}","F");
l2->AddEntry(Nbin_errAu, "syst. uncertainty from N_{bin}^{AuAu}","F");
l2->SetTextSize(0.025);
l2->SetFillColor(0);

TLegend *lband = new TLegend(0.3, 0.65, 0.85, 0.75);
lband->AddEntry(cronin_lower,"I.Vitev; (large Cronin, CNM e-loss)+(QGP e-loss, dissoc.)","L");
lband->AddEntry(eloss_lower,"I.Vitev; (small Cronin, CNM e-loss)+(QGP e-loss, dissoc.)","L");
lband->SetTextSize(0.025);
lband->SetFillColor(0);

//cronin_upper->SetMaximum(2);
//cronin_upper->SetMinimum(0);
//cronin_upper->GetXaxis()->SetRangeUser(1,6);
cronin_upper->Draw("Lsame");
cronin_lower->Draw("Lsame");
eloss_upper->Draw("Lsame");
eloss_lower->Draw("Lsame");
l2->Draw("same");
lband->Draw("same");

graph_can->SaveAs("../final_plots/raa_central05_run09_separateUUErrors.pdf");
graph_can->SaveAs("../final_plots/raa_central05_run09_separateUUErrors.png");
//ccc->SaveAs("../final_plots/theoretical_bands.pdf");
//....................

//..plot relative statistical errors of U+U yield and U+U Raa
TCanvas *canvas = new TCanvas("canvas","canvas",600,450);
customize_histo(relUU, 29, kBlack);
relUU->SetMaximum(100);
relUU->GetYaxis()->SetTitle("relative statistical error [%]");
customize_histo(relraaUU, 29, kBlue+1);
relraaUU->SetMaximum(100);
//relraaUU->SetMarkerSize(1.5);
relraaUU->GetYaxis()->SetTitle("relative statistical error [%]");

relUU->Draw("P");
relraaUU->Draw("Psame");

TLegend *rel_lege = new TLegend(0.2, 0.65, 0.65, 0.85);
rel_lege->AddEntry((TObject*)0, "relative statistical errors of the UU yield and R_{AA}", "");
rel_lege->AddEntry(relUU, "U+U yield","p");
rel_lege->AddEntry(relraaUU, "U+U R_{AA}","p");
rel_lege->SetTextSize(0.035);
rel_lege->SetFillColor(0);
rel_lege->Draw("same");

//canvas->SaveAs("../final_plots/relative_stat_yield_raa_JPsi.pdf");
//.....................

//..plot relative statistical errors of Raa (U+U and Au+Au collisions)
TCanvas *rel_can = new TCanvas("rel_can","rel_can",600,450);
TGraphErrors *relative_error_uu = new TGraphErrors(n_rebin, bin_center, rel_err_uu, 0, 0);
TGraphErrors *relative_error_auau = new TGraphErrors(7, bin_center_auau, rel_err_auau, 0, 0);

customize_graph(relative_error_uu, 29, kBlue+1);
relative_error_uu->SetMaximum(100);
relative_error_uu->GetYaxis()->SetTitle("relative statistical error [%]");
customize_graph(relative_error_auau, 29, kRed+1);
relative_error_auau->SetMaximum(100);
relative_error_auau->GetYaxis()->SetTitle("relative statistical error [%]");

relative_error_uu->Draw("APZ");
relative_error_auau->Draw("PZsame");

TLegend *rel_lege = new TLegend(0.2, 0.65, 0.65, 0.85);
rel_lege->AddEntry((TObject*)0, "relative statistical errors of the R_{AA}", "");
rel_lege->AddEntry(relative_error_uu, "U+U collisions, centrality 0-5%","p");
rel_lege->AddEntry(relative_error_auau, "Au+Au collisions, centrality 0-5%","p");
rel_lege->SetTextSize(0.035);
rel_lege->SetFillColor(0);
rel_lege->Draw("same");

//rel_can->SaveAs("../final_plots/relative_errors_raa_UUemb_JPsi.pdf");
//.......................

//..plot relative systematical errors of Raa (U+U and Au+Au collisions)
TCanvas *rel_syst_can = new TCanvas("rel_syst_can","rel_syst_can",600,450);
//TGraphErrors *relative_syst_uu = new TGraphErrors(n_rebin, bin_center, rel_err_uu_syst, 0, 0);
//TGraphErrors *relative_syst_uu = new TGraphErrors(n_rebin, bin_center, rel_err_pp, 0, 0);
TGraphErrors *relative_syst_uu = new TGraphErrors(7, bin_center_auau, relReb, 0, 0);
TGraphErrors *relative_syst_auau = new TGraphErrors(7, bin_center_auau, rel_err_auau_syst, 0, 0);

customize_graph(relative_syst_uu, 29, kCyan+2);
relative_syst_uu->SetMaximum(100);
relative_syst_uu->GetYaxis()->SetTitle("relative systematical error [%]");
customize_graph(relative_syst_auau, 29, kMagenta+2);
relative_syst_auau->SetMaximum(100);
relative_syst_auau->GetYaxis()->SetTitle("relative systematical error [%]");

relative_syst_uu->Draw("APZ");
relative_syst_auau->Draw("PZsame");

TLegend *rel_syst_lege = new TLegend(0.3, 0.65, 0.65, 0.85);
rel_syst_lege->AddEntry((TObject*)0, "relative systematical errors of the R_{AA}", "");
rel_syst_lege->AddEntry(relative_syst_uu, "U+U collisions, centrality 0-5%","p");
rel_syst_lege->AddEntry(relative_syst_auau, "Au+Au collisions, centrality 0-5%","p");
rel_syst_lege->SetTextSize(0.035);
rel_syst_lege->SetFillColor(0);
rel_syst_lege->Draw("same");

//rel_syst_can->SaveAs("../final_plots/relative_syst_ppbaseline_rebinned.pdf");
//rel_syst_can->SaveAs("../final_plots/relative_syst_ppbaseline.pdf");
//........................

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

