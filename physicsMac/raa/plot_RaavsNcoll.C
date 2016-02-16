void plot_RaavsNcoll(){

 TFile *piFile = TFile::Open("Graph0_yellowline.root");
 TFile *blueFile = TFile::Open("Graph1_bluepoints.root");
 TFile *redFile = TFile::Open("statga1_redpoints.root");
 TFile *rederrFile = TFile::Open("systga1_rederrors.root");

 TFile *auFile = TFile::Open("../../Daniel_AuAu/Spectra_pho_npe_ratio/Data/nph_yield_Raa.root");
 TFile *auCentFile = TFile::Open("../../Daniel_AuAu/Spectra_pho_npe_ratio/Data/nph_yield_Raa_central_trigger_only_0_5pct.root");
 //TFile *uuFile = TFile::Open("../final_plots/yield_file.root");
 TFile *uuFile = TFile::Open("../rebinned_yield/yield_file_rebinned.root");
 //TFile *ppFile = TFile::Open("../nph_yield_pp_subtract_jpsi_DY.root");
 TFile *ppFile = TFile::Open("../rebinned_yield/pp12_yield.root");

 TH1F *auYield[6], *auYield_sys[6];
 
 double Nscale[6] = {1048.11384, 941.23714, 593.66913, 290.87634, 91.33495, 21.57039}; 

 for(int iAuCent=0; iAuCent<5; iAuCent++){
	auYield[iAuCent+1] = (TH1F*)auFile->Get(Form("nph%d",iAuCent))->Clone(Form("auYield%d",iAuCent));
	auYield[iAuCent+1]->Scale(Nscale[iAuCent+1]);
	auYield_sys[iAuCent+1] = (TH1F*)auFile->Get(Form("nph_sys%d",iAuCent))->Clone(Form("auYield_sys%d",iAuCent));	
	auYield_sys[iAuCent+1]->Scale(Nscale[iAuCent+1]);	
 }

 auYield[0] = (TH1F*)auCentFile->Get("nph")->Clone("auYield0");
 auYield[0]->Scale(Nscale[0]);
 auYield_sys[0] = (TH1F*)auCentFile->Get("nph_sys")->Clone("auYield_sys0");
 auYield_sys[0]->Scale(Nscale[0]);

 TH1F *uuYield = (TH1F*)uuFile->Get("yield_npe")->Clone("uuYield");
 TH1F *uuYield_sys = (TH1F*)uuFile->Get("yield_npe_syst")->Clone("uuYield_sys");

 TH1F *ppYield = (TH1F*)ppFile->Get("nph_pp")->Clone("ppYield");

 //do correction to 193 GeV
 double array[7] = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0};
 double nbin = 6;
 TH1F *ppYieldCorrected = new TH1F("ppYieldCorrected","ppYieldCorrected",nbin,array);
 TFile *corrFile = TFile::Open("final_corr_pp_rebin.root");
 TH1F *ppCorr = (TH1F*)corrFile->Get("final_correction")->Clone("ppCorr");
 
 for(int i=0; i<nbin; i++){
	double bin = ppYield->GetXaxis()->FindBin(array[i]+0.1);
	double content = ppYield->GetBinContent(bin);
	double error = ppYield->GetBinError(bin);
	
	double bincorr = ppCorr->GetXaxis()->FindBin(array[i]+0.1);
	double correctionContent = ppCorr->GetBinContent(bincorr);
cout<<"correction = "<<correctionContent<<endl;
	
	content = content*correctionContent;
	error = error*correctionContent;

	ppYieldCorrected->SetBinContent(i+1, content);
	ppYieldCorrected->SetBinError(i+1, error);
 }

 TGraphErrors *piGraph = (TGraphErrors*)piFile->Get("Graph0")->Clone("piGraph");
 TGraphErrors *blueGraph = (TGraphErrors*)blueFile->Get("Graph1")->Clone("blueGraph");
 TGraphErrors *redGraph = (TGraphErrors*)redFile->Get("statga1")->Clone("redGraph");
 TGraphErrors *redErrGraph = (TGraphErrors*)rederrFile->Get("systga1")->Clone("redErrGraph");
 //double xx[4] = {126.7, 325.5, 174.1, 41.8};
 //double xxerr[4] = {7.7, 3.6, 9.9, 7.8};//7.7 MB
 double xx[3] = {325.5, 174.1, 41.8};
 double xxerr[3] = {3.6, 9.9, 7.8};
 //double yy[4] = {0.768, 0.499, 0.581, 0.947};
 //double yyerr[4] = {0.104, 0.0579, 0.109, 0.249};
 double yy[3] = {0.499, 0.581, 0.947};
 double yyerr[3] = {0.0579, 0.109, 0.249};
 TGraphErrors *blueErrGraph = new TGraphErrors(3, xx, yy, xxerr, yyerr);

 double parterr[4] = {14,14,14,13};
 //double parterr[4] = {10,10,10,10};
 for(int i=0; i<4; i++){
	if(i==0){redErrGraph->SetPoint(i,-5, -5 );}
	else{
		double yyErr = redErrGraph->GetErrorY(i);
		redErrGraph->SetPointError(i,parterr[i],yyErr);}
 }

 double bluex[1], bluey[1];
 blueGraph->GetPoint(0,bluex[0], bluey[0]);
 double blueerry[1] = {blueGraph->GetErrorY(0)};
 TGraphErrors *blueMB = new TGraphErrors(1,bluex,bluey,0,blueerry);
 double redx[1], redy[1];
 redGraph->GetPoint(0,redx[0],redy[0]);
 double rederry[1] = {redGraph->GetErrorY(0)};
 TGraphErrors *redMB = new TGraphErrors(1,redx,redy,0,rederry);

 blueGraph->SetPoint(0,-10,-10);
 redGraph->SetPoint(0,-10,-10);

 double auyieldInt[6], auyieldInt410[6];
 double uuyieldInt, ppYieldInt, ppYieldInt410, ppYieldCorrInt;
 double auStatErr[6], auSysErr[6], auStatErr410[6], auSysErr410[6];
 double uuStatErr, uuSysErr, ppErr, ppErr410, ppCorrErr;

 for(int iAuCent=0; iAuCent<6; iAuCent++){
	double start_bin = auYield[iAuCent]->GetXaxis()->FindBin(3.1);
	double start4_bin = auYield[iAuCent]->GetXaxis()->FindBin(4.1);
	double end_bin = auYield[iAuCent]->GetXaxis()->FindBin(5.9);
	double end4_bin = auYield[iAuCent]->GetXaxis()->FindBin(9.9);
	
	auyieldInt[iAuCent] = auYield[iAuCent]->IntegralAndError(start_bin, end_bin, auStatErr[iAuCent]);
	auyieldInt[iAuCent] = auYield_sys[iAuCent]->IntegralAndError(start_bin, end_bin, auSysErr[iAuCent]);
	
	auyieldInt410[iAuCent] = auYield[iAuCent]->IntegralAndError(start4_bin, end4_bin, auStatErr410[iAuCent]);
	auyieldInt410[iAuCent] = auYield_sys[iAuCent]->IntegralAndError(start4_bin, end4_bin, auSysErr410[iAuCent]);
 }

 double start_binpp = ppYield->GetXaxis()->FindBin(3.1);
 double end_binpp = ppYield->GetXaxis()->FindBin(5.9);
 double start4_binpp = ppYield->GetXaxis()->FindBin(4.1); 
 double end4_binpp = ppYield->GetXaxis()->FindBin(9.9);
 ppYieldInt = ppYield->IntegralAndError(start_binpp, end_binpp, ppErr);
cout<<"ppYieldInt = "<<ppYieldInt<<endl;
 double binCorr1 = ppYieldCorrected->GetXaxis()->FindBin(3.1);
 double binCorr2 = ppYieldCorrected->GetXaxis()->FindBin(5.9);
 ppYieldCorrInt = ppYieldCorrected->IntegralAndError(binCorr1,binCorr2,ppCorrErr);
 ppYieldInt410 = ppYield->IntegralAndError(start4_binpp, end4_binpp, ppErr410);

 double start_binUU = uuYield->GetXaxis()->FindBin(3.1);
 double end_binUU = uuYield->GetXaxis()->FindBin(5.9);
 uuyieldInt = uuYield->IntegralAndError(start_binUU, end_binUU, uuStatErr);
 uuyieldInt = uuYield_sys->IntegralAndError(start_binUU, end_binUU, uuSysErr);

 double Nbin[7] = {21.57039, 91.33495, 290.87634, 593.66913, 941.23714, 1048.11384,1156.166};//, 1155.9};
 double NbinErr[7] = {8.038872, 20.00852, 30.46602, 30.17927, 26.27357, 27.46719, 163.32};
 double Npart[7] = {21.03694, 62.44176, 142.68674, 237.26805, 325.47688, 349.64479, 382.022};//..Npart
 double NpartErr[7] = {6.15739, 10.01175, 10.67755, 8.51684, 3.62739, 2.10277, 33.3521};//..Npart error
 double Raa[7], Raa410[7] = {0};
 double RaaStat[7], RaaSys[7],RaaSyspp[7],RaaSysNbin[7] = {0};
 double RaaStat410[7], RaaSys410[7], RaaSyspp410[7] = {0};
 
 for(int iAuCent=0; iAuCent<6; iAuCent++){

 	Raa[iAuCent] = (1/Nbin[iAuCent])*(auyieldInt[5-iAuCent]/ppYieldInt);
	RaaStat[iAuCent] = (1/Nbin[iAuCent])*(auStatErr[5-iAuCent]/ppYieldInt);
	//Raa[iAuCent] = (auyieldInt[5-iAuCent]/ppYieldInt);
        //RaaStat[iAuCent] = (auStatErr[5-iAuCent]/ppYieldInt);
	//RaaSys[iAuCent] = (1/Nbin[iAuCent])*(auSysErr[5-iAuCent]/ppYieldInt);
//cout<<"Raa AuAu syst = "<<RaaSys[iAuCent]<<endl;
	RaaSys[iAuCent] = Raa[iAuCent]*TMath::Sqrt((NbinErr[iAuCent]/Nbin[iAuCent])*(NbinErr[iAuCent]/Nbin[iAuCent]) + (auSysErr[5-iAuCent]/auyieldInt[5-iAuCent])*(auSysErr[5-iAuCent]/auyieldInt[5-iAuCent]));
	RaaSyspp[iAuCent] = Raa[iAuCent]*(ppErr/ppYieldInt);
cout<<"Raa AuAu syst pp = "<<RaaSyspp[iAuCent]/Raa[iAuCent]<<endl;
	//RaaSysNbin[iAuCent] = Raa[iAuCent]*(NbinErr[iAuCent]/Nbin[iAuCent]);

	Raa410[iAuCent] = (1/Nbin[iAuCent])*(auyieldInt410[5-iAuCent]/ppYieldInt410);
	RaaStat410[iAuCent] = (1/Nbin[iAuCent])*(auStatErr410[5-iAuCent]/ppYieldInt410);
	RaaSys410[iAuCent] = (1/Nbin[iAuCent])*(auSysErr410[5-iAuCent]/ppYieldInt410);
	RaaSyspp410[iAuCent] = Raa410[iAuCent]*(ppErr410/ppYieldInt410);

 }

 Raa[6] = -5;
 RaaStat[6] = 0;
 RaaSys[6] = 0;
 RaaSyspp[6] = 0;

 double NbinUU[1], RaaUU[1], RaaStatUU[1], RaaSysUU[1], RaaSysppUU[1], RaaSysNbinUU[1], uuxErr[1], uuxxErr[1], NbinE[1];
 double NpartUU[1] = {382.022};
 double NpartErrUU[1] = {33.3521};
 //NbinUU[0] = 1281;
 NbinUU[0] = 1156.166;
 NbinE[0] = 163.32;
 RaaUU[0] = (1/NbinUU[0])*(uuyieldInt/ppYieldCorrInt);
 RaaStatUU[0] = (1/NbinUU[0])*(uuStatErr/ppYieldCorrInt);
 //RaaSysUU[0] = (1/NbinUU[0])*(uuSysErr/ppYieldCorrInt);
 RaaSysUU[0] = RaaUU[0]*TMath::Sqrt((NbinE[0]/NbinUU[0])*(NbinE[0]/NbinUU[0]) + (uuSysErr/uuyieldInt)*(uuSysErr/uuyieldInt));
 RaaSysppUU[0] = RaaUU[0]*(ppCorrErr/ppYieldCorrInt);
cout<<"pp err UU = "<<RaaSysppUU[0]/RaaUU[0]<<endl;
 //RaaSysNbinUU[0] = RaaUU[0]*(NbinErr[6]/NbinUU[0]);
 uuxErr[0] = 33.3521;
 //uuxErr[0] = 10;
 uuxxErr[0] = 10;

 double xErr[7] = {6.15739, 10.01175, 10.67755, 8.51684, 3.62739, 2.10277, 0};
 //double xErr[7] = {7, 7, 7, 7, 7, 7, 0};
 //double xErr[7] = {10,10,10,10,10,10,0};
 double xErrpp[7] = {10, 10, 10, 10, 10, 10, 0};
 TGraphErrors *g = new TGraphErrors(7, Npart, Raa, 0, RaaStat);
 TGraphErrors *gsys = new TGraphErrors(7, Npart, Raa, NpartErr, RaaSys);
 TGraphErrors *gsyspp = new TGraphErrors(7, Npart, Raa, NpartErr, RaaSyspp);
 TGraphErrors *gsysNbin = new TGraphErrors(7, Npart, Raa, 0, RaaSysNbin);

 TGraphErrors *g410 = new TGraphErrors(7, Nbin, Raa410, 0, RaaStat410);
 TGraphErrors *gsys410 = new TGraphErrors(7, Nbin, Raa410, NbinErr, RaaSys410);
 TGraphErrors *gsyspp410 = new TGraphErrors(7, Nbin, Raa410, NbinErr, RaaSyspp410);

 TGraphErrors *gUU = new TGraphErrors(1, NpartUU, RaaUU, 0, RaaStatUU);
 TGraphErrors *gsysUU = new TGraphErrors(1, NpartUU, RaaUU, uuxErr, RaaSysUU);
 TGraphErrors *gsysppUU = new TGraphErrors(1, NpartUU, RaaUU, uuxxErr, RaaSysppUU);
 TGraphErrors *gsysNbinUU = new TGraphErrors(1, NpartUU, RaaUU, 0, RaaSysNbinUU);

 double XbandD[1] = {415};
 double YbandD[1] = {1};
 double XerrbandD[1] = {5};
 double YerrbandD[1] = {0.096};
 TGraphErrors *bandD = new TGraphErrors(1,XbandD,YbandD,XerrbandD,YerrbandD);
 bandD->SetFillColor(kGreen+3);
 double XbandAu[1] = {405};
 double YbandAu[1] = {1};
 double XerrbandAu[1] = {5};
 //double YerrbandAu[1] = {0.166};
 double YerrbandAu[1] = {0.111};
 TGraphErrors *bandAu = new TGraphErrors(1,XbandAu,YbandAu,XerrbandAu,YerrbandAu);
 bandAu->SetFillColor(kGreen+2);
 double XbandUU[1] = {395};
 double YbandUU[1] = {1};
 double XerrbandUU[1] = {5};
 //double YerrbandUU[1] = {0.085};
 double YerrbandUU[1] = {0.111};
 TGraphErrors *bandUU = new TGraphErrors(1,XbandUU,YbandUU,XerrbandUU,YerrbandUU);
 bandUU->SetFillColor(kGreen+2);

 gsys->GetXaxis()->SetRangeUser(-5,1500);
 gsys->GetXaxis()->SetTitle("N_{part}");
 //gsys->GetXaxis()->SetTitle("N_{coll}");
 gsys->GetXaxis()->CenterTitle(true);
 gsys->GetXaxis()->SetTitleSize(0.045);
 gsys->GetXaxis()->SetLabelSize(0.04);
 gsys->GetYaxis()->SetTitle("R_{AA}");
 gsys->GetYaxis()->CenterTitle(true);
 gsys->GetYaxis()->SetTitleSize(0.045);
 gsys->GetYaxis()->SetLabelSize(0.04);
 gsys->SetTitle(0);
 gsys->SetMaximum(2.5);
 //gsys->SetMaximum(1.2);
 gsys->SetMinimum(0);
 gsys->SetFillColor(kRed-7);
 gsys->SetMarkerStyle(33);
 gsys->SetMarkerSize(1.3);
 gsys->SetMarkerColor(kRed+2);
 gsys->SetLineColor(kRed+2);
 gsys->SetLineWidth(1);

 gsyspp->GetXaxis()->SetRangeUser(-5,1500);
 gsyspp->SetTitle(0);
 gsyspp->SetLineColor(kRed+2);
 gsyspp->SetMarkerStyle(33);
 gsyspp->SetMarkerSize(1.3);
 gsyspp->SetMarkerColor(kRed+2);
 gsyspp->SetFillStyle(0);

 gsysNbin->GetXaxis()->SetRangeUser(-5,1500);
 gsysNbin->SetTitle(0);
 gsysNbin->SetLineColor(kRed+2);
 gsysNbin->SetMarkerStyle(29);
 gsysNbin->SetMarkerColor(kRed+2);

 g->GetXaxis()->SetRangeUser(-5,1500);
 g->SetTitle(0);
 g->SetMarkerStyle(33);
 g->SetMarkerSize(1.5);
 g->SetMarkerColor(kRed+2);
 g->SetLineColor(kRed+2);
 g->SetLineWidth(2);

 g410->SetMarkerStyle(33);
 g410->SetMarkerColor(kBlue);
 g410->SetLineColor(kBlue);

 gsys410->SetMarkerStyle(33);
 gsys410->SetMarkerColor(kBlue);
 gsys410->SetLineColor(kBlue);

 gsyspp410->SetMarkerStyle(33);
 gsyspp410->SetMarkerColor(kBlue);
 gsyspp410->SetLineColor(kBlue);
 gsyspp410->SetFillStyle(0);

 gUU->GetXaxis()->SetRangeUser(-5,1500);
 gUU->SetTitle(0);
 gUU->SetLineColor(kBlue+2);
 gUU->SetLineWidth(2);
 gUU->SetMarkerStyle(29);
 gUU->SetMarkerSize(1.5);
 gUU->SetMarkerColor(kBlue+2);
 
 gsysUU->GetXaxis()->SetRangeUser(-5,1500);
 gsysUU->SetTitle(0);
 gsysUU->SetFillColor(kBlue-7);
 gsysUU->SetMarkerStyle(29);
 gsysUU->SetMarkerSize(1.5);
 gsysUU->SetMarkerColor(kBlue+2);
 gsysUU->SetLineColor(kBlue+2);

 gsysppUU->GetXaxis()->SetRangeUser(-5,1500);
 gsysppUU->SetTitle(0);
 gsysppUU->SetLineColor(kBlue+2);
 gsysppUU->SetMarkerStyle(29);
 gsysppUU->SetMarkerColor(kBlue+2);
 gsysppUU->SetFillStyle(0);

 gsysNbinUU->GetXaxis()->SetRangeUser(-5,1500);
 gsysNbinUU->SetTitle(0);
 gsysNbinUU->SetLineColor(kBlue+2);
 gsysNbinUU->SetMarkerStyle(29);
 gsysNbinUU->SetMarkerColor(kBlue+2);

 blueGraph->SetMarkerStyle(24);//20
 blueGraph->SetMarkerSize(1);
 //blueGraph->SetMarkerColor(kGreen+1);
 //blueGraph->SetLineColor(kGreen+1);
 blueGraph->SetMarkerColor(kRed+2);
 blueGraph->SetLineColor(kRed+2);
 blueGraph->SetLineWidth(1);

 /*blueMB->SetMarkerStyle(24);
 blueMB->SetMarkerSize(1);
 //blueMB->SetMarkerColor(kGreen+1);
 blueMB->SetMarkerColor(kRed+2);
 blueMB->SetLineWidth(1);
 //blueMB->SetLineColor(kGreen+1);
 blueMB->SetLineColor(kRed+2);*/

 blueErrGraph->SetFillColor(kRed+2);
 blueErrGraph->SetMarkerStyle(24);//20
 //blueErrGraph->SetMarkerColor(kGreen+1);
 blueErrGraph->SetMarkerColor(kRed+2);
 blueErrGraph->SetMarkerSize(1);
 //blueErrGraph->SetLineColor(kGreen+1);
 blueErrGraph->SetLineColor(kRed+2);
 blueErrGraph->SetLineWidth(1);
 //blueErrGraph->SetFillStyle(3005);
 blueErrGraph->SetFillStyle(0);

 redGraph->SetMarkerStyle(25);//21
 redGraph->SetMarkerSize(1);
 //redGraph->SetMarkerColor(kOrange-7);
 //redGraph->SetLineColor(kOrange-7);
 redGraph->SetMarkerColor(kBlue+2);
 redGraph->SetLineColor(kBlue+2);
 redGraph->SetLineWidth(1);

 /*redMB->SetMarkerStyle(25);
 redMB->SetMarkerSize(1);
 //redMB->SetMarkerColor(kOrange-7);
 redMB->SetMarkerColor(kBlue+2);
 redMB->SetLineWidth(1);
 //redMB->SetLineColor(kOrange-7);
 redMB->SetLineColor(kBlue+2);*/

 //redErrGraph->SetFillColor(kOrange-9);
 redErrGraph->SetFillColor(kBlue+2);
 redErrGraph->SetMarkerStyle(25);//25
 //redErrGraph->SetMarkerColor(kOrange-7);
 redErrGraph->SetMarkerColor(kBlue+2);
 redErrGraph->SetMarkerSize(1);
 //redErrGraph->SetLineColor(kOrange-7);
 redErrGraph->SetLineColor(kBlue+2);
 redErrGraph->SetLineWidth(1);
 //redErrGraph->SetFillStyle(3005);
 redErrGraph->SetFillStyle(0);

 piGraph->SetFillColor(kGray+1);
 
 gStyle->SetLegendBorderSize(0);
 TLegend *leg = new TLegend(0.55, 0.65, 0.85, 0.85);
 leg->AddEntry((TObject*)0, "3 < p_{T} < 6 (GeV/c)", "");
 leg->SetFillColor(0);

 TLegend *leg2 = new TLegend(0.11, 0.6, 0.7, 0.88);
 leg2->AddEntry(gsys,"NPE Au+Au, |#eta|<0.5, 3<p_{T}<6 GeV/#font[12]{c}","lfp");
 leg2->AddEntry(gsysUU,"NPE U+U, |#eta|<0.7, 3<p_{T}<6 GeV/#font[12]{c}","lpf");
 //leg2->AddEntry(gsyspp,"syst. uncertainty from pp baseline","f");
 //leg2->AddEntry(gsysppUU,"syst. uncertainty from pp baseline","f");
 leg2->AddEntry(blueErrGraph,"D^{0} Au+Au, |y|<1, 3<p_{T}<8 GeV/#font[12]{c}, PRL113, 142301(2014)","lfp");
 leg2->AddEntry(redErrGraph,"D^{0} U+U, |y|<1, 3<p_{T}<5 GeV/#font[12]{c}","lpf");
 leg2->AddEntry(piGraph,"#pi^{#pm} Au+Au, |y|<0.5, p_{T}>6 GeV/#font[12]{c}, PLB655, 104(2007)","f");
 leg2->SetFillColor(0);
 leg2->SetTextSize(0.037);

 TLegend *leg3 = new TLegend(0.35, 0.5, 0.85, 0.6);
 leg3->AddEntry(bandD,"systematic err. from p+p baseline, D^{0}","f"); 
 leg3->AddEntry(bandAu,"systematic err. from p+p basline, NPE","f");
 leg3->SetFillColor(0);
 leg3->SetTextSize(0.037);
 
 TLine *l = new TLine(0, 1, 410, 1);
 //TLine *l = new TLine(0,1,1200,1);
 l->SetLineStyle(2);
 l->SetLineColor(1);
 l->SetLineWidth(1);

 TCanvas *c = new TCanvas("c","c",600,450);
 gsys->GetXaxis()->SetRangeUser(0,430);
 gsys->Draw("A2P");
 /*gsys->Draw("A[]P");
 gsyspp->Draw("2Psame");

 gsyspp410->Draw("2Psame");
 gsys410->Draw("[]Psame");
 g410->Draw("PZsame");
*/
 piGraph->Draw("e3same");
 blueErrGraph->Draw("2same");
 redErrGraph->Draw("2same");
 blueGraph->Draw("pZsame");
 redGraph->Draw("pZsame");
 gsys->Draw("2same");
 //blueMB->Draw("pZsame");
 //redMB->Draw("pZsame");

 g->Draw("PZsame");
 gsysUU->Draw("2same");
 gUU->Draw("PZsame");
 leg2->Draw("same");
 leg3->Draw("same");
 l->Draw("same");

 bandD->Draw("2same");
 bandAu->Draw("2same");
 
 TText *text = new TText();
 text->SetNDC();
 text->SetTextFont(1);
 text->SetTextColor(kRed);
 text->SetTextSize(0.04);
 text->SetTextAlign(22);
 text->SetTextAngle(0);
 text->DrawText(0.23, 0.55, "STAR Preliminary");

 c->SaveAs("RaavsNpart_pp12_preliminary.pdf");
 c->SaveAs("RaavsNpart_pp12_preliminary.png");
 //c->SaveAs("Raa_check_4_10.pdf");
}
