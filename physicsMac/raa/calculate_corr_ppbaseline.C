TGraphErrors *sigma_c, *sigma_b;//, *final_correction;
TH1F *final_correction;
TGraphErrors *sigma_c_up, *sigma_b_up;

const double corr_c = 0.98;
const double corr_c_err = 0.758;
const double corr_b = 0.917;
const double corr_b_err = 0.22;

const double NPoints;

/*double Range_array[15] = {1.2, 1.4, 1.6, 1.8, 2, 2.4, 2.8, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6};
//double Center_array[14] = {1.29342, 1.49393, 1.69437, 1.89476, 2.18104, 2.5831, 2.98477, 3.38613, 3.78727, 4.18824, 4.58907, 4.98979, 5.39043, 5.79099}; //..from bin shift correction
double Center_array[14] = {1.3, 1.5, 1.7, 1.9, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0, 5.4, 5.8};
const double nbin = 14;*/
/*double Range_array[8] = {2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0};
double Center_array[7] = {2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75};
const double nbin = 7;*/
double Range_array[17] = {1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0};
double Center_array[16] = {1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75};
const double nbin = 16;

//..functions
void GetUpperGraphLine(TGraphErrors *input, TGraphErrors *output);
//=====================================
void calculate_corr_ppbaseline(){

 //TFile *result = new TFile("final_corr_pp.root","recreate");
 TFile *result = new TFile("final_corr_pp_rebin_pp12.root","recreate");

 sigma_c = new TGraphErrors("cross_section_c_e.txt","%lg %lg %lg");
 sigma_b = new TGraphErrors("cross_section_b_e.txt","%lg %lg %lg");
 
 sigma_c_up = new TGraphErrors(nbin);
 sigma_b_up = new TGraphErrors(nbin);

 GetUpperGraphLine(sigma_c, sigma_c_up);
 GetUpperGraphLine(sigma_b, sigma_b_up);

 TCanvas *c = new TCanvas("c","c",600,450);
 sigma_c->SetLineColor(kBlue);
 sigma_c->SetFillColor(kBlue-9);
 sigma_c->GetXaxis()->SetTitle("p_{T} (GeV/c)");
 sigma_c->GetYaxis()->SetTitle("d#sigma/dp_{T}");
 sigma_c->SetTitle(0);
 sigma_c->Draw("A3");
 sigma_c->Draw("LXsame");

 sigma_b->SetLineColor(kRed);
 sigma_b->SetFillColor(kRed-9);
 sigma_b->Draw("3same");
 sigma_b->Draw("LXsame");

 gStyle->SetLegendBorderSize(0);
 TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.8);
 leg->AddEntry(sigma_c, "c -> e cross section","LF");
 leg->AddEntry(sigma_b, "b -> e cross section","LF");
 leg->SetFillColor(0);
 leg->Draw("same");

 c->SetLogy();
 //c->SaveAs("fonll_cross_sections.pdf");

 double xx[nbin];
 double yy[nbin];
 double error[nbin];

cout<<"sigma c 200 = "<<sigma_c->Eval(Center_array[0])<<endl;
cout<<"sigma c 200 error = "<<sigma_c_up->Eval(Center_array[0]) - sigma_c->Eval(Center_array[0])<<endl;
cout<<"sigma b 200 = "<<sigma_b->Eval(Center_array[0])<<endl;
cout<<"sigma b 200 error = "<<sigma_b_up->Eval(Center_array[0]) - sigma_b->Eval(Center_array[0])<<endl;

 final_correction = new TH1F("final_correction","final_correction",nbin, Range_array);

 for(int i=0; i<nbin; i++){

	//double c_200 = sigma_c->Eval((Range_array[i] + Range_array[i+1])/2);
	//double b_200 = sigma_b->Eval((Range_array[i] + Range_array[i+1])/2);
	double c_200 = sigma_c->Eval(Center_array[i]);
	double b_200 = sigma_b->Eval(Center_array[i]);

	double c_193 = c_200*corr_c;
	double b_193 = b_200*corr_b;

	//xx[i] = (Range_array[i] + Range_array[i+1])/2;
	xx[i] = Center_array[i];
	yy[i] = (c_193 + b_193)/(c_200 + b_200);

	//..calculate errors
	//double c_200_err = ((sigma_c_up->Eval((Range_array[i] + Range_array[i+1])/2)) - c_200);
	//double b_200_err = ((sigma_b_up->Eval((Range_array[i] + Range_array[i+1])/2)) - b_200);
	double c_200_err = ((sigma_c_up->Eval(Center_array[i])) - c_200);
        double b_200_err = ((sigma_b_up->Eval(Center_array[i])) - b_200);

	//double c_193_err = c_193*TMath::Sqrt((corr_c_err/corr_c)*(corr_c_err/corr_c) + (c_200_err/c_200)*(c_200_err/c_200));
	//double b_193_err = b_193*TMath::Sqrt((corr_b_err/corr_b)*(corr_b_err/corr_b) + (b_200_err/b_200)*(b_200_err/b_200));

	double num = c_193 + b_193;
	double denom = c_200 + b_200;
	
	//double num_err = TMath::Sqrt((c_193_err)*(c_193_err) + (b_193_err)*(b_193_err));
	//double denom_err = TMath::Sqrt((c_200_err)*(c_200_err) + (b_200_err)*(b_200_err));

	double e1 = (c_200/denom)*corr_c_err;
	double e2 = (b_200/denom)*corr_b_err;
	double e3 = ((corr_c*denom - num)/(denom*denom))*c_200_err;
	double e4 = ((corr_b*denom - num)/(denom*denom))*b_200_err;

	//error[i] = yy[i]*TMath::Sqrt((num_err/num)*(num_err/num) + (denom_err/denom)*(denom_err/denom));
	error[i] = TMath::Sqrt(e1*e1 + e2*e2 + e3*e3 + e4*e4);

	final_correction->SetBinContent(i+1, yy[i]);
	//final_correction->SetBinError(i+1, error[i]);
 }
 cout<<"final error = "<<error[0]<<endl;

 //final_correction = new TGraphErrors(nbin, xx, yy, 0, error);
 c->cd();
 final_correction->SetTitle("final correction of pp baseline");
 final_correction->GetXaxis()->SetTitle("p_{T} (GeV/c)");
 final_correction->SetMarkerStyle(29);
 final_correction->SetStats(0);
 final_correction->Draw("P");
 final_correction->SetMinimum(0);

 TLine *l = new TLine(1.2,1,6,1);
 l->SetLineColor(kBlack);
 l->SetLineStyle(2);
 l->Draw("same");

 c->SetLogy(0);
 //c->SaveAs("correction_ppbaseline.pdf");

 result->cd();
 final_correction->Write();

 result->Close();
 delete result;

}
//======================================
void GetUpperGraphLine(TGraphErrors *input, TGraphErrors *output){

 double xx[20];
 double yy[20];

 for(int i=0; i<20; i++){

	input->GetPoint(i, xx[i], yy[i]);
	double err = input->GetErrorY(i);

	yy[i] = yy[i] + err;

	output->SetPoint(i, xx[i], yy[i]);
	
 }

}
