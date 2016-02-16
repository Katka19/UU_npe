void curves(){

 TGraphErrors *cronin_upper = new TGraphErrors("eR-U_cron1_eloss0_S200D+Ball-Y0.0.75.evolve420New","%lg %lg");
 TGraphErrors *cronin_lower = new TGraphErrors("eR-U_cron1_eloss0_S200D+Ball-Y0.1.evolve420New","%lg %lg");
 TGraphErrors *eloss_upper = new TGraphErrors("eR-U_cron1.5_eloss1_S200D+Ball-Y0.0.75.evolve420New","%lg %lg");
 TGraphErrors *eloss_lower = new TGraphErrors("eR-U_cron1.5_eloss1_S200D+Ball-Y0.1.evolve420New","%lg %lg");

 TCanvas *c = new TCanvas("c","c",600,450);
 cronin_upper->SetLineColor(kRed);
 cronin_upper->SetMaximum(2);
 cronin_upper->SetMinimum(0);
 cronin_lower->SetLineColor(kRed);
 cronin_upper->Draw("AL");
 cronin_lower->Draw("Lsame");
 eloss_upper->Draw("Lsame");
 eloss_lower->Draw("Lsame");

 

}
