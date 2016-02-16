double get_correctionFactorcc(){

 double xx[5] = {39, 62.4, 200, 2760, 5000};
 double yy[5] = {0.0317, 0.0758, 0.3628, 3.473, 5.418};

 TGraph *g = new TGraph(5);
 for(int i=0; i<5; i++){
	g->SetPoint(i,xx[i],yy[i]);
 }

 TCanvas *c = new TCanvas("c","c",600,450);
 g->SetMarkerStyle(20);
 g->Draw("AP");
 c->SetLogy();
 c->SetLogx();

 TF1 *f = new TF1("f","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0,5500);
 g->Fit(f);

 double corr = (f->Eval(193))/(f->Eval(200));

 return corr;

}
