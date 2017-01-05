{

  TFile * file = new TFile("TimeBin_1_Efficiency_2_200_phd.root");
  TGraph * g = (TGraph *)file->Get("g_efficiency");
  double x,y;

  for(int i=0; i<g->GetN(); i++) {
     g->GetPoint(i,x,y);
     cout << x << "," << y << endl;
  }

}
