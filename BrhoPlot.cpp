
//-------------------TXT FILE----------------------
void graph40() {
  TCanvas *c1 = new TCanvas();

  TGraph* g = new TGraph("NaF_Na40.txt");
  g->SetTitle("23Na States - 40");
  g->Draw("A*");
  g->GetXaxis()->SetTitle("Position");
  g->GetYaxis()->SetTitle("Bp");
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->CenterTitle();

  TGraph* g2 = new TGraph("NaF_Na40.txt","%*lg %*lg %lg %lg");
  g2->Draw("*");
  g2->SetMarkerColor(3);

  g->Fit("pol2");

    int npoints = 450; //declaring new variable for number of points
    TGraphErrors *grint = new TGraphErrors(npoints); //creating new TGraphError
    for(int i=0;i<npoints;i++) //adding points
      grint->SetPoint(grint->GetN(),6*i,0);

    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,0.6827);
    grint->SetLineColor(kGreen);
    grint->SetFillColor(4);
    //grint->SetFillStyle(3010);
    grint->Draw("A3");
    g->Draw("same P");
}

void graph50() {
  TCanvas *c1 = new TCanvas();

  TGraph* g = new TGraph("NaF_Na50.txt");
  g->SetTitle("23Na States - 50");
  g->Draw("A*");
  g->GetXaxis()->SetTitle("Position");
  g->GetYaxis()->SetTitle("Bp");
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->CenterTitle();

  TGraph* g2 = new TGraph("NaF_Na50.txt","%*lg %*lg %lg %lg");
  g2->Draw("*");
  g2->SetMarkerColor(3);

  g->Fit("pol2");

    int npoints = 450; //declaring new variable for number of points
    TGraphErrors *grint = new TGraphErrors(npoints); //creating new TGraphError
    for(int i=0;i<npoints;i++) //adding points
      grint->SetPoint(grint->GetN(),6*i,0);

    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,0.6827);
    grint->SetLineColor(kGreen);
    grint->SetFillColor(4);
    //grint->SetFillStyle(3010);
    grint->Draw("A3");
    g->Draw("same P");
}

void graph35() {
  TCanvas *c1 = new TCanvas();

  TGraph* g = new TGraph("NaF_Na35.txt");
  g->SetTitle("23Na States - 35");
  g->Draw("A*");
  g->GetXaxis()->SetTitle("Position");
  g->GetYaxis()->SetTitle("Bp");
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->CenterTitle();

  TGraph* g2 = new TGraph("NaF_Na35.txt","%*lg %*lg %lg %lg");
  g2->Draw("*");
  g2->SetMarkerColor(3);

  g->Fit("pol2");

    int npoints = 450; //declaring new variable for number of points
    TGraphErrors *grint = new TGraphErrors(npoints); //creating new TGraphError
    for(int i=0;i<npoints;i++) //adding points
      grint->SetPoint(grint->GetN(),6*i,0);

    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,0.6827);
    grint->SetLineColor(kGreen);
    grint->SetFillColor(4);
    //grint->SetFillStyle(3010);
    grint->Draw("A3");
    g->Draw("same P");
}

void graph25() {
  TCanvas *c1 = new TCanvas();

  TGraph* g = new TGraph("NaF_Na25.txt");
  g->SetTitle("23Na States - 25");
  g->Draw("A*");
  g->GetXaxis()->SetTitle("Position");
  g->GetYaxis()->SetTitle("Bp");
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->CenterTitle();

  TGraph* g2 = new TGraph("NaF_Na25.txt","%*lg %*lg %lg %lg");
  g2->Draw("*");
  g2->SetMarkerColor(3);

  g->Fit("pol2");

    int npoints = 450; //declaring new variable for number of points
    TGraphErrors *grint = new TGraphErrors(npoints); //creating new TGraphError
    for(int i=0;i<npoints;i++) //adding points
      grint->SetPoint(grint->GetN(),6*i,0);

    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,0.6827);
    grint->SetLineColor(kGreen);
    grint->SetFillColor(4);
    //grint->SetFillStyle(3010);
    grint->Draw("A3");
    g->Draw("same P");
}

void graph70() {
  TCanvas *c1 = new TCanvas();

  TGraph* g = new TGraph("NaF_Na70.txt");
  g->SetTitle("23Na States - 70");
  g->Draw("A*");
  g->GetXaxis()->SetTitle("Position");
  g->GetYaxis()->SetTitle("Bp");
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->CenterTitle();

  TGraph* g2 = new TGraph("NaF_Na70.txt","%*lg %*lg %lg %lg");
  g2->Draw("*");
  g2->SetMarkerColor(3);

  g->Fit("pol2");

    int npoints = 450; //declaring new variable for number of points
    TGraphErrors *grint = new TGraphErrors(npoints); //creating new TGraphError
    for(int i=0;i<npoints;i++) //adding points
      grint->SetPoint(grint->GetN(),6*i,0);

    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,0.6827);
    grint->SetLineColor(kGreen);
    grint->SetFillColor(4);
    //grint->SetFillStyle(3010);
    grint->Draw("A3");
    g->Draw("same P");
}
