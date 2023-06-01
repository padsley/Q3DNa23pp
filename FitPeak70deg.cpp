#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TH1.h>

void FitRigidities();
void Draw23NaExSpectrum();

double p0 = 3.00029e-01, p1 = 6.65746e-06, p2 = -5.30683e-10;//pos->Brho calibration

double BrhoToExNaF(double Brho, double mass, double theta);

double mNa = 21414.8345021;
double mF = 17696.9005009;
double mSi = 26060.3420706;
double mO = 14899.1686366;
double mLi = 6535.36582157;

TGraphErrors *gExErrors1sig;

//------------Functions to calculate Energy loss after an indicated target------------------//
double ELossNaF(double Ep, double Thickness)
{
  //declaring variables
  double p0 = -21.865; //-40.92 (previous values commented)
  double p1 = 19.975;  //49.45
  double p2 = 4.5031;  //10.89

  //calculating ranges
  double R1 = p0 + p1*Ep + p2*pow(Ep,2.);//initial range
  double R2 = R1 - Thickness;//final range

  //go the other way using the quadratic formula
  double c = p0 - R2;
  double Ep2 = (-p1 + sqrt(pow(p1,2)-4*c*p2))/(2*p2);

  return Ep2;

}

double ELossLiF(double Ep, double Thickness)
{
  //declaring variables
  double p0 = -19.316; //-46.67
  double p1 = 17.279;  //55.79
  double p2 = 4.2728;  //13.47

  //calculating ranges
  double R1 = p0 + p1*Ep + p2*pow(Ep,2.);//initial range
  double R2 = R1 - Thickness;//final range

  //go the other way using the quadratic formula
  double c = p0 - R2;
  double Ep2 = (-p1 + sqrt(pow(p1,2)-4*c*p2))/(2*p2);

  return Ep2;
}

double ELossSiO2(double Ep, double Thickness)
{
  //declaring variables
  double p0 = -22.902; //-16.7
  double p1 = 19.272;  //20.76
  double p2 = 4.0374;  //4.627

  //calculating ranges
  double R1 = p0 + p1*Ep + p2*pow(Ep,2.);//initial range
  double R2 = R1 - Thickness;//final range

  //go the other way using the quadratic formula
  double c = p0 - R2;
  double Ep2 = (-p1 + sqrt(pow(p1,2)-4*c*p2))/(2*p2);

  return Ep2;

}

double ELossC(double Ep, double Thickness)
{
  //declaring variables
  double p0 = -21.833; //-15.57
  double p1 = 18.04;   //17.64
  double p2 = 4.4693;  //4.44

  //calculating ranges
  double R1 = p0 + p1*Ep + p2*pow(Ep,2.); //initial range
  double R2 = R1 - Thickness; //final range

  //go the other way using the quadratic formula
  double c = p0 - R2;
  double Ep2 = (-p1 + sqrt(pow(p1,2)-4*c*p2))/(2*p2);

  return Ep2;
}

double TailedGaussian(double *x, double *pars)
{
    
    double result = 0;
    
    if(pars[2] < (x[0] - pars[0])/pars[1])//this is the Gaussian part
    {
        result += pars[3] * exp ( - 0.5 * pow(x[0] - pars[0],2.)/pow(pars[1],2.));
    
    }
    else //this is the tail part
    {
        result += pars[3] * exp ( 0.5 * pow(pars[2],2.) - pars[2] * (x[0] - pars[0]) / pars[1]);
    }
    result += pars[4] + pars[5] * x[0];//constant+linear background
    
    return result;
}


void FitPeak70deg()
{
    TCanvas *c1 = new TCanvas();
    
//     gROOT->ProcessLine(".x run49.C");
    
//     TH1F *histo = (TH1F*)gROOT->FindObjectAny("hPos__49");
    TFile *fNaF = TFile::Open("run049_an.root");
    TTree *tNaF = (TTree*)fNaF->Get("readout");
    tNaF->Draw("pos>>histoNaF(620,0,2500)","pos>0","");
    TH1F *histoNaF = (TH1F*)gROOT->FindObjectAny("histoNaF");
    histoNaF->Draw("E");
    
    TF1 *fitty = new TF1("fitty",TailedGaussian,600,750,6);
    
    std::cout << "Ex = 19F: 8864 keV" << std::endl;
    fitty->SetParameters(705,27.6,-3,156,101,0);
    fitty->FixParameter(5,0.);
    
    histoNaF->Fit(fitty,"BRME");
    
    c1->SaveAs("figures/FitPeak0_70.png");
    
    fitty->SetRange(2200,2350);
    fitty->SetParameters(2284,40.,-2.1,400,140,0);
//     fitty->SetParLimits(5,0.,0.);
    std::cout << "Ex = 19F: 8629 keV" << std::endl;
    histoNaF->Fit(fitty,"BRME");
    
    c1->SaveAs("figures/FitPeak1_70.png");
    
    fitty->SetRange(760,860);
    fitty->SetParameters(810,25,-3,250,100);
    std::cout << "Ex = 8975.3(7)-keV state in 23Na" << std::endl;
    histoNaF->Fit(fitty,"BRME");
    
    c1->SaveAs("figures/FitPeak2_70.png");
    
    fitty->SetRange(900,1050);
    fitty->SetParameters(975,23,-3,900,100);
    std::cout << "Ex = 8945.1(8)-keV state in 23Na" << std::endl;
    histoNaF->Fit(fitty,"BRME");
    
    c1->SaveAs("figures/FitPeak3_70.png");
    
    fitty->SetRange(1860,2000);
    fitty->SetParameters(1950,23,-3,300,75);
    std::cout << "Ex = 8798.7(8)-keV state in 23Na" << std::endl;
    histoNaF->Fit(fitty,"BRME");
    
    c1->SaveAs("figures/FitPeak4_70.png");
    
//     std::cout << histoNaF->GetNbinsX() << "\t" << histoNaF->GetBinLowEdge(1) << "\t" << histoNaF->GetBinLowEdge(histoNaF->GetNbinsX()+1) << std::endl;
//     TH1F *histoLiF = new TH1F("histoLiF","",625,0,2500);
    TFile *fLiF = TFile::Open("run051_an.root");
    TTree *tLiF = (TTree*)fLiF->Get("readout");
//     tLiF->Draw("pos>>histoLiF(625,0,2500)","pos>0","");
    TH1F *histoLiF = (TH1F*)gROOT->FindObjectAny("histoLiF");
   
    FitRigidities();
}

void FitRigidities()
{
    TGraphErrors *ge = new TGraphErrors();
    
    ge->SetPoint(ge->GetN(),9.82088e+02,0.306073);//8945.1(8)-keV state in 23Na
    ge->SetPointError(ge->GetN()-1,4.27897e-01,2.70931e-5);
    
    ge->SetPoint(ge->GetN(),8.09844e+02,0.305053);//8975.3(7)-keV state in 23Na
    ge->SetPointError(ge->GetN()-1,7.53423e-01,2.37919e-5);
    
    ge->SetPoint(ge->GetN(),2.28506e+03,0.312559);//8629(4)-keV state in 19F
    ge->SetPointError(ge->GetN()-1,7.43117e-01,1.32159e-04);
    
    ge->SetPoint(ge->GetN(),7.11899e+02,0.304729);//8864(4)-keV state in 19F
    ge->SetPointError(ge->GetN()-1,8.95086e-01,1.35611e-04);
    
    ge->SetPoint(ge->GetN(),1.94649e+03,0.31097);//8798.7-keV state in 23Na
    ge->SetPointError(ge->GetN()-1,8.41249e-01,2.66624e-05);
    
//     ge->SetPoint(ge->GetN(),2284.59,0.31178090);//19F 8650-keV state
//     ge->SetPoint(ge->GetN()-1,1260,0.30927948);//28Si 8905-keV state
    ge->SetMarkerStyle(8);
 
    TCanvas *c2 = new TCanvas();
    
    ge->Draw("AP");
    
    TF1 *fPol2 = new TF1("fPol2","[0] + [1]*x + [2]*x*x",0,2560);
    fPol2->SetParameters(0.3,4.e-6,-4.e-10);
    
    ge->Fit(fPol2,"BRME");
    
    int npoints = 2500;
    TGraphErrors *grint = new TGraphErrors(npoints);grint->SetName("g1sigma");
    TGraphErrors *grint2 = new TGraphErrors(npoints);grint2->SetName("g2sigma");
    TGraphErrors *grint3 = new TGraphErrors(npoints);grint2->SetName("g3sigma");
    
    for(int i=0;i<npoints;i++)
    {
      grint->SetPoint(grint->GetN(),i,0);
      grint2->SetPoint(grint2->GetN(),i,0);
      grint3->SetPoint(grint3->GetN(),i,0);
    }

    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,0.6827);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint2,0.9545);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint3,0.9943);
//     grint->SetLineColor(kGreen);
    grint->SetFillColorAlpha(kBlue,1.0);
    grint2->SetFillColorAlpha(kBlue,0.66);
    grint3->SetFillColorAlpha(kBlue,0.33);
    //grint->SetFillStyle(3010);
    
    c2->SetLeftMargin(0.136);
    grint3->Draw("A3");
    grint3->SetTitle("");
//     grint->Draw("same3");
//     grint2->Draw("same 3");
    
    ge->Draw("same P");
    
    
    
    grint3->GetXaxis()->SetTitle("Focal-plane position [a.u.]");
    grint3->GetXaxis()->CenterTitle();
    grint3->GetXaxis()->SetTitleSize(0.05);
    grint3->GetXaxis()->SetLabelSize(0.05);
    grint3->GetXaxis()->SetTitleOffset(0.9);
    grint3->GetYaxis()->SetTitle("B#rho [Tm]");
    grint3->GetYaxis()->CenterTitle();
    grint3->GetYaxis()->SetTitleOffset(1.43);
    grint3->GetYaxis()->SetTitleSize(0.05);
    grint3->GetYaxis()->SetLabelSize(0.05);    
    
    c2->Update();
    
    c2->SaveAs("BrhoCalibration70deg.eps");
    c2->SaveAs("BrhoCalibration70deg.pdf");
    c2->SaveAs("BrhoCalibration70deg.png");
    c2->SaveAs("BrhoCalibration70deg.C");
    
    TCanvas *c3 = new TCanvas();
    
    TGraphErrors *gExErrors = new TGraphErrors();
    gExErrors1sig = new TGraphErrors();gExErrors1sig->SetName("gExErrors1sig");
    TGraphErrors *gExErrors2sig = new TGraphErrors();gExErrors2sig->SetName("gExErrors2sig");
    TGraphErrors *gExErrors3sig = new TGraphErrors();gExErrors3sig->SetName("gExErrors3sig");
    
    TGraph *gEx3sig_ActualValues = new TGraph();gEx3sig_ActualValues->SetName("gEx3sig_ActualValues");
    
//     gROOT->ProcessLine(".L NaProtonKinematicsCalculationNew.cpp");
    
    for(int i=0;i<grint->GetN();i++)
    {
        double x = 0, y= 0, sig_y = 0;
        x = grint->GetPointX(i);
        y = grint->GetPointY(i);
        
//         (double Brho, double mass, double theta)
        
        gExErrors->SetPoint(i,y,0);
        
        x = grint->GetPointX(i);
        y = grint->GetPointY(i);
        sig_y = grint->GetErrorYhigh(i);
        
        gExErrors1sig->SetPoint(i,BrhoToExNaF(y,938.782980,70.),0);
        gExErrors1sig->SetPointError(i,0, 1000.*(BrhoToExNaF(y,938.782980,70.) - BrhoToExNaF(y + sig_y,938.782980,70.))); 
        
        
        x = grint2->GetPointX(i);
        y = grint2->GetPointY(i);
        sig_y = grint2->GetErrorYhigh(i);
        
        gExErrors2sig->SetPoint(i,BrhoToExNaF(y,938.782980,70.),0);
        gExErrors2sig->SetPointError(i,0, 1000.*(BrhoToExNaF(y,938.782980,70.) - BrhoToExNaF(y + sig_y,938.782980,70.))); 
        
        x = grint3->GetPointX(i);
        y = grint3->GetPointY(i);
        sig_y = grint3->GetErrorYhigh(i);
        
        gExErrors3sig->SetPoint(i,BrhoToExNaF(y,938.782980,70.),0);
        gExErrors3sig->SetPointError(i,0, 1000.*(BrhoToExNaF(y,938.782980,70.) - BrhoToExNaF(y + sig_y,938.782980,70.))); 
        
//         std::cout << "i = " << i << " y = " << y << "\t sig_y = " << sig_y << std::endl;
        
        gEx3sig_ActualValues->SetPoint(i,BrhoToExNaF(y,938.782980,70.),1000.*(BrhoToExNaF(y,938.782980,70.) - BrhoToExNaF(y + sig_y,938.782980,70.)));
    }
    
    
        gExErrors1sig->SetFillColorAlpha(kBlue,1.0);
        gExErrors2sig->SetFillColorAlpha(kBlue,0.66);
        gExErrors3sig->SetFillColorAlpha(kBlue,0.33);
        
        gExErrors3sig->Draw("A3");
        gExErrors2sig->Draw("same3");
        gExErrors1sig->Draw("same3");
    
        gExErrors3sig->GetXaxis()->SetTitle("E_{x} [MeV]");
        gExErrors3sig->GetYaxis()->SetTitle("#Delta(E_{x}) [keV]");
        gExErrors3sig->GetXaxis()->CenterTitle();
        gExErrors3sig->GetYaxis()->CenterTitle();
        gExErrors3sig->GetXaxis()->SetTitleSize(0.05);
        gExErrors3sig->GetXaxis()->SetLabelSize(0.05);
        gExErrors3sig->GetXaxis()->SetTitleOffset(0.9);
        gExErrors3sig->GetYaxis()->SetTitleSize(0.05);
        gExErrors3sig->GetYaxis()->SetLabelSize(0.05);
        gExErrors3sig->GetYaxis()->SetTitleOffset(0.75);
        
        std::vector<double> vExValues;
        vExValues.push_back(8.7987);
        vExValues.push_back(8.8208);
        vExValues.push_back(8827.1e-3);
        vExValues.push_back(8862e-3);
        vExValues.push_back(8894e-3);
        vExValues.push_back(8945.1e-3);
        vExValues.push_back(8946.8e-3);
        vExValues.push_back(8963.9e-3);
        vExValues.push_back(8975.3e-3);
        vExValues.push_back(8976.e-3);
        vExValues.push_back(9000.e-3);
        vExValues.push_back(9039.5e-3);
        vExValues.push_back(9042.6e-3);
        vExValues.push_back(9072e-3);
        vExValues.push_back(9101.5e-3);
        vExValues.push_back(9113.e-3);
        
        std::cout << "**************" << std::endl << std::endl << std::endl << std::endl << std::endl;
        
        for(unsigned int i=0;i<vExValues.size();i++)
        {
            std::cout << vExValues.at(i) << "\t" << gEx3sig_ActualValues->Eval(vExValues.at(i))*1.e-3 << std::endl;
        }
        
        std::cout << "**************" << std::endl << std::endl << std::endl << std::endl << std::endl;
        
        c3->SaveAs("ExUncertaintyBand70deg.eps");
        c3->SaveAs("ExUncertaintyBand70deg.pdf");
        c3->SaveAs("ExUncertaintyBand70deg.png");
        c3->SaveAs("ExUncertaintyBand70deg.C");
        Draw23NaExSpectrum();
}

void Draw23NaExSpectrum()
{
    TCanvas *c4 = new TCanvas();
    c4->SetLeftMargin(0.136);

    
//     gROOT->ProcessLine(".L NaProtonKinematicsCalculationNew.cpp");
    
    TFile *fNaF = TFile::Open("run049_an.root");
    TTree *tNaF = (TTree*)fNaF->Get("readout");
//     tNaF->Draw("pos>>histoNaF(620,0,2500)","pos>0","");
//     TH1F *histoNaF = (TH1F*)gROOT->FindObjectAny("histoNaF");
    
    tNaF->SetAlias("Brho",Form("%e + %e * pos + %e * pos*pos",p0,p1,p2));
    
    std::cout << "Brho = " << Form("%e + %e * pos + %e * pos*pos",p0,p1,p2) << std::endl;
    tNaF->Draw("BrhoToExNaF(Brho,938.782980,70.)>>hEx23Na(450,8.70,9.15)","pos>0 && pos<2500","");
    
    TH1F *hEx23Na = (TH1F*)gROOT->FindObjectAny("hEx23Na");
    hEx23Na->SetTitle("");
    hEx23Na->SetStats(0);
    hEx23Na->GetXaxis()->SetTitle("E_{x} [MeV]");
    hEx23Na->GetXaxis()->CenterTitle();
    hEx23Na->GetYaxis()->CenterTitle();
    hEx23Na->GetYaxis()->SetTitle(Form("Counts per %d keV",(int)(hEx23Na->GetBinWidth(1)*1000)));
    hEx23Na->SetLineColor(1);
    hEx23Na->GetYaxis()->SetTitleSize(0.05);
    hEx23Na->GetYaxis()->SetLabelSize(0.05);
    hEx23Na->GetYaxis()->SetTitleOffset(1.43);
    hEx23Na->GetXaxis()->SetTitleSize(0.05);
    hEx23Na->GetXaxis()->SetLabelSize(0.05);
    hEx23Na->GetXaxis()->SetTitleOffset(0.9);
    hEx23Na->GetYaxis()->SetRangeUser(1,1550);
    
    TFile *fLiF = TFile::Open("run051_an.root");
    TTree *tLiF = (TTree*)fLiF->Get("readout");
    tLiF->SetAlias("Brho",Form("%e + %e * pos + %e * pos*pos",p0,p1,p2));
    tLiF->Draw("BrhoToExNaF(Brho,938.782980,70.)>>hEx19F(450,8.70,9.15)","pos>0 && pos<2500","goff");
    
    TH1F *hEx19F = (TH1F*)gROOT->FindObjectAny("hEx19F");
    hEx19F->SetLineColor(2);
    hEx19F->SetFillColorAlpha(kRed,0.35);
    hEx19F->Scale(0.5);
    hEx19F->Draw("same hist");
    
    double Sp = 8.79410;
    
    double axisLength = 350; //in keV
    
    TGaxis *EcmAxis = new TGaxis(Sp,1550,Sp+1.e-3*axisLength,1550,0,axisLength,510,"L-");
    EcmAxis->SetTextFont(42);
    EcmAxis->SetLabelFont(42);
    EcmAxis->SetLabelSize(0.05);
    EcmAxis->SetTitle("E_{cm} [keV]");
    EcmAxis->SetTitleSize(0.05);
    EcmAxis->SetTitleOffset(0.9);
    EcmAxis->CenterTitle();
    EcmAxis->Draw();
    
    TLine *l1 = new TLine(8.7987,1,8.7987,1400);
    l1->Draw("same");
    TBox *b1 = new TBox(8.7987-0.8e-3,1,8.7987+0.8e-3,1400);
    b1->SetFillColorAlpha(kBlack,0.25);
    b1->Draw("same");
    
    double FitUncertaintyAtPoint = 0.;
    
    double dummy = 0.1;
    int ClosestPoint = 0;
    
    for(int i=0;i<gExErrors1sig->GetN();i++)
    {
        
        if(abs(gExErrors1sig->GetPointX(i)-8.7987)<dummy)
        {
//             std::cout << abs(gExErrors1sig->GetPointX(i)-8.7987) << "\t" << dummy << std::endl;;
            
//             std::cout << "found a new smaller value: " << i << std::endl;
            dummy = abs(gExErrors1sig->GetPointX(i)-8.7987);
            FitUncertaintyAtPoint = gExErrors1sig->GetErrorYhigh(i);
//             std::cout << "FitUncertaintyAtPoint: " << FitUncertaintyAtPoint << std::endl;
            ClosestPoint = i;
        }
    }
    TBox *b1_2 = new TBox(8.7987-FitUncertaintyAtPoint*1.e-3,1,8.7987+FitUncertaintyAtPoint*1.e-3,1400);
    b1_2->SetFillColorAlpha(kBlue,0.25);
    b1_2->Draw("same");
    
    
    TLine *l2 = new TLine(8.8208,1,8.8208,1400);
    l2->Draw("same");
    TBox *b2 = new TBox(8.8208-0.7e-3,1,8.8208+0.7e-3,1400);
    b2->SetFillColorAlpha(kBlack,0.25);
    b2->Draw("same");
    
    dummy = 0.1;
    for(int i=0;i<gExErrors1sig->GetN();i++)
    {
        
        if(abs(gExErrors1sig->GetPointX(i)-8.8208)<dummy)
        {
//             std::cout << abs/*(*/gExErrors1sig->GetPointX(i)-8.8208) << "\t" << dummy << std::endl;;
            
//             std::cout << "found a new smaller value: " << i << std::endl;
            dummy = abs(gExErrors1sig->GetPointX(i)-8.8208);
            FitUncertaintyAtPoint = gExErrors1sig->GetErrorYhigh(i);
//             std::cout << "FitUncertaintyAtPoint: " << FitUncertaintyAtPoint << std::endl;
            ClosestPoint = i;
        }
    }
    TBox *b2_2 = new TBox(8.8208-FitUncertaintyAtPoint*1.e-3,1,8.8208+FitUncertaintyAtPoint*1.e-3,1400);
    b2_2->SetFillColorAlpha(kBlue,0.25);
    b2_2->Draw("same");
    
    TLine *l2a = new TLine(8.8271,1,8.8271,1400);
    l2a->Draw("same");
    TBox *b2a = new TBox(8.8271-1.1e-3,1,8.8279+1.1e-3,1400);
    b2a->SetFillColorAlpha(kBlack,0.25);
    b2a->Draw("same");
    
    dummy = 0.1;
    for(int i=0;i<gExErrors1sig->GetN();i++)
    {
        
        if(abs(gExErrors1sig->GetPointX(i)-8.8271)<dummy)
        {
//             std::cout << abs(gExErrors1sig->GetPointX(i)-8.8271) << "\t" << dummy << std::endl;;
            
//             std::cout << "found a new smaller value: " << i << std::endl;
            dummy = abs(gExErrors1sig->GetPointX(i)-8.8271);
            FitUncertaintyAtPoint = gExErrors1sig->GetErrorYhigh(i);
//             std::cout << "FitUncertaintyAtPoint: " << FitUncertaintyAtPoint << std::endl;
            ClosestPoint = i;
        }
    }
    TBox *b2a_2 = new TBox(8.8271-FitUncertaintyAtPoint*1.e-3,1,8.8271+FitUncertaintyAtPoint*1.e-3,1400);
    b2a_2->SetFillColorAlpha(kBlue,0.25);
    b2a_2->Draw("same");
    
    TLine *l6 = new TLine(8.862,1,8.862,1400);
    l6->SetLineColor(3);
    l6->SetLineStyle(2);
    l6->SetLineWidth(4);
    l6->Draw("same");
    
    TLine *l5 = new TLine(8.894,1,8.894,1400);
    l5->SetLineColor(3);
    l5->SetLineStyle(2);
    l5->SetLineWidth(4);
    l5->Draw("same");
    
    TLine *l8 = new TLine(8.9451,1,8.9451,1400);//is this a LUNA resonance energy?
    l8->Draw("same");
    TBox *b8 = new TBox(8.9451-0.8e-3,1,8.9451+0.8e-3,1400);
    b8->SetFillColorAlpha(kBlack,0.25);
    b8->Draw("same");
    
    dummy = 0.1;
    for(int i=0;i<gExErrors1sig->GetN();i++)
    {
        
        if(abs(gExErrors1sig->GetPointX(i)-8.9451)<dummy)
        {
            std::cout << abs(gExErrors1sig->GetPointX(i)-8.9451) << "\t" << dummy << std::endl;;
            
            std::cout << "found a new smaller value: " << i << std::endl;
            dummy = abs(gExErrors1sig->GetPointX(i)-8.9451);
            FitUncertaintyAtPoint = gExErrors1sig->GetErrorYhigh(i);
            std::cout << "FitUncertaintyAtPoint: " << FitUncertaintyAtPoint << std::endl;
            ClosestPoint = i;
        }
    }
    TBox *b8_2 = new TBox(8.9451-FitUncertaintyAtPoint*1.e-3,1,8.9451+FitUncertaintyAtPoint*1.e-3,1400);
    b8_2->SetFillColorAlpha(kBlue,0.25);
    b8_2->Draw("same");
    
    TLine *l81 = new TLine(8.9468,1,8.9468,1400);//is this a LUNA resonance energy?
    l81->Draw("same");
    TBox *b81 = new TBox(8.9468-0.6e-3,1,8.9468+0.6e-3,1400);
    b81->SetFillColorAlpha(kBlack,0.25);
    b81->Draw("same");
    
    dummy = 0.1;
    for(int i=0;i<gExErrors1sig->GetN();i++)
    {
        
        if(abs(gExErrors1sig->GetPointX(i)-8.9468)<dummy)
        {
//             std::cout << abs(gExErrors1sig->GetPointX(i)-8.9468) << "\t" << dummy << std::endl;;
            
//             std::cout << "found a new smaller value: " << i << std::endl;
            dummy = abs(gExErrors1sig->GetPointX(i)-8.9468);
            FitUncertaintyAtPoint = gExErrors1sig->GetErrorYhigh(i);
//             std::cout << "FitUncertaintyAtPoint: " << FitUncertaintyAtPoint << std::endl;
            ClosestPoint = i;
        }
    }
    TBox *b81_2 = new TBox(8.9468-FitUncertaintyAtPoint*1.e-3,1,8.9468+FitUncertaintyAtPoint*1.e-3,1400);
    b81_2->SetFillColorAlpha(kBlue,0.25);
    b81_2->Draw("same");
    
    TLine *l4 = new TLine(8.9639,1,8.9639,1400);
    l4->Draw("same");
    TBox *b4 = new TBox(8.9639-1.1e-3,1,8.9639+1.1e-3,1400);
    b4->SetFillColorAlpha(kBlack,0.25);
    b4->Draw("same");
    
    dummy = 0.1;
    for(int i=0;i<gExErrors1sig->GetN();i++)
    {
        
        if(abs(gExErrors1sig->GetPointX(i)-8.9639)<dummy)
        {
//             std::cout << abs(gExErrors1sig->GetPointX(i)-8.9639) << "\t" << dummy << std::endl;;
            
//             std::cout << "found a new smaller value: " << i << std::endl;
            dummy = abs(gExErrors1sig->GetPointX(i)-8.9639);
            FitUncertaintyAtPoint = gExErrors1sig->GetErrorYhigh(i);
//             std::cout << "FitUncertaintyAtPoint: " << FitUncertaintyAtPoint << std::endl;
            ClosestPoint = i;
        }
    }
    TBox *b4_2 = new TBox(8.9639-FitUncertaintyAtPoint*1.e-3,1,8.9639+FitUncertaintyAtPoint*1.e-3,1400);
    b4_2->SetFillColorAlpha(kBlue,0.25);
    b4_2->Draw("same");
    
    TLine *l7 = new TLine(8.9753,1,8.9753,1400);
    l7->Draw("same");
    TBox *b7 = new TBox(8.9753-0.7e-3,1,8.9753+0.7e-3,1400);
    b7->SetFillColorAlpha(kBlack,0.25);
    b7->Draw("same");
    
    dummy = 0.1;
    for(int i=0;i<gExErrors1sig->GetN();i++)
    {
        
        if(abs(gExErrors1sig->GetPointX(i)-8.9753)<dummy)
        {
//             std::cout << abs(gExErrors1sig->GetPointX(i)-8.9753) << "\t" << dummy << std::endl;;
            
//             std::cout << "found a new smaller value: " << i << std::endl;
            dummy = abs(gExErrors1sig->GetPointX(i)-8.9753);
            FitUncertaintyAtPoint = gExErrors1sig->GetErrorYhigh(i);
//             std::cout << "FitUncertaintyAtPoint: " << FitUncertaintyAtPoint << std::endl;
            ClosestPoint = i;
        }
    }
    TBox *b7_2 = new TBox(8.9753-FitUncertaintyAtPoint*1.e-3,1,8.9753+FitUncertaintyAtPoint*1.e-3,1400);
    b7_2->SetFillColorAlpha(kBlue,0.25);
    b7_2->Draw("same");
    
    TLine *l71 = new TLine(8.976,1,8.976,1400);
    l71->Draw("same");
    TBox *b71 = new TBox(8.976-1.e-3,1,8.976+1.0e-3,1400);
    b71->SetFillColorAlpha(kBlack,0.25);
    b71->Draw("same");
    
    dummy = 0.1;
    for(int i=0;i<gExErrors1sig->GetN();i++)
    {
        
        if(abs(gExErrors1sig->GetPointX(i)-8.976)<dummy)
        {
//             std::cout << abs(gExErrors1sig->GetPointX(i)-8.976) << "\t" << dummy << std::endl;;
            
//             std::cout << "found a new smaller value: " << i << std::endl;
            dummy = abs(gExErrors1sig->GetPointX(i)-8.976);
            FitUncertaintyAtPoint = gExErrors1sig->GetErrorYhigh(i);
//             std::cout << "FitUncertaintyAtPoint: " << FitUncertaintyAtPoint << std::endl;
            ClosestPoint = i;
        }
    }
    TBox *b71_2 = new TBox(8.976-FitUncertaintyAtPoint*1.e-3,1,8.976+FitUncertaintyAtPoint*1.e-3,1400);
    b71_2->SetFillColorAlpha(kBlue,0.25);
    b71_2->Draw("same");
    
    TLine *l6a = new TLine(9.000,1,9.000,1400);
    l6a->SetLineColor(3);
    l6a->SetLineStyle(2);
    l6a->SetLineWidth(4);
    l6a->Draw("same");
    
    TLine *l91 = new TLine(9.0395,1,9.0395,1400);//is this a LUNA resonance energy?
    l91->Draw("same");
    TBox *b91 = new TBox(9.0395-0.8e-3,1,9.0395+0.8e-3,1400);
    b91->SetFillColorAlpha(kBlack,0.25);
    b91->Draw("same");
    
    dummy = 0.1;
    for(int i=0;i<gExErrors1sig->GetN();i++)
    {
        
        if(abs(gExErrors1sig->GetPointX(i)-9.0395)<dummy)
        {
//             std::cout << abs(gExErrors1sig->GetPointX(i)-9.0395) << "\t" << dummy << std::endl;;
            
//             std::cout << "found a new smaller value: " << i << std::endl;
            dummy = abs(gExErrors1sig->GetPointX(i)-9.0395);
            FitUncertaintyAtPoint = gExErrors1sig->GetErrorYhigh(i);
//             std::cout << "FitUncertaintyAtPoint: " << FitUncertaintyAtPoint << std::endl;
            ClosestPoint = i;
        }
    }
    TBox *b91_2 = new TBox(9.0395-FitUncertaintyAtPoint*1.e-3,1,9.0395+FitUncertaintyAtPoint*1.e-3,1400);
    b91_2->SetFillColorAlpha(kBlue,0.25);
    b91_2->Draw("same");
    
    TLine *l9 = new TLine(9.0426,1,9.0426,1400);//is this a LUNA resonance energy?
    l9->Draw("same");
    TBox *b9 = new TBox(9.0426-0.8e-3,1,9.0426+0.8e-3,1400);
    b9->SetFillColorAlpha(kBlack,0.25);
    b9->Draw("same");
    
    dummy = 0.1;
    for(int i=0;i<gExErrors1sig->GetN();i++)
    {
        
        if(abs(gExErrors1sig->GetPointX(i)-9.0426)<dummy)
        {
//             std::cout << abs(gExErrors1sig->GetPointX(i)-9.0426) << "\t" << dummy << std::endl;;
            
//             std::cout << "found a new smaller value: " << i << std::endl;
            dummy = abs(gExErrors1sig->GetPointX(i)-9.0426);
            FitUncertaintyAtPoint = gExErrors1sig->GetErrorYhigh(i);
//             std::cout << "FitUncertaintyAtPoint: " << FitUncertaintyAtPoint << std::endl;
            ClosestPoint = i;
        }
    }
    TBox *b9_2 = new TBox(9.0426-FitUncertaintyAtPoint*1.e-3,1,9.0426+FitUncertaintyAtPoint*1.e-3,1400);
    b9_2->SetFillColorAlpha(kBlue,0.25);
    b9_2->Draw("same");
    
    TLine *l3 = new TLine(9.072,1,9.072,1400);
    l3->Draw("same");
    TBox *b3 = new TBox(9.072-3.e-3,1,9.072+3.e-3,1400);
    b3->SetFillColorAlpha(kBlack,0.25);
    b3->Draw("same");
    
    dummy = 0.1;
    for(int i=0;i<gExErrors1sig->GetN();i++)
    {
        
        if(abs(gExErrors1sig->GetPointX(i)-9.072)<dummy)
        {
//             std::cout << abs(gExErrors1sig->GetPointX(i)-9.072) << "\t" << dummy << std::endl;;
            
//             std::cout << "found a new smaller value: " << i << std::endl;
            dummy = abs(gExErrors1sig->GetPointX(i)-9.072);
            FitUncertaintyAtPoint = gExErrors1sig->GetErrorYhigh(i);
//             std::cout << "FitUncertaintyAtPoint: " << FitUncertaintyAtPoint << std::endl;
            ClosestPoint = i;
        }
    }
    TBox *b3_2 = new TBox(9.072-FitUncertaintyAtPoint*1.e-3,1,9.072+FitUncertaintyAtPoint*1.e-3,1400);
    b3_2->SetFillColorAlpha(kBlue,0.25);
    b3_2->Draw("same");
    
//     TLine *l10 = new TLine(9.072,1,9.072,1400);
//     l10->Draw("same");
//     TBox *b10 = new TBox(9.072-3.e-3,1,9.072+3.e-3,1400);
//     b10->SetFillColorAlpha(kBlack,0.25);
//     b10->Draw("same");
    
//     TLine *l11 = new TLine(9.072,1,9.072,1400);
//     l11->Draw("same");
//     TBox *b11 = new TBox(9.072-3.e-3,1,9.072+3.e-3,1400);
//     b11->SetFillColorAlpha(kBlack,0.25);
//     b11->Draw("same");
    
    TLine *l12 = new TLine(9.1015,1,9.1015,1400);
    l12->Draw("same");
    TBox *b12 = new TBox(9.1015-0.7e-3,1,9.1015+0.7e-3,1400);
    b12->SetFillColorAlpha(kBlack,0.25);
    b12->Draw("same");
    
    dummy = 0.1;
    for(int i=0;i<gExErrors1sig->GetN();i++)
    {
        
        if(abs(gExErrors1sig->GetPointX(i)-9.1015)<dummy)
        {
//             std::cout << abs(gExErrors1sig->GetPointX(i)-9.1015) << "\t" << dummy << std::endl;;
            
//             std::cout << "found a new smaller value: " << i << std::endl;
            dummy = abs(gExErrors1sig->GetPointX(i)-9.1015);
            FitUncertaintyAtPoint = gExErrors1sig->GetErrorYhigh(i);
//             std::cout << "FitUncertaintyAtPoint: " << FitUncertaintyAtPoint << std::endl;
            ClosestPoint = i;
        }
    }
    TBox *b12_2 = new TBox(9.1015-FitUncertaintyAtPoint*1.e-3,1,9.1015+FitUncertaintyAtPoint*1.e-3,1400);
    b12_2->SetFillColorAlpha(kBlue,0.25);
    b12_2->Draw("same");
    
    TLine *l13 = new TLine(9.113,1,9.113,1400);
    l13->Draw("same");
    TBox *b13 = new TBox(9.113-3.e-3,1,9.113+3.e-3,1400);
    b13->SetFillColorAlpha(kBlack,0.25);
    b13->Draw("same");
    
    dummy = 0.1;
    for(int i=0;i<gExErrors1sig->GetN();i++)
    {
        
        if(abs(gExErrors1sig->GetPointX(i)-9.113)<dummy)
        {
//             std::cout << abs(gExErrors1sig->GetPointX(i)-9.113) << "\t" << dummy << std::endl;;
            
//             std::cout << "found a new smaller value: " << i << std::endl;
            dummy = abs(gExErrors1sig->GetPointX(i)-9.113);
            FitUncertaintyAtPoint = gExErrors1sig->GetErrorYhigh(i);
//             std::cout << "FitUncertaintyAtPoint: " << FitUncertaintyAtPoint << std::endl;
            ClosestPoint = i;
        }
    }
    TBox *b13_2 = new TBox(9.113-FitUncertaintyAtPoint*1.e-3,1,9.113+FitUncertaintyAtPoint*1.e-3,1400);
    b13_2->SetFillColorAlpha(kBlue,0.25);
    b13_2->Draw("same");
    
    /*
    
    TLine *l2 = new TLine(8.8208,1,8.8208,1400);
    l2->Draw("same");
    TBox *b2 = new TBox(8.8208-0.7e-3,1,8.8208+0.7e-3,1400);
    b2->SetFillColorAlpha(kBlack,0.25);
    b2->Draw("same");
    
    TLine *l2a = new TLine(8.8271,1,8.8271,1400);
    l2a->Draw("same");
    TBox *b2a = new TBox(8.8271-1.1e-3,1,8.8279+1.1e-3,1400);
    b2a->SetFillColorAlpha(kBlack,0.25);
    b2a->Draw("same");
    
    TLine *l6 = new TLine(8.862,1,8.862,1400);
    l6->SetLineColor(3);
    l6->SetLineStyle(2);
    l6->SetLineWidth(4);
    l6->Draw("same");
    
    TLine *l5 = new TLine(8.894,1,8.894,1400);
    l5->SetLineColor(3);
    l5->SetLineStyle(2);
    l5->SetLineWidth(4);
    l5->Draw("same");
    
    TLine *l8 = new TLine(8.9451,1,8.9451,1400);//is this a LUNA resonance energy?
    l8->Draw("same");
    TBox *b8 = new TBox(8.9451-0.8e-3,1,8.9451+0.8e-3,1400);
    b8->SetFillColorAlpha(kBlack,0.25);
    b8->Draw("same");
    
    TLine *l81 = new TLine(8.9468,1,8.9468,1400);//is this a LUNA resonance energy?
    l81->Draw("same");
    TBox *b81 = new TBox(8.9468-0.6e-3,1,8.9468+0.6e-3,1400);
    b81->SetFillColorAlpha(kBlack,0.25);
    b81->Draw("same");
    
    TLine *l4 = new TLine(8.9639,1,8.9639,1400);
    l4->Draw("same");
    TBox *b4 = new TBox(8.9639-1.1e-3,1,8.9639+1.1e-3,1400);
    b4->SetFillColorAlpha(kBlack,0.25);
    b4->Draw("same");
    
    TLine *l7 = new TLine(8.9753,1,8.9753,1400);
    l7->Draw("same");
    TBox *b7 = new TBox(8.9753-0.7e-3,1,8.9753+0.7e-3,1400);
    b7->SetFillColorAlpha(kBlack,0.25);
    b7->Draw("same");
    
    TLine *l71 = new TLine(8.976,1,8.976,1400);
    l71->Draw("same");
    TBox *b71 = new TBox(8.976-1.e-3,1,8.976+1.0e-3,1400);
    b71->SetFillColorAlpha(kBlack,0.25);
    b71->Draw("same");
    
    TLine *l6a = new TLine(9.000,1,9.000,1400);
    l6a->SetLineColor(3);
    l6a->SetLineStyle(2);
    l6a->SetLineWidth(4);
    l6a->Draw("same");
    
    TLine *l91 = new TLine(9.0395,1,9.0395,1400);//is this a LUNA resonance energy?
    l91->Draw("same");
    TBox *b91 = new TBox(9.0395-0.8e-3,1,9.0395+0.8e-3,1400);
    b91->SetFillColorAlpha(kBlack,0.25);
    b91->Draw("same");
    
    TLine *l9 = new TLine(9.0426,1,9.0426,1400);//is this a LUNA resonance energy?
    l9->Draw("same");
    TBox *b9 = new TBox(9.0426-0.8e-3,1,9.0426+0.8e-3,1400);
    b9->SetFillColorAlpha(kBlack,0.25);
    b9->Draw("same");
    
    TLine *l3 = new TLine(9.072,1,9.072,1400);
    l3->Draw("same");
    TBox *b3 = new TBox(9.072-3.e-3,1,9.072+3.e-3,1400);
    b3->SetFillColorAlpha(kBlack,0.25);
    b3->Draw("same");
    
//     TLine *l10 = new TLine(9.072,1,9.072,1400);
//     l10->Draw("same");
//     TBox *b10 = new TBox(9.072-3.e-3,1,9.072+3.e-3,1400);
//     b10->SetFillColorAlpha(kBlack,0.25);
//     b10->Draw("same");
    
//     TLine *l11 = new TLine(9.072,1,9.072,1400);
//     l11->Draw("same");
//     TBox *b11 = new TBox(9.072-3.e-3,1,9.072+3.e-3,1400);
//     b11->SetFillColorAlpha(kBlack,0.25);
//     b11->Draw("same");
    
    TLine *l12 = new TLine(9.1015,1,9.1015,1400);
    l12->Draw("same");
    TBox *b12 = new TBox(9.1015-0.7e-3,1,9.1015+0.7e-3,1400);
    b12->SetFillColorAlpha(kBlack,0.25);
    b12->Draw("same");
    
    TLine *l13 = new TLine(9.113,1,9.113,1400);
    l13->Draw("same");
    TBox *b13 = new TBox(9.113-3.e-3,1,9.113+3.e-3,1400);
    b13->SetFillColorAlpha(kBlack,0.25);
    b13->Draw("same");
    */
    c4->Update();
    c4->SaveAs("Na23_Ex_70deg_PA.eps");
    c4->SaveAs("Na23_Ex_70deg_PA.pdf");
    c4->SaveAs("Na23_Ex_70deg_PA.C");
}

void Draw19FExSpectrum()
{
    TFile *fNaF = TFile::Open("run049_an.root");
    TTree *tNaF = (TTree*)fNaF->Get("readout");
    
    gROOT->ProcessLine(".L NaProtonKinematicsCalculationNew.cpp");
    
    TFile *fLiF = TFile::Open("run051_an.root");
    TTree *tLiF = (TTree*)fLiF->Get("readout");
//     tNaF->Draw("pos>>histoNaF(620,0,2500)","pos>0","");
//     TH1F *histoNaF = (TH1F*)gROOT->FindObjectAny("histoNaF");
    
    tNaF->SetAlias("Brho",Form("%e + %e * pos + %e * pos*pos",p0,p1,p2));
    tLiF->SetAlias("Brho",Form("%e + %e * pos + %e * pos*pos",p0,p1,p2));
    tNaF->Draw("BrhoToExLiF(Brho,938.782980,70.)>>hEx23Na(450,8.55,9.00)","pos>0","goff");
    tLiF->Draw("BrhoToExLiF(Brho,938.782980,70.)>>hEx19F(450,8.55,9.00)","(pos>0)","goff");
    
    
    TH1F *hEx23Na = (TH1F*)gROOT->FindObjectAny("hEx23Na");
    hEx23Na->Draw();
    hEx23Na->SetTitle("");
    hEx23Na->SetStats(0);
    hEx23Na->GetXaxis()->SetTitle("E_{x} [MeV]");
    hEx23Na->GetXaxis()->CenterTitle();
    hEx23Na->GetYaxis()->CenterTitle();
    hEx23Na->GetYaxis()->SetTitle(Form("Counts per %d keV",(int)(hEx23Na->GetBinWidth(1)*1000)));
    
    TH1F *hEx19F = (TH1F*)gROOT->FindObjectAny("hEx19F");
    hEx19F->Draw("same");
    hEx19F->SetLineColor(2);
}

void Draw28SiExSpectrum()
{
    gROOT->ProcessLine(".L NaProtonKinematicsCalculationNew.cpp");
    
    TFile *fSi = TFile::Open("run055_an.root");
    TTree *tSi = (TTree*)fSi->Get("readout");
    
    tSi->SetAlias("Brho","0.300223 + 6.35942e-06 * pos -4.42157e-10 * pos*pos");
    tSi->Draw("BrhoToExSiO2(Brho,938.782980,70.)>>hEx28Si(450,8.8,9.25)","pos>0","goff");
    
    TH1F *hEx28Si = (TH1F*)gROOT->FindObjectAny("hEx28Si");
    hEx28Si->Draw();
    hEx28Si->SetTitle("");
    hEx28Si->SetStats(0);
    hEx28Si->GetXaxis()->SetTitle("E_{x} [MeV]");
    hEx28Si->GetXaxis()->CenterTitle();
    hEx28Si->GetYaxis()->CenterTitle();
    hEx28Si->GetYaxis()->SetTitle(Form("Counts per %d keV",(int)(hEx28Si->GetBinWidth(1)*1000)));
}

double BrhoToExNaF(double Brho, double mass, double theta)
{
  //Calculating proton energy throughout targets
  double Ep = sqrt(pow(TMath::C()/1e6*Brho,2.) + pow(mass,2.)) - mass;
   //cout << "Proton energy from focal plane: " << Ep << endl;
   double CarbonThickness = 0.08849557522;
   double NaFThickness = 0.1953125;
   Ep = ELossC(Ep,-CarbonThickness);
   Ep = ELossNaF(Ep,-0.5*NaFThickness);
   //cout << "Proton energy after target energy-loss correction: " << Ep << endl;

  double p3 = sqrt(Ep * (Ep + 2*mass));

  double T1 = 14.000;
  theta *= TMath::Pi()/180.;
  T1 = ELossNaF(T1,0.5*NaFThickness/cos(theta));

  double p1 = sqrt(T1 * (T1 + 2*mass));

  double p4 = sqrt(p1*p1 - 2*p1*p3*cos(theta) + p3*p3);
//   cout << "p4: " << p4 << endl;

  //double m_28Si = 26060.33946; //MeV/c/c // silicon
  double m_23Na = mNa; //MeV/c/c // sodium
  //double m_28Si = 17696.90035; //MeV/c/c // fluorine

  double T4 = sqrt(p4*p4 + m_23Na*m_23Na) - m_23Na;
//   cout << "T4: " << T4 << endl;

  double Ex = sqrt(pow(T1 - Ep + m_23Na,2.) - p4*p4) - m_23Na;

  return Ex;
}
