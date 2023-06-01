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

double p0 = 3.13067e-01, p1 = 6.76993e-06, p2 = -4.46741e-10;//pos->Brho calibration

double BrhoToExNaF(double Brho, double mass, double theta);

double mNa = 21414.8345021;
double mF = 17696.9005009;
double mSi = 26060.3420706;
double mO = 14899.1686366;
double mLi = 6535.36582157;

TGraphErrors *gExErrors1sig;

//lots of hacky changes to put the uncertainties in Ex from the fits onto the plots.

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


void FitPeak25deg()
{
    TCanvas *c1 = new TCanvas();
    
//     gROOT->ProcessLine(".x run49.C");
    
//     TH1F *histo = (TH1F*)gROOT->FindObjectAny("hPos__49");
    TFile *fNaF = TFile::Open("/home/padsley/data/Q3DNa23pp/rootfiles/run032_an.root");
    TTree *tNaF = (TTree*)fNaF->Get("readout");
    tNaF->Draw("pos>>histoNaF(620,0,2500)","pos>0","");
    TH1F *histoNaF = (TH1F*)gROOT->FindObjectAny("histoNaF");
    histoNaF->Draw("E");
    
    TF1 *fitty = new TF1("fitty",TailedGaussian,1860,1980,6);
    
    std::cout << "Ex = 23Na: 8798.7 keV" << std::endl;
    fitty->SetParameters(1940,27.6,-3,800,200,0);
    fitty->FixParameter(5,0.);
    
    histoNaF->Fit(fitty,"BRME");
    
    c1->SaveAs("figures/FitPeak0.png");
    
    fitty->SetRange(280,385);
    fitty->SetParameters(325,40.,-2.1,600,400,0);
//     fitty->SetParLimits(5,0.,0.);
    std::cout << "Ex = 23Na: 9072 keV" << std::endl;
    histoNaF->Fit(fitty,"BRME");
    
    c1->SaveAs("figures/FitPeak1.png");
    
    fitty->SetRange(850,920);
    fitty->SetParameters(870,25,-3,800,600,0);
//     fitty->SetParLimits(5,0.,0.);
    std::cout << "Ex = 8975.3(7)-keV state in 23Na" << std::endl;
    histoNaF->Fit(fitty,"BRME");
    
    c1->SaveAs("figures/FitPeak2.png");
//     
    fitty->SetRange(1200,1325);
    fitty->SetParameters(1275,23,-3,500,200);
    std::cout << "Ex = 8925.7(28)-keV state in 19F" << std::endl;
    histoNaF->Fit(fitty,"BRME");
    
    c1->SaveAs("figures/FitPeak3.png");
    
    fitty->SetRange(75,141);
    fitty->SetParameters(125,17,-3,800,400);
    std::cout << "Ex = 9113(3)-keV state in 23Na" << std::endl;
    histoNaF->Fit(fitty,"BRME");
    
    c1->SaveAs("figures/FitPeak4.png");
    
//     std::cout << histoNaF->GetNbinsX() << "\t" << histoNaF->GetBinLowEdge(1) << "\t" << histoNaF->GetBinLowEdge(histoNaF->GetNbinsX()+1) << std::endl;
//     TH1F *histoLiF = new TH1F("histoLiF","",625,0,2500);
//     TFile *fLiF = TFile::Open("run051_an.root");
//     TTree *tLiF = (TTree*)fLiF->Get("readout");
//     tLiF->Draw("pos>>histoLiF(625,0,2500)","pos>0","");
//     TH1F *histoLiF = (TH1F*)gROOT->FindObjectAny("histoLiF");
   
    FitRigidities();
}

void FitRigidities()
{
    TGraphErrors *ge = new TGraphErrors();
    
    ge->SetPoint(ge->GetN(),1.29109e+02,0.313955);//9113(3)-keV state in 23Na
    ge->SetPointError(ge->GetN()-1,6.37165e-01,0.000103375);
    
    ge->SetPoint(ge->GetN(),1.27755e+03,0.321173);//8864(4)-keV state in 19F
    ge->SetPointError(ge->GetN()-1,8.51401e-01,0.000134911);
    
    ge->SetPoint(ge->GetN(),8.74915e+02,0.318642);//8975.3(7)-keV state in 23Na
    ge->SetPointError(ge->GetN()-1,7.43117e-01,2.37459e-05);
    
//     ge->SetPoint(ge->GetN(),5.05310e+02,0.315358);//9072(3)-keV state in 23Na - wrong state!
//     ge->SetPointError(ge->GetN()-1,2.58634e-01,0.000102881);
    
    ge->SetPoint(ge->GetN(),1.94614e+03,0.324548);//8798.7-keV state in 23Na
    ge->SetPointError(ge->GetN()-1,1.09353e+00,2.6618e-05);
    
//     ge->SetPoint(ge->GetN(),2284.59,0.31178090);//19F 8650-keV state
//     ge->SetPoint(ge->GetN()-1,1260,0.30927948);//28Si 8905-keV state
    ge->SetMarkerStyle(8);
 
    TCanvas *c2 = new TCanvas();
    
    ge->Draw("AP");
    
    TF1 *fPol2 = new TF1("fPol2","[0] + [1]*x + [2]*x*x",0,2560);
    fPol2->SetParameters(0.3,4.e-6,-4.e-10);
    
    ge->Fit(fPol2,"BRME");
    
    int npoints = 2500;
    TGraphErrors *grint = new TGraphErrors(npoints);
    TGraphErrors *grint2 = new TGraphErrors(npoints);
    TGraphErrors *grint3 = new TGraphErrors(npoints);
    
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
    grint->Draw("same3");
    grint2->Draw("same 3");
    
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
    
    c2->SaveAs("BrhoCalibration25deg.eps");
    c2->SaveAs("BrhoCalibration25deg.pdf");
    c2->SaveAs("BrhoCalibration25deg.png");
    c2->SaveAs("BrhoCalibration25deg.C");
    
    TCanvas *c3 = new TCanvas();
    
    TGraphErrors *gExErrors = new TGraphErrors();
    gExErrors1sig = new TGraphErrors();
    TGraphErrors *gExErrors2sig = new TGraphErrors();
    TGraphErrors *gExErrors3sig = new TGraphErrors();
    
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
        
        gExErrors1sig->SetPoint(i,BrhoToExNaF(y,938.782980,25.),0);
        gExErrors1sig->SetPointError(i,0, 1000.*(BrhoToExNaF(y,938.782980,25.) - BrhoToExNaF(y + sig_y,938.782980,25.))); 
        
        
        x = grint2->GetPointX(i);
        y = grint2->GetPointY(i);
        sig_y = grint2->GetErrorYhigh(i);
        
        gExErrors2sig->SetPoint(i,BrhoToExNaF(y,938.782980,25.),0);
        gExErrors2sig->SetPointError(i,0, 1000.*(BrhoToExNaF(y,938.782980,25.) - BrhoToExNaF(y + sig_y,938.782980,25.))); 
        
        x = grint3->GetPointX(i);
        y = grint3->GetPointY(i);
        sig_y = grint3->GetErrorYhigh(i);
        
        gExErrors3sig->SetPoint(i,BrhoToExNaF(y,938.782980,25.),0);
        gExErrors3sig->SetPointError(i,0, 1000.*(BrhoToExNaF(y,938.782980,25.) - BrhoToExNaF(y + sig_y,938.782980,25.))); 
        
//         std::cout << "i = " << i << " y = " << y << "\t sig_y = " << sig_y << std::endl;
        gEx3sig_ActualValues->SetPoint(i,BrhoToExNaF(y,938.782980,25.),1000.*(BrhoToExNaF(y,938.782980,25.) - BrhoToExNaF(y + sig_y,938.782980,25.)));
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
        
        c3->SaveAs("ExUncertaintyBand25deg.eps");
        c3->SaveAs("ExUncertaintyBand25deg.pdf");
        c3->SaveAs("ExUncertaintyBand25deg.png");
        c3->SaveAs("ExUncertaintyBand25deg.C");
        Draw23NaExSpectrum();
}

void Draw23NaExSpectrum()
{
    TCanvas *c4 = new TCanvas();
    c4->SetLeftMargin(0.136);
    
//     gROOT->ProcessLine(".L NaProtonKinematicsCalculationNew.cpp");
    
    TFile *fNaF = TFile::Open("/home/padsley/data/Q3DNa23pp/rootfiles/run032_an.root");
    TTree *tNaF = (TTree*)fNaF->Get("readout");
//     tNaF->Draw("pos>>histoNaF(620,0,2500)","pos>0","");
//     TH1F *histoNaF = (TH1F*)gROOT->FindObjectAny("histoNaF");
    
    tNaF->SetAlias("Brho",Form("%e + %e * pos + %e * pos*pos",p0,p1,p2));
    
    std::cout << "Brho = " << Form("%e + %e * pos + %e * pos*pos",p0,p1,p2) << std::endl;
    tNaF->Draw("BrhoToExNaF(Brho,938.782980,25.)>>hEx23Na(450,8.7,9.15)","pos>0 && pos<2500","");
    
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
    hEx23Na->GetYaxis()->SetRangeUser(1,5800);
    
    TFile *fC = TFile::Open("/home/padsley/data/Q3DNa23pp/rootfiles/run034_an.root");
    TTree *tC = (TTree*)fC->Get("readout");
    tC->SetAlias("Brho",Form("%e + %e * pos + %e * pos*pos",p0,p1,p2));
    tC->Draw("BrhoToExNaF(Brho,938.782980,25.)>>hEx19F(450,8.7,9.15)","pos>0 && pos<2500","goff");
    
    TH1F *hEx19F = (TH1F*)gROOT->FindObjectAny("hEx19F");
    hEx19F->SetLineColor(2);
    hEx19F->SetFillColorAlpha(kRed,0.35);
    hEx19F->Scale(1.0);
    hEx19F->Draw("same hist");
    
    double Sp = 8.79410;
    
    double axisLength = 350; //in keV
    
    TGaxis *EcmAxis = new TGaxis(Sp,5800,Sp+1.e-3*axisLength,5800,0,axisLength,510,"L-");
    EcmAxis->SetTextFont(42);
    EcmAxis->SetLabelFont(42);
    EcmAxis->SetLabelSize(0.05);
    EcmAxis->SetTitle("E_{cm} [keV]");
    EcmAxis->SetTitleSize(0.05);
    EcmAxis->SetTitleOffset(0.9);
    EcmAxis->CenterTitle();
    EcmAxis->Draw();
    
    TLine *l1 = new TLine(8.7987,1,8.7987,5000);
    l1->Draw("same");
    TBox *b1 = new TBox(8.7987-0.8e-3,1,8.7987+0.8e-3,5000);
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
    TBox *b1_2 = new TBox(8.7987-FitUncertaintyAtPoint*1.e-3,1,8.7987+FitUncertaintyAtPoint*1.e-3,5000);
    b1_2->SetFillColorAlpha(kBlue,0.25);
    b1_2->Draw("same");
    
    
    TLine *l2 = new TLine(8.8208,1,8.8208,5000);
    l2->Draw("same");
    TBox *b2 = new TBox(8.8208-0.7e-3,1,8.8208+0.7e-3,5000);
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
    TBox *b2_2 = new TBox(8.8208-FitUncertaintyAtPoint*1.e-3,1,8.8208+FitUncertaintyAtPoint*1.e-3,5000);
    b2_2->SetFillColorAlpha(kBlue,0.25);
    b2_2->Draw("same");
    
    TLine *l2a = new TLine(8.8271,1,8.8271,5000);
    l2a->Draw("same");
    TBox *b2a = new TBox(8.8271-1.1e-3,1,8.8279+1.1e-3,5000);
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
    TBox *b2a_2 = new TBox(8.8271-FitUncertaintyAtPoint*1.e-3,1,8.8271+FitUncertaintyAtPoint*1.e-3,5000);
    b2a_2->SetFillColorAlpha(kBlue,0.25);
    b2a_2->Draw("same");
    
    TLine *l6 = new TLine(8.862,1,8.862,5000);
    l6->SetLineColor(3);
    l6->SetLineStyle(2);
    l6->SetLineWidth(4);
    l6->Draw("same");
    
    TLine *l5 = new TLine(8.894,1,8.894,5000);
    l5->SetLineColor(3);
    l5->SetLineStyle(2);
    l5->SetLineWidth(4);
    l5->Draw("same");
    
    TLine *l8 = new TLine(8.9451,1,8.9451,5000);//is this a LUNA resonance energy?
    l8->Draw("same");
    TBox *b8 = new TBox(8.9451-0.8e-3,1,8.9451+0.8e-3,5000);
    b8->SetFillColorAlpha(kBlack,0.25);
    b8->Draw("same");
    
    dummy = 0.1;
    for(int i=0;i<gExErrors1sig->GetN();i++)
    {
        
        if(abs(gExErrors1sig->GetPointX(i)-8.9451)<dummy)
        {
//             std::cout << abs(gExErrors1sig->GetPointX(i)-8.9451) << "\t" << dummy << std::endl;;
            
//             std::cout << "found a new smaller value: " << i << std::endl;
            dummy = abs(gExErrors1sig->GetPointX(i)-8.9451);
            FitUncertaintyAtPoint = gExErrors1sig->GetErrorYhigh(i);
//             std::cout << "FitUncertaintyAtPoint: " << FitUncertaintyAtPoint << std::endl;
            ClosestPoint = i;
        }
    }
    TBox *b8_2 = new TBox(8.9451-FitUncertaintyAtPoint*1.e-3,1,8.9451+FitUncertaintyAtPoint*1.e-3,5000);
    b8_2->SetFillColorAlpha(kBlue,0.25);
    b8_2->Draw("same");
    
    TLine *l81 = new TLine(8.9468,1,8.9468,5000);//is this a LUNA resonance energy?
    l81->Draw("same");
    TBox *b81 = new TBox(8.9468-0.6e-3,1,8.9468+0.6e-3,5000);
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
    TBox *b81_2 = new TBox(8.9468-FitUncertaintyAtPoint*1.e-3,1,8.9468+FitUncertaintyAtPoint*1.e-3,5000);
    b81_2->SetFillColorAlpha(kBlue,0.25);
    b81_2->Draw("same");
    
    TLine *l4 = new TLine(8.9639,1,8.9639,5000);
    l4->Draw("same");
    TBox *b4 = new TBox(8.9639-1.1e-3,1,8.9639+1.1e-3,5000);
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
    TBox *b4_2 = new TBox(8.9639-FitUncertaintyAtPoint*1.e-3,1,8.9639+FitUncertaintyAtPoint*1.e-3,5000);
    b4_2->SetFillColorAlpha(kBlue,0.25);
    b4_2->Draw("same");
    
    TLine *l7 = new TLine(8.9753,1,8.9753,5000);
    l7->Draw("same");
    TBox *b7 = new TBox(8.9753-0.7e-3,1,8.9753+0.7e-3,5000);
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
    TBox *b7_2 = new TBox(8.9753-FitUncertaintyAtPoint*1.e-3,1,8.9753+FitUncertaintyAtPoint*1.e-3,5000);
    b7_2->SetFillColorAlpha(kBlue,0.25);
    b7_2->Draw("same");
    
    TLine *l71 = new TLine(8.976,1,8.976,5000);
    l71->Draw("same");
    TBox *b71 = new TBox(8.976-1.e-3,1,8.976+1.0e-3,5000);
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
    TBox *b71_2 = new TBox(8.976-FitUncertaintyAtPoint*1.e-3,1,8.976+FitUncertaintyAtPoint*1.e-3,5000);
    b71_2->SetFillColorAlpha(kBlue,0.25);
    b71_2->Draw("same");
    
    TLine *l6a = new TLine(9.000,1,9.000,5000);
    l6a->SetLineColor(3);
    l6a->SetLineStyle(2);
    l6a->SetLineWidth(4);
    l6a->Draw("same");
    
    TLine *l91 = new TLine(9.0395,1,9.0395,5000);//is this a LUNA resonance energy?
    l91->Draw("same");
    TBox *b91 = new TBox(9.0395-0.8e-3,1,9.0395+0.8e-3,5000);
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
    TBox *b91_2 = new TBox(9.0395-FitUncertaintyAtPoint*1.e-3,1,9.0395+FitUncertaintyAtPoint*1.e-3,5000);
    b91_2->SetFillColorAlpha(kBlue,0.25);
    b91_2->Draw("same");
    
    TLine *l9 = new TLine(9.0426,1,9.0426,5000);//is this a LUNA resonance energy?
    l9->Draw("same");
    TBox *b9 = new TBox(9.0426-0.8e-3,1,9.0426+0.8e-3,5000);
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
    TBox *b9_2 = new TBox(9.0426-FitUncertaintyAtPoint*1.e-3,1,9.0426+FitUncertaintyAtPoint*1.e-3,5000);
    b9_2->SetFillColorAlpha(kBlue,0.25);
    b9_2->Draw("same");
    
    TLine *l3 = new TLine(9.072,1,9.072,5000);
    l3->Draw("same");
    TBox *b3 = new TBox(9.072-3.e-3,1,9.072+3.e-3,5000);
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
    TBox *b3_2 = new TBox(9.072-FitUncertaintyAtPoint*1.e-3,1,9.072+FitUncertaintyAtPoint*1.e-3,5000);
    b3_2->SetFillColorAlpha(kBlue,0.25);
    b3_2->Draw("same");
    
//     TLine *l10 = new TLine(9.072,1,9.072,5000);
//     l10->Draw("same");
//     TBox *b10 = new TBox(9.072-3.e-3,1,9.072+3.e-3,5000);
//     b10->SetFillColorAlpha(kBlack,0.25);
//     b10->Draw("same");
    
//     TLine *l11 = new TLine(9.072,1,9.072,5000);
//     l11->Draw("same");
//     TBox *b11 = new TBox(9.072-3.e-3,1,9.072+3.e-3,5000);
//     b11->SetFillColorAlpha(kBlack,0.25);
//     b11->Draw("same");
    
    TLine *l12 = new TLine(9.1015,1,9.1015,5000);
    l12->Draw("same");
    TBox *b12 = new TBox(9.1015-0.7e-3,1,9.1015+0.7e-3,5000);
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
    TBox *b12_2 = new TBox(9.1015-FitUncertaintyAtPoint*1.e-3,1,9.1015+FitUncertaintyAtPoint*1.e-3,5000);
    b12_2->SetFillColorAlpha(kBlue,0.25);
    b12_2->Draw("same");
    
    TLine *l13 = new TLine(9.113,1,9.113,5000);
    l13->Draw("same");
    TBox *b13 = new TBox(9.113-3.e-3,1,9.113+3.e-3,5000);
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
    TBox *b13_2 = new TBox(9.113-FitUncertaintyAtPoint*1.e-3,1,9.113+FitUncertaintyAtPoint*1.e-3,5000);
    b13_2->SetFillColorAlpha(kBlue,0.25);
    b13_2->Draw("same");
    
    c4->Update();
    c4->SaveAs("Na23_Ex_25deg_PA.eps");
    c4->SaveAs("Na23_Ex_25deg_PA.pdf");
    c4->SaveAs("Na23_Ex_25deg_PA.C");
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

// double ExToBrhoNaF_Na(double fBeamEnergy, double Ex, double thetaLab)
// {
//   //print initial beam energy, calculate beam energy after target and print out
//   //cout << "Initial Beam Energy: " << fBeamEnergy << endl;
//     thetaLab *= TMath::Pi()/180.; //converting angle into radians
//     double CarbonThickness = 0.08849557522;
//     double NaFThickness = 0.1953125;
//     fBeamEnergy = ELossNaF(fBeamEnergy,0.5*NaFThickness/cos(thetaLab)); //Proton energy after going through SiO2 target
//   //cout << "Beam energy after target: " << fBeamEnergy << endl;
// 
//   //defining masses
//   double m1 = 938.782980; //proton
//   double m2 = mNa;        // Sodium
//   double m3 = m1;         //proton
//   double m4 = m2 + Ex;    //SiO2 + excitation energy
//   double s = m1*m1 + m2*m2 + 2*m2*(fBeamEnergy + m1);
// 
//   TLorentzVector fTotalEnergyImpulsionCM = TLorentzVector(0,0,0,sqrt(s));
// 
//   //Calculating energies in center of mass reference frame
//   double ECM_1 = (s + m1*m1 - m2*m2)/(2*sqrt(s));
//   double ECM_2 = (s + m2*m2 - m1*m1)/(2*sqrt(s));
//   double ECM_3 = (s + m3*m3 - m4*m4)/(2*sqrt(s));
//   double ECM_4 = (s + m4*m4 - m3*m3)/(2*sqrt(s));
// 
//   //Calculating momentum in center of mass frame with formula p^2 = E^2 - m^2
//   double pCM_1 = sqrt(ECM_1*ECM_1 - m1*m1);
//   double pCM_2 = sqrt(ECM_2*ECM_2 - m2*m2);
//   double pCM_3 = sqrt(ECM_3*ECM_3 - m3*m3);
//   double pCM_4 = sqrt(ECM_4*ECM_4 - m4*m4);
// 
//   TVector3 fImpulsionLab_1 = TVector3(0,0,sqrt(fBeamEnergy*fBeamEnergy + 2*fBeamEnergy*m1));
//   TVector3 fImpulsionLab_2 = TVector3(0,0,0);
// 
//   TLorentzVector fEnergyImpulsionLab_1 = TLorentzVector(fImpulsionLab_1,m1+fBeamEnergy);
//   TLorentzVector fEnergyImpulsionLab_2 = TLorentzVector(fImpulsionLab_2,m2);
// 
//   TLorentzVector fTotalEnergyImpulsionLab = fEnergyImpulsionLab_1 + fEnergyImpulsionLab_2;
// 
//   double BetaCM = fTotalEnergyImpulsionLab.Beta();
// 
//   TLorentzVector fEnergyImpulsionCM_1 = fEnergyImpulsionLab_1;
//   fEnergyImpulsionCM_1.Boost(0,0,-BetaCM);
// 
//   TLorentzVector fEnergyImpulsionCM_2 = fEnergyImpulsionLab_2;
//   fEnergyImpulsionCM_2.Boost(0,0,-BetaCM);
// 
//   double thetaCM = 0;
//   double rapidity = TMath::ATanH(BetaCM);
//   if(pCM_3 - pCM_4 > 1e-6)cout << "I don't understand this: " << pCM_3 << "\t" << pCM_4 << endl;
//   double p3 = (sqrt(m3*m3 + pCM_3*pCM_3) * cos(thetaLab) * sinh(rapidity) + cosh(rapidity) * sqrt(pCM_3*pCM_3 - m3*m3*sin(thetaLab)*sin(thetaLab)*sinh(rapidity)*sinh(rapidity))) / (1.0 + sin(thetaLab)*sin(thetaLab) * sinh(rapidity)*sinh(rapidity));
//   thetaCM = asin(p3 * sin(thetaLab) / pCM_3);
// //   cout << "thetaCM: " << thetaCM*180./TMath::Pi() << endl;
// 
//   TLorentzVector fEnergyImpulsionCM_3  = TLorentzVector(pCM_3*sin(thetaCM),0,pCM_3*cos(thetaCM),ECM_3);
//   TLorentzVector fEnergyImpulsionCM_4  = fTotalEnergyImpulsionCM - fEnergyImpulsionCM_3;
// 
//   TLorentzVector fEnergyImpulsionLab_3 = fEnergyImpulsionCM_3;
//   fEnergyImpulsionLab_3.Boost(0,0,BetaCM);
//   TLorentzVector fEnergyImpulsionLab_4 = fEnergyImpulsionCM_4;
//   fEnergyImpulsionLab_4.Boost(0,0,BetaCM);
// 
//   // Angle in the lab frame
//   double ThetaLab3 = fEnergyImpulsionLab_3.Angle(fEnergyImpulsionLab_1.Vect());
//   if (ThetaLab3 < 0) ThetaLab3 += M_PI;
// //   cout << "ThetaLab3: " << ThetaLab3*180./TMath::Pi() << endl;
// 
//   double ThetaLab4 = fEnergyImpulsionLab_4.Angle(fEnergyImpulsionLab_1.Vect());
//   if (fabs(ThetaLab3) < 1e-6) ThetaLab3 = 0;
//   ThetaLab4 = fabs(ThetaLab4);
//   if (fabs(ThetaLab4) < 1e-6) ThetaLab4 = 0;
// 
//   // Kinetic Energy in the lab frame
//   double KineticEnergyLab3 = fEnergyImpulsionLab_3.E() - m3;
//   double KineticEnergyLab4 = fEnergyImpulsionLab_4.E() - m4;
// //   cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;
// 
//   // test for total energy conversion
//   if (fabs(fTotalEnergyImpulsionLab.E() - (fEnergyImpulsionLab_3.E()+fEnergyImpulsionLab_4.E())) > 1e-6)
//     cout << "Problem for energy conservation" << endl;
// 
//   KineticEnergyLab3 = ELossNaF(KineticEnergyLab3,0.5*NaFThickness);
//   KineticEnergyLab3 = ELossC(KineticEnergyLab3,CarbonThickness);
// //   cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;
// //   cout << "Brho: " << EpToBrho(KineticEnergyLab3,m1) << endl;
// 
//   return EpToBrho(KineticEnergyLab3,m1);
// }
