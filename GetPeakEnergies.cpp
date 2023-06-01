#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

double TailedGaussian(double *x, double *pars);
double DoubleTailedGaussian(double *x, double *pars);

void GetPeakEnergies()
{
    //70 degrees
    gROOT->ProcessLine(".x Na23_Ex_70deg_PA.C");
    
    TCanvas *c_70 = (TCanvas*)gROOT->FindObjectAny("c1_n4");
//     c_70->ls();
    TH1F *h70 = (TH1F*)c_70->FindObject("hEx23Na__1");
    h70->ls();
    
    TF1 *f1peak = new TF1("f1peak",&DoubleTailedGaussian,9.095,9.125,8);
//     TF1 *f1peak = new TF1("f1peak",&TailedGaussian,9.065,9.085,6);
//     f1peak->SetParameters(8.975,3.7e-3,-26,5000,100,0);
    
//     f1peak->FixParameter(5,0);
//     h70->Fit(f1peak,"BRLME");
    
    
//     gROOT->ProcessLine(".x ExUncertaintyBand70deg.C");
    
//     TGraphErrors *g1sigma = (TGraphErrors*)((TCanvas*)gROOT->FindObjectAny("c1_n4"));
    
//     ofstream mOutputState1;
//     mOutputState1.open("ExState1_8798keV.dat");
//     mOutputState1 << "70\t" << f1peak->GetParameter(1) << "\t" << f1peak->GetParErrors(1) << 
    
    
}

double TailedGaussian(double *x, double *pars)
{
    //[-1]: I can't count :(
    //[0]: mu
    //[1]: sigma
    //[2]: kappa
    //[3]: amplitude
    //[4]: const bg
    //[5]: linear bg
    
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

double DoubleTailedGaussian(double *x, double *pars)
{
    //[-1]: I can't count :(
    //[0]: mu
    //[1]: sigma
    //[2]: kappa
    //[3]: amplitude
    //[4]: const bg
    //[5]: linear bg
    //[6]: mu too
    //[7]: amplitude too
    
    double result = 0;
    
    if(pars[2] < (x[0] - pars[0])/pars[1])//this is the Gaussian part
    {
        result += pars[3] * exp ( - 0.5 * pow(x[0] - pars[0],2.)/pow(pars[1],2.));
    }
    else //this is the tail part
    {
        result += pars[3] * exp ( 0.5 * pow(pars[2],2.) - pars[2] * (x[0] - pars[0]) / pars[1]);
    }
    if(pars[2] < (x[0] - pars[6])/pars[1])
    {
        result += pars[7] * exp ( - 0.5 * pow(x[0] - pars[6],2.)/pow(pars[1],2.));
    }
    else
    {
        result += pars[7] * exp ( 0.5 * pow(pars[2],2.) - pars[2] * (x[0] - pars[6]) / pars[1]);
    }
    
    result += pars[4] + pars[5] * x[0];//constant+linear background
    
    return result;
}
