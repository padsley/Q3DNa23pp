{
    std::cout << "Start" << std::endl;
    
    double Ex = 8975.3e-3, sigmaEx = 0.7e-3, thetaQ3D = 50.;
    
    std::cout << "Load kinematics calculator" << std::endl;
    gROOT->ProcessLine(".L NaProtonKinematicsCalculationNew.cpp");

    double CentroidBrho = ExToBrhoNaF_Na(14.,Ex,thetaQ3D);
    
    std::cout << "Centroid Brho value: " << CentroidBrho << " Tm" << std::endl;
   
    TF1 *fEx = new TF1("fEx","TMath::Gaus(x,[0],[1],0)",Ex - 5*sigmaEx, Ex + 5*sigmaEx);
    fEx->SetParameters(Ex,sigmaEx);
    fEx->SetNpx(1.e6);
    
    TCanvas *c1 = new TCanvas();
    
    fEx->Draw();
    fEx->SetTitle("");
    fEx->GetXaxis()->SetTitle("E_{x} [MeV]");
    fEx->GetXaxis()->CenterTitle();
 
    int nSamples = 10000;
    
    
    TH1F *hEx = new TH1F("hEx","",1000,Ex-10*sigmaEx,Ex+10*sigmaEx);
    hEx->SetStats(0);
    
    double minRange = ExToBrhoNaF_Na(14.,Ex,thetaQ3D) - 0.0005;
    double maxRange = ExToBrhoNaF_Na(14.,Ex,thetaQ3D) + 0.0005;

    TH1F *hBrho = new TH1F("hBrho","",1000,minRange,maxRange);
    hBrho->SetStats(0);
    
    TCanvas *c2 = new TCanvas();
    TCanvas *c3 = new TCanvas();
    
    for(int i=0;i<nSamples;i++)
    {
        double tempEx = fEx->GetRandom(Ex - 5*sigmaEx,Ex + 5*sigmaEx);
        
        hEx->Fill(tempEx);
        
        hBrho->Fill(ExToBrhoNaF_Na(14.,tempEx,thetaQ3D));
        
    }
    
    c2->cd();
    hEx->Draw();
    hEx->GetXaxis()->SetTitle("E_{x} [MeV]");
    hEx->GetXaxis()->CenterTitle();
    
    TF1 *fitEx = new TF1("fitEx","[0] * TMath::Gaus(x,[1],[2],0)",5*sigmaEx,Ex + 5*sigmaEx);
    fitEx->SetParameters(200,Ex,sigmaEx);
    hEx->Fit(fitEx,"BRMLE");
    
    c3->cd();
    hBrho->Draw();
    hBrho->GetXaxis()->SetTitle("B#rho [Tm]");
    hBrho->GetXaxis()->CenterTitle();
    TF1 *fitBrho = new TF1("fitBrho","[0] * TMath::Gaus(x,[1],[2],0)",minRange,maxRange);
    fitBrho->SetParameters(200,CentroidBrho,1.e-6);
    hBrho->Fit(fitBrho,"BRMLE");
    
    std::cout << fitBrho->GetParameter(1) << "\t" << fitBrho->GetParameter(2) << std::endl;
}
