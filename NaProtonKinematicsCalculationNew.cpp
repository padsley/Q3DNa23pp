#include <Math/IFunction.h> //library files to be used in this code

//declaring masses
double mNa = 21414.8345021;
double mF = 17696.9005009;
double mSi = 26060.3420706;
double mO = 14899.1686366;
double mLi = 6535.36582157;

double EpToBrho(double Ep, double mass) //function to calculate mag rigidity (Brho) from the proton energy
{
  double Brho = 1e6/TMath::C() * sqrt(Ep * (Ep + 2*mass));
  return Brho;
}

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

//--------------------------------------------------------------------------------------------------------------//

double ExToBrhoSiO2_Si(double fBeamEnergy, double Ex, double thetaLab)
{
  //print initial beam energy, calculate beam energy after target and print out
//   cout << "Initial Beam Energy: " << fBeamEnergy << endl;
  thetaLab *= TMath::Pi()/180.; //converting angle in radians
  double CarbonThickness = 0.08849557522;
  double SiO2Thickness = 0.15094339622;
  fBeamEnergy = ELossSiO2(fBeamEnergy,0.5*SiO2Thickness/cos(thetaLab));
//   cout << "Beam energy after SiO2 target: " << fBeamEnergy << endl;

  //defining masses
  double m1 = 938.782980;  //proton
  double m2 = mSi;         //silicon
  double m3 = m1;          //proton
  double m4 = m2 + Ex;     //silicon oxide + excitation energy
  double s = m1*m1 + m2*m2 + 2*m2*(fBeamEnergy + m1);
  TLorentzVector fTotalEnergyImpulsionCM = TLorentzVector(0,0,0,sqrt(s)); //impulsion = momentum

  //Calculating energies in center of mass reference frame
  double ECM_1 = (s + m1*m1 - m2*m2)/(2*sqrt(s));
  double ECM_2 = (s + m2*m2 - m1*m1)/(2*sqrt(s));
  double ECM_3 = (s + m3*m3 - m4*m4)/(2*sqrt(s));
  double ECM_4 = (s + m4*m4 - m3*m3)/(2*sqrt(s));

  //Calculating momentum in center of mass frame with formula p^2 = E^2 - m^2
  double pCM_1 = sqrt(ECM_1*ECM_1 - m1*m1);
  double pCM_2 = sqrt(ECM_2*ECM_2 - m2*m2);
  double pCM_3 = sqrt(ECM_3*ECM_3 - m3*m3);
  double pCM_4 = sqrt(ECM_4*ECM_4 - m4*m4);

  TVector3 fImpulsionLab_1 = TVector3(0,0,sqrt(fBeamEnergy*fBeamEnergy + 2*fBeamEnergy*m1));
  TVector3 fImpulsionLab_2 = TVector3(0,0,0);

  TLorentzVector fEnergyImpulsionLab_1 = TLorentzVector(fImpulsionLab_1,m1+fBeamEnergy);
  TLorentzVector fEnergyImpulsionLab_2 = TLorentzVector(fImpulsionLab_2,m2);

  TLorentzVector fTotalEnergyImpulsionLab = fEnergyImpulsionLab_1 + fEnergyImpulsionLab_2;

  double BetaCM = fTotalEnergyImpulsionLab.Beta();

  TLorentzVector fEnergyImpulsionCM_1 = fEnergyImpulsionLab_1;
  fEnergyImpulsionCM_1.Boost(0,0,-BetaCM);

  TLorentzVector fEnergyImpulsionCM_2 = fEnergyImpulsionLab_2;
  fEnergyImpulsionCM_2.Boost(0,0,-BetaCM);

  double thetaCM = 0;
  double rapidity = TMath::ATanH(BetaCM);
  if(pCM_3 - pCM_4 > 1e-6)cout << "I don't understand this: " << pCM_3 << "\t" << pCM_4 << endl;
  double p3 = (sqrt(m3*m3 + pCM_3*pCM_3) * cos(thetaLab) * sinh(rapidity) + cosh(rapidity) * sqrt(pCM_3*pCM_3 - m3*m3*sin(thetaLab)*sin(thetaLab)*sinh(rapidity)*sinh(rapidity))) / (1.0 + sin(thetaLab)*sin(thetaLab) * sinh(rapidity)*sinh(rapidity));
  thetaCM = asin(p3 * sin(thetaLab) / pCM_3);
  //cout << "thetaCM: " << TMath::Pi() << endl;
//   cout << "thetaCM: " << thetaCM*180./TMath::Pi() << endl;

  TLorentzVector fEnergyImpulsionCM_3  = TLorentzVector(pCM_3*sin(thetaCM),0,pCM_3*cos(thetaCM),ECM_3);
  TLorentzVector fEnergyImpulsionCM_4  = fTotalEnergyImpulsionCM - fEnergyImpulsionCM_3;

  TLorentzVector fEnergyImpulsionLab_3 = fEnergyImpulsionCM_3;
  fEnergyImpulsionLab_3.Boost(0,0,BetaCM);
  TLorentzVector fEnergyImpulsionLab_4 = fEnergyImpulsionCM_4;
  fEnergyImpulsionLab_4.Boost(0,0,BetaCM);

  // Angle in the lab frame
  double ThetaLab3 = fEnergyImpulsionLab_3.Angle(fEnergyImpulsionLab_1.Vect());
  if (ThetaLab3 < 0) ThetaLab3 += M_PI;
//   cout << "ThetaLab3: " << ThetaLab3*180./TMath::Pi() << endl;

  double ThetaLab4 = fEnergyImpulsionLab_4.Angle(fEnergyImpulsionLab_1.Vect());
  if (fabs(ThetaLab3) < 1e-6) ThetaLab3 = 0;
  ThetaLab4 = fabs(ThetaLab4);
  if (fabs(ThetaLab4) < 1e-6) ThetaLab4 = 0;

  // Kinetic Energy in the lab frame
  double KineticEnergyLab3 = fEnergyImpulsionLab_3.E() - m3;
  double KineticEnergyLab4 = fEnergyImpulsionLab_4.E() - m4;
//   cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;

  // test for total energy conversion
  if (fabs(fTotalEnergyImpulsionLab.E() - (fEnergyImpulsionLab_3.E()+fEnergyImpulsionLab_4.E())) > 1e-6)
    cout << "Problem for energy conservation" << endl;

  KineticEnergyLab3 = ELossSiO2(KineticEnergyLab3,0.5*SiO2Thickness);
  KineticEnergyLab3 = ELossC(KineticEnergyLab3,CarbonThickness);
//   cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;
//   cout << "Brho: " << EpToBrho(KineticEnergyLab3,m1) << endl;

  return EpToBrho(KineticEnergyLab3,m1);
}

double ExToBrhoSiO2_O(double fBeamEnergy, double Ex, double thetaLab)
{
  //print initial beam energy, calculate beam energy after target and print out
//   cout << "Initial Beam Energy: " << fBeamEnergy << endl;
  thetaLab *= TMath::Pi()/180.; //converting angle in radians
  double CarbonThickness = 0.08849557522;
  double SiO2Thickness = 0.15094339622;
  fBeamEnergy = ELossSiO2(fBeamEnergy,0.5*SiO2Thickness/cos(thetaLab));
//   cout << "Beam energy after SiO2 target: " << fBeamEnergy << endl;

  //defining masses
  double m1 = 938.782980;  //proton
  double m2 = mO;          //oxygen
  double m3 = m1;          //proton
  double m4 = m2 + Ex;     //silicon oxide + excitation energy
  double s = m1*m1 + m2*m2 + 2*m2*(fBeamEnergy + m1);
  TLorentzVector fTotalEnergyImpulsionCM = TLorentzVector(0,0,0,sqrt(s)); //impulsion = momentum

  //Calculating energies in center of mass reference frame
  double ECM_1 = (s + m1*m1 - m2*m2)/(2*sqrt(s));
  double ECM_2 = (s + m2*m2 - m1*m1)/(2*sqrt(s));
  double ECM_3 = (s + m3*m3 - m4*m4)/(2*sqrt(s));
  double ECM_4 = (s + m4*m4 - m3*m3)/(2*sqrt(s));

  //Calculating momentum in center of mass frame with formula p^2 = E^2 - m^2
  double pCM_1 = sqrt(ECM_1*ECM_1 - m1*m1);
  double pCM_2 = sqrt(ECM_2*ECM_2 - m2*m2);
  double pCM_3 = sqrt(ECM_3*ECM_3 - m3*m3);
  double pCM_4 = sqrt(ECM_4*ECM_4 - m4*m4);

  TVector3 fImpulsionLab_1 = TVector3(0,0,sqrt(fBeamEnergy*fBeamEnergy + 2*fBeamEnergy*m1));
  TVector3 fImpulsionLab_2 = TVector3(0,0,0);

  TLorentzVector fEnergyImpulsionLab_1 = TLorentzVector(fImpulsionLab_1,m1+fBeamEnergy);
  TLorentzVector fEnergyImpulsionLab_2 = TLorentzVector(fImpulsionLab_2,m2);

  TLorentzVector fTotalEnergyImpulsionLab = fEnergyImpulsionLab_1 + fEnergyImpulsionLab_2;

  double BetaCM = fTotalEnergyImpulsionLab.Beta();

  TLorentzVector fEnergyImpulsionCM_1 = fEnergyImpulsionLab_1;
  fEnergyImpulsionCM_1.Boost(0,0,-BetaCM);

  TLorentzVector fEnergyImpulsionCM_2 = fEnergyImpulsionLab_2;
  fEnergyImpulsionCM_2.Boost(0,0,-BetaCM);

  double thetaCM = 0;
  double rapidity = TMath::ATanH(BetaCM);
  if(pCM_3 - pCM_4 > 1e-6)cout << "I don't understand this: " << pCM_3 << "\t" << pCM_4 << endl;
  double p3 = (sqrt(m3*m3 + pCM_3*pCM_3) * cos(thetaLab) * sinh(rapidity) + cosh(rapidity) * sqrt(pCM_3*pCM_3 - m3*m3*sin(thetaLab)*sin(thetaLab)*sinh(rapidity)*sinh(rapidity))) / (1.0 + sin(thetaLab)*sin(thetaLab) * sinh(rapidity)*sinh(rapidity));
  thetaCM = asin(p3 * sin(thetaLab) / pCM_3);
  //cout << "thetaCM: " << TMath::Pi() << endl;
//   cout << "thetaCM: " << thetaCM*180./TMath::Pi() << endl;

  TLorentzVector fEnergyImpulsionCM_3  = TLorentzVector(pCM_3*sin(thetaCM),0,pCM_3*cos(thetaCM),ECM_3);
  TLorentzVector fEnergyImpulsionCM_4  = fTotalEnergyImpulsionCM - fEnergyImpulsionCM_3;

  TLorentzVector fEnergyImpulsionLab_3 = fEnergyImpulsionCM_3;
  fEnergyImpulsionLab_3.Boost(0,0,BetaCM);
  TLorentzVector fEnergyImpulsionLab_4 = fEnergyImpulsionCM_4;
  fEnergyImpulsionLab_4.Boost(0,0,BetaCM);

  // Angle in the lab frame
  double ThetaLab3 = fEnergyImpulsionLab_3.Angle(fEnergyImpulsionLab_1.Vect());
  if (ThetaLab3 < 0) ThetaLab3 += M_PI;
//   cout << "ThetaLab3: " << ThetaLab3*180./TMath::Pi() << endl;

  double ThetaLab4 = fEnergyImpulsionLab_4.Angle(fEnergyImpulsionLab_1.Vect());
  if (fabs(ThetaLab3) < 1e-6) ThetaLab3 = 0;
  ThetaLab4 = fabs(ThetaLab4);
  if (fabs(ThetaLab4) < 1e-6) ThetaLab4 = 0;

  // Kinetic Energy in the lab frame
  double KineticEnergyLab3 = fEnergyImpulsionLab_3.E() - m3;
  double KineticEnergyLab4 = fEnergyImpulsionLab_4.E() - m4;
//   cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;

  // test for total energy conversion
  if (fabs(fTotalEnergyImpulsionLab.E() - (fEnergyImpulsionLab_3.E()+fEnergyImpulsionLab_4.E())) > 1e-6)
    cout << "Problem for energy conservation" << endl;

  KineticEnergyLab3 = ELossSiO2(KineticEnergyLab3,0.5*SiO2Thickness);
  KineticEnergyLab3 = ELossC(KineticEnergyLab3,CarbonThickness);
//   cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;
//   cout << "Brho: " << EpToBrho(KineticEnergyLab3,m1) << endl;

  return EpToBrho(KineticEnergyLab3,m1);
}

double ExToBrhoNaF_Na(double fBeamEnergy, double Ex, double thetaLab)
{
  //print initial beam energy, calculate beam energy after target and print out
  //cout << "Initial Beam Energy: " << fBeamEnergy << endl;
    thetaLab *= TMath::Pi()/180.; //converting angle into radians
    double CarbonThickness = 0.08849557522;
    double NaFThickness = 0.1953125;
    fBeamEnergy = ELossNaF(fBeamEnergy,0.5*NaFThickness/cos(thetaLab)); //Proton energy after going through SiO2 target
  //cout << "Beam energy after target: " << fBeamEnergy << endl;

  //defining masses
  double m1 = 938.782980; //proton
  double m2 = mNa;        // Sodium
  double m3 = m1;         //proton
  double m4 = m2 + Ex;    //SiO2 + excitation energy
  double s = m1*m1 + m2*m2 + 2*m2*(fBeamEnergy + m1);

  TLorentzVector fTotalEnergyImpulsionCM = TLorentzVector(0,0,0,sqrt(s));

  //Calculating energies in center of mass reference frame
  double ECM_1 = (s + m1*m1 - m2*m2)/(2*sqrt(s));
  double ECM_2 = (s + m2*m2 - m1*m1)/(2*sqrt(s));
  double ECM_3 = (s + m3*m3 - m4*m4)/(2*sqrt(s));
  double ECM_4 = (s + m4*m4 - m3*m3)/(2*sqrt(s));

  //Calculating momentum in center of mass frame with formula p^2 = E^2 - m^2
  double pCM_1 = sqrt(ECM_1*ECM_1 - m1*m1);
  double pCM_2 = sqrt(ECM_2*ECM_2 - m2*m2);
  double pCM_3 = sqrt(ECM_3*ECM_3 - m3*m3);
  double pCM_4 = sqrt(ECM_4*ECM_4 - m4*m4);

  TVector3 fImpulsionLab_1 = TVector3(0,0,sqrt(fBeamEnergy*fBeamEnergy + 2*fBeamEnergy*m1));
  TVector3 fImpulsionLab_2 = TVector3(0,0,0);

  TLorentzVector fEnergyImpulsionLab_1 = TLorentzVector(fImpulsionLab_1,m1+fBeamEnergy);
  TLorentzVector fEnergyImpulsionLab_2 = TLorentzVector(fImpulsionLab_2,m2);

  TLorentzVector fTotalEnergyImpulsionLab = fEnergyImpulsionLab_1 + fEnergyImpulsionLab_2;

  double BetaCM = fTotalEnergyImpulsionLab.Beta();

  TLorentzVector fEnergyImpulsionCM_1 = fEnergyImpulsionLab_1;
  fEnergyImpulsionCM_1.Boost(0,0,-BetaCM);

  TLorentzVector fEnergyImpulsionCM_2 = fEnergyImpulsionLab_2;
  fEnergyImpulsionCM_2.Boost(0,0,-BetaCM);

  double thetaCM = 0;
  double rapidity = TMath::ATanH(BetaCM);
  if(pCM_3 - pCM_4 > 1e-6)cout << "I don't understand this: " << pCM_3 << "\t" << pCM_4 << endl;
  double p3 = (sqrt(m3*m3 + pCM_3*pCM_3) * cos(thetaLab) * sinh(rapidity) + cosh(rapidity) * sqrt(pCM_3*pCM_3 - m3*m3*sin(thetaLab)*sin(thetaLab)*sinh(rapidity)*sinh(rapidity))) / (1.0 + sin(thetaLab)*sin(thetaLab) * sinh(rapidity)*sinh(rapidity));
  thetaCM = asin(p3 * sin(thetaLab) / pCM_3);
//   cout << "thetaCM: " << thetaCM*180./TMath::Pi() << endl;

  TLorentzVector fEnergyImpulsionCM_3  = TLorentzVector(pCM_3*sin(thetaCM),0,pCM_3*cos(thetaCM),ECM_3);
  TLorentzVector fEnergyImpulsionCM_4  = fTotalEnergyImpulsionCM - fEnergyImpulsionCM_3;

  TLorentzVector fEnergyImpulsionLab_3 = fEnergyImpulsionCM_3;
  fEnergyImpulsionLab_3.Boost(0,0,BetaCM);
  TLorentzVector fEnergyImpulsionLab_4 = fEnergyImpulsionCM_4;
  fEnergyImpulsionLab_4.Boost(0,0,BetaCM);

  // Angle in the lab frame
  double ThetaLab3 = fEnergyImpulsionLab_3.Angle(fEnergyImpulsionLab_1.Vect());
  if (ThetaLab3 < 0) ThetaLab3 += M_PI;
//   cout << "ThetaLab3: " << ThetaLab3*180./TMath::Pi() << endl;

  double ThetaLab4 = fEnergyImpulsionLab_4.Angle(fEnergyImpulsionLab_1.Vect());
  if (fabs(ThetaLab3) < 1e-6) ThetaLab3 = 0;
  ThetaLab4 = fabs(ThetaLab4);
  if (fabs(ThetaLab4) < 1e-6) ThetaLab4 = 0;

  // Kinetic Energy in the lab frame
  double KineticEnergyLab3 = fEnergyImpulsionLab_3.E() - m3;
  double KineticEnergyLab4 = fEnergyImpulsionLab_4.E() - m4;
//   cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;

  // test for total energy conversion
  if (fabs(fTotalEnergyImpulsionLab.E() - (fEnergyImpulsionLab_3.E()+fEnergyImpulsionLab_4.E())) > 1e-6)
    cout << "Problem for energy conservation" << endl;

  KineticEnergyLab3 = ELossNaF(KineticEnergyLab3,0.5*NaFThickness);
  KineticEnergyLab3 = ELossC(KineticEnergyLab3,CarbonThickness);
//   cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;
//   cout << "Brho: " << EpToBrho(KineticEnergyLab3,m1) << endl;

  return EpToBrho(KineticEnergyLab3,m1);
}

double ExToBrhoNaF_F(double fBeamEnergy, double Ex, double thetaLab)
{
  //print initial beam energy, calculate beam energy after target and print out
//   cout << "Initial Beam Energy: " << fBeamEnergy << endl;
    thetaLab *= TMath::Pi()/180.; //converting angle into radians
    double CarbonThickness = 0.08849557522;
    double NaFThickness = 0.1953125;
    fBeamEnergy = ELossNaF(fBeamEnergy,0.5*NaFThickness/cos(thetaLab)); //Proton energy after going through SiO2 target
//   cout << "Beam energy after target: " << fBeamEnergy << endl;

  //defining masses
  double m1 = 938.782980; //proton
  double m2 = mF;         //fluorine
  double m3 = m1;         //proton
  double m4 = m2 + Ex;    //SiO2 + excitation energy
  double s = m1*m1 + m2*m2 + 2*m2*(fBeamEnergy + m1);

  TLorentzVector fTotalEnergyImpulsionCM = TLorentzVector(0,0,0,sqrt(s));

  //Calculating energies in center of mass reference frame
  double ECM_1 = (s + m1*m1 - m2*m2)/(2*sqrt(s));
  double ECM_2 = (s + m2*m2 - m1*m1)/(2*sqrt(s));
  double ECM_3 = (s + m3*m3 - m4*m4)/(2*sqrt(s));
  double ECM_4 = (s + m4*m4 - m3*m3)/(2*sqrt(s));

  //Calculating momentum in center of mass frame with formula p^2 = E^2 - m^2
  double pCM_1 = sqrt(ECM_1*ECM_1 - m1*m1);
  double pCM_2 = sqrt(ECM_2*ECM_2 - m2*m2);
  double pCM_3 = sqrt(ECM_3*ECM_3 - m3*m3);
  double pCM_4 = sqrt(ECM_4*ECM_4 - m4*m4);

  TVector3 fImpulsionLab_1 = TVector3(0,0,sqrt(fBeamEnergy*fBeamEnergy + 2*fBeamEnergy*m1));
  TVector3 fImpulsionLab_2 = TVector3(0,0,0);

  TLorentzVector fEnergyImpulsionLab_1 = TLorentzVector(fImpulsionLab_1,m1+fBeamEnergy);
  TLorentzVector fEnergyImpulsionLab_2 = TLorentzVector(fImpulsionLab_2,m2);

  TLorentzVector fTotalEnergyImpulsionLab = fEnergyImpulsionLab_1 + fEnergyImpulsionLab_2;

  double BetaCM = fTotalEnergyImpulsionLab.Beta();

  TLorentzVector fEnergyImpulsionCM_1 = fEnergyImpulsionLab_1;
  fEnergyImpulsionCM_1.Boost(0,0,-BetaCM);

  TLorentzVector fEnergyImpulsionCM_2 = fEnergyImpulsionLab_2;
  fEnergyImpulsionCM_2.Boost(0,0,-BetaCM);

  double thetaCM = 0;
  double rapidity = TMath::ATanH(BetaCM);
  if(pCM_3 - pCM_4 > 1e-6)cout << "I don't understand this: " << pCM_3 << "\t" << pCM_4 << endl;
  double p3 = (sqrt(m3*m3 + pCM_3*pCM_3) * cos(thetaLab) * sinh(rapidity) + cosh(rapidity) * sqrt(pCM_3*pCM_3 - m3*m3*sin(thetaLab)*sin(thetaLab)*sinh(rapidity)*sinh(rapidity))) / (1.0 + sin(thetaLab)*sin(thetaLab) * sinh(rapidity)*sinh(rapidity));
  thetaCM = asin(p3 * sin(thetaLab) / pCM_3);
//   cout << "thetaCM: " << thetaCM*180./TMath::Pi() << endl;

  TLorentzVector fEnergyImpulsionCM_3  = TLorentzVector(pCM_3*sin(thetaCM),0,pCM_3*cos(thetaCM),ECM_3);
  TLorentzVector fEnergyImpulsionCM_4  = fTotalEnergyImpulsionCM - fEnergyImpulsionCM_3;

  TLorentzVector fEnergyImpulsionLab_3 = fEnergyImpulsionCM_3;
  fEnergyImpulsionLab_3.Boost(0,0,BetaCM);
  TLorentzVector fEnergyImpulsionLab_4 = fEnergyImpulsionCM_4;
  fEnergyImpulsionLab_4.Boost(0,0,BetaCM);

  // Angle in the lab frame
  double ThetaLab3 = fEnergyImpulsionLab_3.Angle(fEnergyImpulsionLab_1.Vect());
  if (ThetaLab3 < 0) ThetaLab3 += M_PI;
//   cout << "ThetaLab3: " << ThetaLab3*180./TMath::Pi() << endl;

  double ThetaLab4 = fEnergyImpulsionLab_4.Angle(fEnergyImpulsionLab_1.Vect());
  if (fabs(ThetaLab3) < 1e-6) ThetaLab3 = 0;
  ThetaLab4 = fabs(ThetaLab4);
  if (fabs(ThetaLab4) < 1e-6) ThetaLab4 = 0;

  // Kinetic Energy in the lab frame
  double KineticEnergyLab3 = fEnergyImpulsionLab_3.E() - m3;
  double KineticEnergyLab4 = fEnergyImpulsionLab_4.E() - m4;
//   cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;

  // test for total energy conversion
  if (fabs(fTotalEnergyImpulsionLab.E() - (fEnergyImpulsionLab_3.E()+fEnergyImpulsionLab_4.E())) > 1e-6)
    cout << "Problem for energy conservation" << endl;

  KineticEnergyLab3 = ELossNaF(KineticEnergyLab3,0.5*NaFThickness);
  KineticEnergyLab3 = ELossC(KineticEnergyLab3,CarbonThickness);
//   cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;
//   cout << "Brho: " << EpToBrho(KineticEnergyLab3,m1) << endl;

  return EpToBrho(KineticEnergyLab3,m1);
}

double ExToBrhoLiF_Li(double fBeamEnergy, double Ex, double thetaLab)
{
  //print initial beam energy, calculate beam energy after target and print out
//   cout << "Initial Beam Energy: " << fBeamEnergy << endl;
  thetaLab *= TMath::Pi()/180.; //angle in radians
  double CarbonThickness = 0.08849557522;
  double LiFThickness = 0.31818181818;
  fBeamEnergy = ELossLiF(fBeamEnergy,0.5*LiFThickness/cos(thetaLab)); //Proton energy after going through SiO2 target
//   cout << "Beam energy after target: " << fBeamEnergy << endl;

  //declaring masses
  double m1 = 938.782980;  //proton
  double m2 = mLi;         //lithium
  double m3 = m1;          //Proton
  double m4 = m2 + Ex;     //LiF + Excitation energy
  double s = m1*m1 + m2*m2 + 2*m2*(fBeamEnergy + m1);
  TLorentzVector fTotalEnergyImpulsionCM = TLorentzVector(0,0,0,sqrt(s));

  //Energy in Center of Mass reference frame
  double ECM_1 = (s + m1*m1 - m2*m2)/(2*sqrt(s));
  double ECM_2 = (s + m2*m2 - m1*m1)/(2*sqrt(s));
  double ECM_3 = (s + m3*m3 - m4*m4)/(2*sqrt(s));
  double ECM_4 = (s + m4*m4 - m3*m3)/(2*sqrt(s));

  //Momentum in Center of Mass reference frame
  double pCM_1 = sqrt(ECM_1*ECM_1 - m1*m1);
  double pCM_2 = sqrt(ECM_2*ECM_2 - m2*m2);
  double pCM_3 = sqrt(ECM_3*ECM_3 - m3*m3);
  double pCM_4 = sqrt(ECM_4*ECM_4 - m4*m4);

  TVector3 fImpulsionLab_1 = TVector3(0,0,sqrt(fBeamEnergy*fBeamEnergy + 2*fBeamEnergy*m1));
  TVector3 fImpulsionLab_2 = TVector3(0,0,0);

  TLorentzVector fEnergyImpulsionLab_1 = TLorentzVector(fImpulsionLab_1,m1+fBeamEnergy);
  TLorentzVector fEnergyImpulsionLab_2 = TLorentzVector(fImpulsionLab_2,m2);

  TLorentzVector fTotalEnergyImpulsionLab = fEnergyImpulsionLab_1 + fEnergyImpulsionLab_2;

  double BetaCM = fTotalEnergyImpulsionLab.Beta();

  TLorentzVector fEnergyImpulsionCM_1 = fEnergyImpulsionLab_1;
  fEnergyImpulsionCM_1.Boost(0,0,-BetaCM);

  TLorentzVector fEnergyImpulsionCM_2 = fEnergyImpulsionLab_2;
  fEnergyImpulsionCM_2.Boost(0,0,-BetaCM);

  double thetaCM = 0;
  double rapidity = TMath::ATanH(BetaCM);
  if(pCM_3 - pCM_4 > 1e-6)cout << "I don't understand this: " << pCM_3 << "\t" << pCM_4 << endl;
  double p3 = (sqrt(m3*m3 + pCM_3*pCM_3) * cos(thetaLab) * sinh(rapidity) + cosh(rapidity) * sqrt(pCM_3*pCM_3 - m3*m3*sin(thetaLab)*sin(thetaLab)*sinh(rapidity)*sinh(rapidity))) / (1.0 + sin(thetaLab)*sin(thetaLab) * sinh(rapidity)*sinh(rapidity));
  thetaCM = asin(p3 * sin(thetaLab) / pCM_3);
//   cout << "thetaCM: " << thetaCM*180./TMath::Pi() << endl;

  TLorentzVector fEnergyImpulsionCM_3  = TLorentzVector(pCM_3*sin(thetaCM),0,pCM_3*cos(thetaCM),ECM_3);
  TLorentzVector fEnergyImpulsionCM_4  = fTotalEnergyImpulsionCM - fEnergyImpulsionCM_3;

  TLorentzVector fEnergyImpulsionLab_3 = fEnergyImpulsionCM_3;
  fEnergyImpulsionLab_3.Boost(0,0,BetaCM);
  TLorentzVector fEnergyImpulsionLab_4 = fEnergyImpulsionCM_4;
  fEnergyImpulsionLab_4.Boost(0,0,BetaCM);

  // Angle in the lab frame
  double ThetaLab3 = fEnergyImpulsionLab_3.Angle(fEnergyImpulsionLab_1.Vect());
  if (ThetaLab3 < 0) ThetaLab3 += M_PI;
//   cout << "ThetaLab3: " << ThetaLab3*180./TMath::Pi() << endl;

  double ThetaLab4 = fEnergyImpulsionLab_4.Angle(fEnergyImpulsionLab_1.Vect());
  if (fabs(ThetaLab3) < 1e-6) ThetaLab3 = 0;
  ThetaLab4 = fabs(ThetaLab4);
  if (fabs(ThetaLab4) < 1e-6) ThetaLab4 = 0;

  // Kinetic Energy in the lab frame
  double KineticEnergyLab3 = fEnergyImpulsionLab_3.E() - m3;
  double KineticEnergyLab4 = fEnergyImpulsionLab_4.E() - m4;
//   cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;

  // test for total energy conversion
  if (fabs(fTotalEnergyImpulsionLab.E() - (fEnergyImpulsionLab_3.E()+fEnergyImpulsionLab_4.E())) > 1e-6)
    cout << "Problem for energy conservation" << endl;

  KineticEnergyLab3 = ELossLiF(KineticEnergyLab3,0.5*LiFThickness);
  KineticEnergyLab3 = ELossC(KineticEnergyLab3,CarbonThickness);
//   cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;
//   cout << "Brho: " << EpToBrho(KineticEnergyLab3,m1) << endl;

  return EpToBrho(KineticEnergyLab3,m1);
}

double ExToBrhoLiF_F(double fBeamEnergy, double Ex, double thetaLab)
{
  //print initial beam energy, calculate beam energy after target and print out
//   cout << "Initial Beam Energy: " << fBeamEnergy << endl;
  thetaLab *= TMath::Pi()/180.; //angle in radians
  double CarbonThickness = 0.08849557522;
  double LiFThickness = 0.31818181818;
  fBeamEnergy = ELossLiF(fBeamEnergy,0.5*LiFThickness/cos(thetaLab)); //Proton energy after going through SiO2 target
//   cout << "Beam energy after target: " << fBeamEnergy << endl;

  //declaring masses
  double m1 = 938.782980;  //proton
  double m2 = mF;          //Fluorine
  double m3 = m1;          //Proton
  double m4 = m2 + Ex;     //LiF + Excitation energy
  double s = m1*m1 + m2*m2 + 2*m2*(fBeamEnergy + m1);
  TLorentzVector fTotalEnergyImpulsionCM = TLorentzVector(0,0,0,sqrt(s));

  //Energy in Center of Mass reference frame
  double ECM_1 = (s + m1*m1 - m2*m2)/(2*sqrt(s));
  double ECM_2 = (s + m2*m2 - m1*m1)/(2*sqrt(s));
  double ECM_3 = (s + m3*m3 - m4*m4)/(2*sqrt(s));
  double ECM_4 = (s + m4*m4 - m3*m3)/(2*sqrt(s));

  //Momentum in Center of Mass reference frame
  double pCM_1 = sqrt(ECM_1*ECM_1 - m1*m1);
  double pCM_2 = sqrt(ECM_2*ECM_2 - m2*m2);
  double pCM_3 = sqrt(ECM_3*ECM_3 - m3*m3);
  double pCM_4 = sqrt(ECM_4*ECM_4 - m4*m4);

  TVector3 fImpulsionLab_1 = TVector3(0,0,sqrt(fBeamEnergy*fBeamEnergy + 2*fBeamEnergy*m1));
  TVector3 fImpulsionLab_2 = TVector3(0,0,0);

  TLorentzVector fEnergyImpulsionLab_1 = TLorentzVector(fImpulsionLab_1,m1+fBeamEnergy);
  TLorentzVector fEnergyImpulsionLab_2 = TLorentzVector(fImpulsionLab_2,m2);

  TLorentzVector fTotalEnergyImpulsionLab = fEnergyImpulsionLab_1 + fEnergyImpulsionLab_2;

  double BetaCM = fTotalEnergyImpulsionLab.Beta();

  TLorentzVector fEnergyImpulsionCM_1 = fEnergyImpulsionLab_1;
  fEnergyImpulsionCM_1.Boost(0,0,-BetaCM);

  TLorentzVector fEnergyImpulsionCM_2 = fEnergyImpulsionLab_2;
  fEnergyImpulsionCM_2.Boost(0,0,-BetaCM);

  double thetaCM = 0;
  double rapidity = TMath::ATanH(BetaCM);
  if(pCM_3 - pCM_4 > 1e-6)cout << "I don't understand this: " << pCM_3 << "\t" << pCM_4 << endl;
  double p3 = (sqrt(m3*m3 + pCM_3*pCM_3) * cos(thetaLab) * sinh(rapidity) + cosh(rapidity) * sqrt(pCM_3*pCM_3 - m3*m3*sin(thetaLab)*sin(thetaLab)*sinh(rapidity)*sinh(rapidity))) / (1.0 + sin(thetaLab)*sin(thetaLab) * sinh(rapidity)*sinh(rapidity));
  thetaCM = asin(p3 * sin(thetaLab) / pCM_3);
  cout << "thetaCM: " << thetaCM*180./TMath::Pi() << endl;

  TLorentzVector fEnergyImpulsionCM_3  = TLorentzVector(pCM_3*sin(thetaCM),0,pCM_3*cos(thetaCM),ECM_3);
  TLorentzVector fEnergyImpulsionCM_4  = fTotalEnergyImpulsionCM - fEnergyImpulsionCM_3;

  TLorentzVector fEnergyImpulsionLab_3 = fEnergyImpulsionCM_3;
  fEnergyImpulsionLab_3.Boost(0,0,BetaCM);
  TLorentzVector fEnergyImpulsionLab_4 = fEnergyImpulsionCM_4;
  fEnergyImpulsionLab_4.Boost(0,0,BetaCM);

  // Angle in the lab frame
  double ThetaLab3 = fEnergyImpulsionLab_3.Angle(fEnergyImpulsionLab_1.Vect());
  if (ThetaLab3 < 0) ThetaLab3 += M_PI;
  cout << "ThetaLab3: " << ThetaLab3*180./TMath::Pi() << endl;

  double ThetaLab4 = fEnergyImpulsionLab_4.Angle(fEnergyImpulsionLab_1.Vect());
  if (fabs(ThetaLab3) < 1e-6) ThetaLab3 = 0;
  ThetaLab4 = fabs(ThetaLab4);
  if (fabs(ThetaLab4) < 1e-6) ThetaLab4 = 0;

  // Kinetic Energy in the lab frame
  double KineticEnergyLab3 = fEnergyImpulsionLab_3.E() - m3;
  double KineticEnergyLab4 = fEnergyImpulsionLab_4.E() - m4;
//   cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;

  // test for total energy conversion
  if (fabs(fTotalEnergyImpulsionLab.E() - (fEnergyImpulsionLab_3.E()+fEnergyImpulsionLab_4.E())) > 1e-6)
    cout << "Problem for energy conservation" << endl;

  KineticEnergyLab3 = ELossLiF(KineticEnergyLab3,0.5*LiFThickness);
  KineticEnergyLab3 = ELossC(KineticEnergyLab3,CarbonThickness);
//   cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;
//   cout << "Brho: " << EpToBrho(KineticEnergyLab3,m1) << endl;

  return EpToBrho(KineticEnergyLab3,m1);
}


//-------------------------------------------------------------------------------------------------------------//

double BrhoToExSiO2(double Brho, double mass, double theta)
{
  //Calculating proton energy throughout targets
  double Ep = sqrt(pow(TMath::C()/1e6*Brho,2.) + pow(mass,2.)) - mass;
   //cout << "Proton energy from focal plane: " << Ep << endl;
  double CarbonThickness = 0.08849557522;
  double SiO2Thickness = 0.15094339622;
  Ep = ELossC(Ep,-CarbonThickness); //Proton energy after going through Carbon target (- because opposite direction)
  Ep = ELossSiO2(Ep,-0.5*SiO2Thickness); //Proton energy after going through SiO2 target
   //cout << "Proton energy after target energy-loss correction: " << Ep << endl;

  double p3 = sqrt(Ep * (Ep + 2*mass));
  double T1 = 14.000;
  theta *= TMath::Pi()/180.;
  T1 = ELossSiO2(T1,0.5*SiO2Thickness/cos(theta));

  double p1 = sqrt(T1 * (T1 + 2*mass));
  double p4 = sqrt(p1*p1 - 2*p1*p3*cos(theta) + p3*p3);
//   cout << "p4: " << p4 << endl;

  double m_28Si = mSi; //MeV/c/c // silicon
  //double m_23Na = mNa;   //MeV/c/c // sodium
  //double m_28F = mF;  //MeV/c/c // fluorine

  double T4 = sqrt(p4*p4 + m_28Si*m_28Si) - m_28Si;
  //cout << "T4: " << T4 << endl;

  double Ex = sqrt(pow(T1 - Ep + m_28Si,2.) - p4*p4) - m_28Si;

  return Ex;
}


double BrhoToExNaF(double Brho, double mass, double theta)
{
  //Calculating proton energy throughout targets
  double Ep = sqrt(pow(TMath::C()/1e6*Brho,2.) + pow(mass,2.)) - mass;
//    cout << "Proton energy from focal plane: " << Ep << endl;
   double CarbonThickness = 0.08849557522;
   double NaFThickness = 0.1953125;
   Ep = ELossC(Ep,-CarbonThickness);
   Ep = ELossNaF(Ep,-0.5*NaFThickness);
//    cout << "Proton energy after target energy-loss correction: " << Ep << endl;

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


double BrhoToExLiF(double Brho, double mass, double theta)
{
  double Ep = sqrt(pow(TMath::C()/1e6*Brho,2.) + pow(mass,2.)) - mass;
//    cout << "Proton energy from focal plane: " << Ep << endl;
   double CarbonThickness = 0.08849557522;
   double LiFThickness = 0.31818181818;
   Ep = ELossC(Ep,-CarbonThickness);
   Ep = ELossLiF(Ep,-0.5*LiFThickness);
//    cout << "Proton energy after target energy-loss correction: " << Ep << endl;

  double p3 = sqrt(Ep * (Ep + 2*mass));

  double T1 = 14.000;
  theta *= TMath::Pi()/180.;
  T1 = ELossLiF(T1,0.5*LiFThickness/cos(theta));

  double p1 = sqrt(T1 * (T1 + 2*mass));

  double p4 = sqrt(p1*p1 - 2*p1*p3*cos(theta) + p3*p3);
//   cout << "p4: " << p4 << endl;

  //double m_28Si = 26060.33946; //MeV/c/c // silicon
  double m_23Na = mF; //MeV/c/c // sodium
  //double m_recoil = 17696.90035; //MeV/c/c // fluorine
  //m_recoil = m_23Na;

  double T4 = sqrt(p4*p4 + m_23Na*m_23Na) - m_23Na;
//   cout << "T4: " << T4 << endl;

  double Ex = sqrt(pow(T1 - Ep + m_23Na,2.) - p4*p4) - m_23Na;

  return Ex;
}
