//----------------------------------------------------------------------------
// Implementation of the KFParticle class
// .
// @author  I.Kisel, I.Kulakov, M.Zyzak
// @version 1.0
// @since   20.08.13
// 
// 
//  -= Copyright &copy ALICE HLT and CBM L1 Groups =-
//____________________________________________________________________________

#ifndef KFPartEfficiencies_H
#define KFPartEfficiencies_H

#ifndef HLTCA_STANDALONE
#include "TNamed.h"
#endif

#include <map>
#include <iomanip>
#include "KFMCCounter.h"

class KFPartEfficiencies: public TNamed
{
 public:

  KFPartEfficiencies():
    names(),
    indices(),
    ratio_reco1(),
    ratio_reco2(),
    ratio_reco3(),
    mc1(),
    mc2(),
    mc3(),
    reco(),
    ratio_ghost(),
    ratio_bg(),
    ratio_clone(),
    ghost(),
    bg(),
    clone()
  {
    int mPartPDG[nParticles] = {310,3122,-3122,3312,-3312,3322,-3322,3334,-3334,3212,-3212,3222,-3222, //strange meson and hyperons
                                313,-313,323,-323, 100313, 100323, -100323,//K* resonances
                                3224,3114,-3114,-3224, 3214, -3214,//sigma resonances
                                3124,-3124, //Lambda resonances
                                3324, -3324, 1003314, -1003314, 3314, -3314, //Xi resonances
                                1003334, -1003334, //Omega resonances
                                3000, //exotics
                                333,113, //vector mesons, hadron chanel
                                100113, 200113, //light vector mesons
                                22, //dielectrons
                                111,221, //pi0, eta
                                443,100443, // J/Psi
                                421,-421,100421,-100421, //D0
                                411,-411, //D+, D-
                                431,-431, //Ds+, Ds-
                                4122,-4122, //Lambdac
                                10421, -10421, 10411, -10411, 20411, -20411,
                                3001, //H->Lambda p pi
                                11, -11, 13, -13, 211,  -211, 321, -321, 2212, -2212, // stable particles
                                123456789 //V0
                               };
    TString mPartName[nParticles] = {"ks","lambda","lambdab","xi-","xi+","xi0","xi0b","omega-","omega+","#Sigma^0","#Sigma^0b", "#Sigma^+", "#Sigma^+b",
                                     "k*0","k*0b","k*+","k*-", "k*0_{K0,#pi0}", "k*+_{K+,#pi0}", "k*-_{K-,#pi0}",
                                     "sigma*+","sigma*-","sigma*+b","sigma*-b","sigma*0","sigma*0b",
                                     "lambda*","lambda*b",
                                     "xi*0", "xi*0b", "xi*-_{#Lambda,K}", "xi*+_{#Lambda,K}", "xi*-_{#xi-,#pi0}", "xi*+_{#xi+,#pi0}",
                                     "omega*-","omega*+",
                                     "Hdb",
                                     "phi_{KK}", "rho_{#pi#pi}",
                                     "rho_{ee}", "rho_{#mu#mu}",
                                     "gamma",
                                     "#pi^{0}","eta",
                                     "J#Psi_ee","J#Psi_#mu#mu",
                                     "D0","D0b","D0_4","D0b_4",
                                     "D+","D-",
                                     "Ds+","Ds-",
                                     "lambdac", "lambdacb",
                                     "D*0", "D*0b", "D*+", "D*-", "D*+_4", "D*-_4",
                                     "H0",
                                     "e-", "e+", "mu-", "mu+", "pi+", "pi-", "K+", "K-", "p+", "p-",
                                     "V0"
                                    };
    TString mPartTitle[nParticles] = {"KShort   ", //0
                                      "Lambda   ", //1
                                      "Lambda b ", //2
                                      "Xi-      ", //3
                                      "Xi+      ", //4
                                      "Xi0      ", //5
                                      "Xi0 b    ", //6
                                      "Omega-   ", //7
                                      "Omega+   ", //8
                                      "Sigma0   ", //9
                                      "Sigma0 b ", //10
                                      "Sigma+   ", //11
                                      "Sigma+ b ", //12
                                      "K*0      ", //13
                                      "K*0 b    ", //14
                                      "K*+      ", //15
                                      "K*-      ", //16
                                      "K*0_K0pi0", //17
                                      "K*+_K+pi0", //18
                                      "K*-_K-pi0", //19
                                      "Sigma*+  ", //20
                                      "Sigma*-  ", //21
                                      "Sigma*+ b", //22
                                      "Sigma*- b", //23
                                      "Sigma*0  ", //24
                                      "Sigma*0 b", //25
                                      "Lambda*  ", //26
                                      "Lambda* b", //27
                                      "Xi*0     ", //28
                                      "Xi*0 b   ", //29
                                      "Xi*-_lk  ", //30
                                      "Xi*+_lk  ", //31
                                      "Xi*-_XiPi", //32
                                      "Xi*+_XiPi", //33
                                      "Omega*-  ", //34
                                      "Omega*+  ", //35
                                      "Hdb      ", //36
                                      "phi_kk   ", //37
                                      "rho_pipi ", //38
                                      "rho_ee   ", //39
                                      "rho_mm   ", //40
                                      "gamma    ", //41
                                      "Pi0      ", //42
                                      "eta      ", //43
                                      "JPsi_ee  ", //44
                                      "JPsi_mm  ", //45
                                      "D0       ", //46
                                      "D0b      ", //47
                                      "D0_4     ", //48
                                      "D0b_4    ", //49
                                      "D+       ", //50
                                      "D-       ", //51
                                      "Ds+      ", //52
                                      "Ds-      ", //53
                                      "Lambdac  ", //54
                                      "Lambdac b", //55
                                      "D*0      ", //56
                                      "D*0 b    ", //57
                                      "D*+      ", //58
                                      "D*-      ", //59
                                      "D*+_4    ", //60
                                      "D*-_4    ", //61
                                      "H0       ", //62
                                      "e-       ", //63
                                      "e+       ", //64
                                      "mu-      ", //65
                                      "mu+      ", //66
                                      "pi+      ", //67
                                      "pi-      ", //68
                                      "K+       ", //69
                                      "K-       ", //70
                                      "p+       ", //71
                                      "p-       ", //72
                                      "V0       " //73
                                     };

    float mPartMHistoMin[nParticles] = {0.3, 1., 1., 1., 1., 1., 1.,1.,1.,1.,1.,1.,1.,
                                        0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
                                        1.,1.,1.,1.,1.,1.,
                                        1.4, 1.4,
                                        1.4, 1.4, 1.4, 1.4, 1.4, 1.4,
                                        1.8,1.8,
                                        1.,
                                        0.8, 0.1,
                                        0.1, 0.1,
                                        0.,
                                        0.,0.,
                                        1.,1.,
                                        1.,1.,1.,1.,
                                        1.,1.,
                                        1.,1.,
                                        1.8,1.8,
                                        1.8,1.8,1.8,1.8,1.8,1.8,
                                        1.,
                                        0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                        0.3 };
    float mPartMHistoMax[nParticles] = {1.3, 2., 2., 3., 3., 3., 3., 3., 3.,3.,3.,3.,3.,
                                        2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6,
                                        3., 3., 3., 3.,3.,3.,
                                        3.4, 3.4,
                                        3.4, 3.4, 3.4, 3.4, 3.4, 3.4,
                                        3.8, 3.8,
                                        3.,
                                        2.8, 2.1,
                                        2.1, 2.1,
                                        3.,
                                        3.,3.,
                                        4.,4.,
                                        3.,3.,3.,3.,
                                        3.,3.,
                                        3.,3.,
                                        3.8,3.8,
                                        3.8,3.8,3.8,3.8,3.8,3.8,
                                        3.,
                                        0.01, 0.01, 1., 1., 1., 1., 1., 1., 1.5, 1.5,
                                        1.3};
                                        
    int mPartMaxMult[nParticles];
    for(int i=0; i<71; i++)
      mPartMaxMult[i] = 20;
    mPartMaxMult[63] = 20;
    mPartMaxMult[64] = 20;
    mPartMaxMult[65] = 20;
    mPartMaxMult[66] = 20;
    mPartMaxMult[67] = 500;
    mPartMaxMult[68] = 500;
    mPartMaxMult[69] = 50;
    mPartMaxMult[70] = 50;
    mPartMaxMult[71] = 500;
    mPartMaxMult[72] = 20;
        
    //set decay mode
    partDaughterPdg.resize(nParticles);

    int curPart = 0;
    
    partDaughterPdg[curPart].push_back(  211); //K0s -> pi+ pi-
    partDaughterPdg[curPart].push_back( -211);
    curPart++;
    
    partDaughterPdg[curPart].push_back( 2212); //Lambda -> p pi-
    partDaughterPdg[curPart].push_back( -211);
    curPart++;
    
    partDaughterPdg[curPart].push_back(-2212); //Lambda_bar -> p- pi+
    partDaughterPdg[curPart].push_back(  211);
    curPart++;
    
    partDaughterPdg[curPart].push_back( 3122); //Xi- -> Lambda pi-
    partDaughterPdg[curPart].push_back( -211);
    curPart++;
    
    partDaughterPdg[curPart].push_back(-3122); //Xi+ -> Lambda_bar pi+
    partDaughterPdg[curPart].push_back(  211);
    curPart++;
    
    partDaughterPdg[curPart].push_back( 3122); //Xi0 -> Lambda pi0
    partDaughterPdg[curPart].push_back(  111);
    curPart++;
    
    partDaughterPdg[curPart].push_back(-3122); //Xi0_bar -> Lambda_bar pi0
    partDaughterPdg[curPart].push_back(  111);
    curPart++;
    
    partDaughterPdg[curPart].push_back( 3122); //Omega- -> Lambda K-
    partDaughterPdg[curPart].push_back( -321);
    curPart++;
    
    partDaughterPdg[curPart].push_back(-3122); //Omega+ -> Lambda_bar K+
    partDaughterPdg[curPart].push_back(  321);
    curPart++;
    
    partDaughterPdg[curPart].push_back(   22); //Sigma0 -> Lambda Gamma
    partDaughterPdg[curPart].push_back( 3122);
    curPart++;
    
    partDaughterPdg[curPart].push_back(   22); //Sigma0_bar -> Lambda_bar Gamma
    partDaughterPdg[curPart].push_back(-3122);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  111); //Sigma+ -> p pi0
    partDaughterPdg[curPart].push_back( 2212);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  111); //Sigma+_bar -> p- pi0
    partDaughterPdg[curPart].push_back(-2212);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  321); //K*0 -> K+ pi-
    partDaughterPdg[curPart].push_back( -211);
    curPart++;
    
    partDaughterPdg[curPart].push_back( -321); //K*0_bar -> K- pi+
    partDaughterPdg[curPart].push_back(  211);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  310); //K*+ -> K0s pi+
    partDaughterPdg[curPart].push_back(  211);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  310); //K*- -> K0s pi-
    partDaughterPdg[curPart].push_back( -211);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  310); //K*0 -> K0 pi0
    partDaughterPdg[curPart].push_back(  111);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  321); //K*+ -> K+ pi0
    partDaughterPdg[curPart].push_back(  111);
    curPart++;
    
    partDaughterPdg[curPart].push_back( -321); //K*- -> K- pi0
    partDaughterPdg[curPart].push_back(  111);
    curPart++;
    
    partDaughterPdg[curPart].push_back( 3122); //Sigma+ -> Lambda pi+
    partDaughterPdg[curPart].push_back(  211);
    curPart++;
    
    partDaughterPdg[curPart].push_back( 3122); //Sigma- -> Lambda pi-
    partDaughterPdg[curPart].push_back( -211);
    curPart++;
    
    partDaughterPdg[curPart].push_back(-3122); //Sigma+_bar -> Lambda_bar pi+
    partDaughterPdg[curPart].push_back(  211);
    curPart++;
    
    partDaughterPdg[curPart].push_back(-3122); //Sigma-_bar -> Lambda_bar pi-
    partDaughterPdg[curPart].push_back( -211);
    curPart++;
    
    partDaughterPdg[curPart].push_back( 3122); //Sigma*0 -> Lambda pi0
    partDaughterPdg[curPart].push_back(  111);
    curPart++;
    
    partDaughterPdg[curPart].push_back(-3122); //Sigma*0_bar -> Lambda_bar pi0
    partDaughterPdg[curPart].push_back(  111);
    curPart++;
    
    partDaughterPdg[curPart].push_back( 2212); //Lambda* -> p K-
    partDaughterPdg[curPart].push_back( -321);
    curPart++;
    
    partDaughterPdg[curPart].push_back(-2212); //Lambda*_bar -> p- K+
    partDaughterPdg[curPart].push_back(  321);
    curPart++;
    
    partDaughterPdg[curPart].push_back( 3312); //Xi*0 -> Xi- pi+
    partDaughterPdg[curPart].push_back(  211);
    curPart++;
    
    partDaughterPdg[curPart].push_back(-3312); //Xi*0_bar -> Xi+ pi-
    partDaughterPdg[curPart].push_back( -211);
    curPart++;
    
    partDaughterPdg[curPart].push_back( 3122); //Xi*- -> Lambda K-
    partDaughterPdg[curPart].push_back( -321);
    curPart++;
    
    partDaughterPdg[curPart].push_back(-3122); //Xi*+ -> Lambda_bar K+
    partDaughterPdg[curPart].push_back(  321);
    curPart++;
    
    partDaughterPdg[curPart].push_back( 3312); //Xi*- -> Xi- pi0
    partDaughterPdg[curPart].push_back(  111);
    curPart++;
    
    partDaughterPdg[curPart].push_back(-3312); //Xi*+ -> Xi+ pi0
    partDaughterPdg[curPart].push_back(  111);
    curPart++;
    
    partDaughterPdg[curPart].push_back( 3312); //Omega*- -> Xi- pi+ K-
    partDaughterPdg[curPart].push_back(  211);
    partDaughterPdg[curPart].push_back( -321);
    curPart++;
    
    partDaughterPdg[curPart].push_back(-3312); //Omega*- -> Xi+ pi- K+
    partDaughterPdg[curPart].push_back( -211);
    partDaughterPdg[curPart].push_back(  321);
    curPart++;
    
    partDaughterPdg[curPart].push_back( 3122); //H-dibar -> Lambda Lambda
    partDaughterPdg[curPart].push_back( 3122);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  321); //phi -> K+ K-
    partDaughterPdg[curPart].push_back( -321);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  211); //rho, omega, phi -> pi+ pi-
    partDaughterPdg[curPart].push_back( -211);
    curPart++;
    
    partDaughterPdg[curPart].push_back(   11); //rho, omega, phi -> e+ e-
    partDaughterPdg[curPart].push_back(  -11);
    curPart++;
    
    partDaughterPdg[curPart].push_back(   13); //rho, omega, phi -> mu+ mu-
    partDaughterPdg[curPart].push_back(  -13);
    curPart++;
    
    partDaughterPdg[curPart].push_back(   11); //gamma -> e+ e-
    partDaughterPdg[curPart].push_back(  -11);
    curPart++;
    
    partDaughterPdg[curPart].push_back(   22); //pi0 -> gamma gamma
    partDaughterPdg[curPart].push_back(   22);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  111); //eta -> pi0 pi0
    partDaughterPdg[curPart].push_back(  111);
    partDaughterPdg[curPart].push_back(  111);
    curPart++;
    
    partDaughterPdg[curPart].push_back(   11); //JPsi -> e+ e-
    partDaughterPdg[curPart].push_back(  -11);
    curPart++;
    
    partDaughterPdg[curPart].push_back(   13); //JPsi -> mu+ mu-
    partDaughterPdg[curPart].push_back(  -13);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  211); //D0 -> pi+ K-
    partDaughterPdg[curPart].push_back( -321);
    curPart++;
    
    partDaughterPdg[curPart].push_back( -211); //D0_bar -> K+ pi-
    partDaughterPdg[curPart].push_back(  321);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  211); //D0 -> pi+ pi+ pi- K-
    partDaughterPdg[curPart].push_back(  211);
    partDaughterPdg[curPart].push_back( -211);
    partDaughterPdg[curPart].push_back( -321);
    curPart++;
    
    partDaughterPdg[curPart].push_back( -211); //D0_bar -> pi- pi- pi+ K+
    partDaughterPdg[curPart].push_back( -211);
    partDaughterPdg[curPart].push_back(  211);
    partDaughterPdg[curPart].push_back(  321);
    curPart++;
    
    partDaughterPdg[curPart].push_back( -321); //D+ -> K- pi+ pi+
    partDaughterPdg[curPart].push_back(  211);
    partDaughterPdg[curPart].push_back(  211);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  321); //D- -> K+ pi- pi-
    partDaughterPdg[curPart].push_back( -211);
    partDaughterPdg[curPart].push_back( -211);
    curPart++;
    
    partDaughterPdg[curPart].push_back( -321); //Ds+ -> K- K+ pi+
    partDaughterPdg[curPart].push_back(  321);
    partDaughterPdg[curPart].push_back(  211);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  321); //Ds- -> K+ K- pi-
    partDaughterPdg[curPart].push_back( -321);
    partDaughterPdg[curPart].push_back( -211);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  211); //Lambdac -> pi+ K- p
    partDaughterPdg[curPart].push_back( -321);
    partDaughterPdg[curPart].push_back( 2212);
    curPart++;
    
    partDaughterPdg[curPart].push_back( -211); //Lambdac_bar -> pi- K+ p-
    partDaughterPdg[curPart].push_back(  321);
    partDaughterPdg[curPart].push_back(-2212);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  411); //D*0 -> D+ pi-
    partDaughterPdg[curPart].push_back( -211);
    curPart++;
    
    partDaughterPdg[curPart].push_back( -411); //D*0_bar -> D- pi+
    partDaughterPdg[curPart].push_back(  211);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  421); //D*+ -> D0 pi+
    partDaughterPdg[curPart].push_back(  211);
    curPart++;
    
    partDaughterPdg[curPart].push_back( -421); //D*- -> D0_bar pi-
    partDaughterPdg[curPart].push_back( -211);
    curPart++;
    
    partDaughterPdg[curPart].push_back(  421); //D*+ -> D04 pi+
    partDaughterPdg[curPart].push_back(  211);
    curPart++;
    
    partDaughterPdg[curPart].push_back( -421); //D*- -> D04_bar pi-
    partDaughterPdg[curPart].push_back( -211);
    curPart++;
    
    partDaughterPdg[curPart].push_back( 3122); //H0-> Lambda pi- p
    partDaughterPdg[curPart].push_back( -211);
    partDaughterPdg[curPart].push_back( 2212);
    curPart++;
    
    for(int iP=0; iP<nParticles; iP++)
    {
      partPDG[iP] = mPartPDG[iP];
      partName[iP] = mPartName[iP];
      partTitle[iP] = mPartTitle[iP];
      partMHistoMin[iP] = mPartMHistoMin[iP];
      partMHistoMax[iP] = mPartMHistoMax[iP];
      partMaxMult[iP] = mPartMaxMult[iP];
    }

    for(int iP=0; iP<nParticles; iP++)
    {
      AddCounter(Form("%s",partName[iP].Data()), Form("%-*s",14,partTitle[iP].Data()));
      AddCounter(Form("%s_prim",partName[iP].Data()), Form("%s Prim",partTitle[iP].Data()));
      AddCounter(Form("%s_sec",partName[iP].Data()), Form("%s Sec ",partTitle[iP].Data()));
    }

    for(int iP=0; iP<nParticles; iP++)
      fPdgToIndex[mPartPDG[iP]] = iP;
  }

  virtual ~KFPartEfficiencies(){};

  int GetParticleIndex(int pdg)
  {
    std::map<int, int>::iterator it;
    it=fPdgToIndex.find(pdg);
    if(it != fPdgToIndex.end()) return it->second;
    else return -1;
  }

  virtual void AddCounter(TString shortname, TString name){
    indices[shortname] = names.size();
    names.push_back(name);

    ratio_reco1.AddCounter();
    ratio_reco2.AddCounter();
    ratio_reco3.AddCounter();

    mc1.AddCounter();
    mc2.AddCounter();
    mc3.AddCounter();
    
    reco.AddCounter();

    ratio_ghost.AddCounter();
    ratio_bg.AddCounter();
    ratio_clone.AddCounter();
    ghost.AddCounter();
    bg.AddCounter();
    clone.AddCounter();
  };

  KFPartEfficiencies& operator+=(KFPartEfficiencies& a){
    mc1 += a.mc1; mc2 += a.mc2; mc3 += a.mc3; reco += a.reco;
    ghost += a.ghost; bg += a.bg; clone += a.clone;
    return *this;
  };
  
  void CalcEff(){
    ratio_reco1 = reco/mc1;
    ratio_reco2 = reco/mc2;
    ratio_reco3 = reco/mc3;

    KFMCCounter<int> allReco = reco + ghost + bg;
    ratio_ghost = ghost/allReco;
    ratio_bg  = bg/allReco;
    ratio_clone  = clone/allReco;
  };
  

  void Inc(bool isReco, int nClones, bool isMC1, bool isMC2, bool isMC3, TString name)
  {
    const int index = indices[name];
    
    if(isMC1) mc1.counters[index]++;
    if(isMC2) mc2.counters[index]++;
    if(isMC3) mc3.counters[index]++;
    
    if(isReco) reco.counters[index]++;
    if(nClones > 0)
      clone.counters[index] += nClones;
  };

  void IncReco(bool isGhost, bool isBg, TString name){
    const int index = indices[name];

    if (isGhost) ghost.     counters[index]++;
    if (isBg)    bg.counters[index]++;
  };

  void PrintEff(){
    std::cout.setf(std::ios::fixed);
    std::cout.setf(std::ios::showpoint);
    std::cout.precision(3);
    std::cout << "Particle        : "
         << "  Eff1 "
         <<" / "<< "  Eff2 "
         <<" / "<< "  Eff3 "
         <<" / "<< " Ghost "
         <<" / "<< "BackGr "
         <<" / "<< "N Ghost"
         <<" / "<< "N BackGr"
         <<" / "<< "N Reco "
         <<" / "<< "N Clone "
         <<" | "<< " N MC1 " 
         <<" | "<< " N MC2 " 
         <<" | "<< " N MC3 "  << std::endl;
    
    int NCounters = mc1.NCounters;
    for (int iC = 0; iC < NCounters; iC++){
        std::cout << names[iC]
             << "  : " << std::setw(6) << ratio_reco1.counters[iC]    
             << "  / " << std::setw(6) << ratio_reco2.counters[iC]
             << "  / " << std::setw(6) << ratio_reco3.counters[iC]
             << "  / " << std::setw(6) << ratio_ghost.counters[iC]  // particles w\o MCParticle
             << "  / " << std::setw(6) << ratio_bg.counters[iC]     // particles with incorrect MCParticle
             << "  / " << std::setw(6) << ghost.counters[iC]
             << "  / " << std::setw(7) << bg.counters[iC]
             << "  / " << std::setw(6) << reco.counters[iC]
             << "  / " << std::setw(7) << clone.counters[iC]
             << "  | " << std::setw(6) << mc1.counters[iC] 
             << "  | " << std::setw(6) << mc2.counters[iC]
             << "  | " << std::setw(6) << mc3.counters[iC]  << std::endl;
    }
  };

  friend std::fstream & operator<<(std::fstream &strm, KFPartEfficiencies &a) {

    strm << a.ratio_reco1;
    strm << a.ratio_reco2;
    strm << a.ratio_reco3;
    strm << a.mc1;
    strm << a.mc2;
    strm << a.mc3;
    strm << a.reco;
    strm << a.ratio_ghost;
    strm << a.ratio_bg;
    strm << a.ratio_clone;
    strm << a.ghost;
    strm << a.bg;
    strm << a.clone;

    return strm;
  }

  friend std::fstream & operator>>(std::fstream &strm, KFPartEfficiencies &a){

    strm >> a.ratio_reco1;
    strm >> a.ratio_reco2;
    strm >> a.ratio_reco3;
    strm >> a.mc1;
    strm >> a.mc2;
    strm >> a.mc3;
    strm >> a.reco;
    strm >> a.ratio_ghost;
    strm >> a.ratio_bg;
    strm >> a.ratio_clone;
    strm >> a.ghost;
    strm >> a.bg;
    strm >> a.clone;

    return strm;
  }

  void AddFromFile(TString fileName)
  {
    std::fstream file(fileName.Data(),fstream::in);
    file >> *this;
  }

  static const int nParticles = 74;
  int partPDG[nParticles];
  TString partName[nParticles];
  TString partTitle[nParticles];
  std::vector<std::vector<int> > partDaughterPdg;
  float partMHistoMin[nParticles];
  float partMHistoMax[nParticles];
  int partMaxMult[nParticles];

 private:
  std::vector<TString> names; // names counters indexed by index of counter
  std::map<TString, int> indices; // indices of counters indexed by a counter shortname

  std::map<int, int> fPdgToIndex;

  KFMCCounter<double> ratio_reco1;
  KFMCCounter<double> ratio_reco2;
  KFMCCounter<double> ratio_reco3;

  KFMCCounter<int> mc1;
  KFMCCounter<int> mc2;
  KFMCCounter<int> mc3;

  KFMCCounter<int> reco;

  KFMCCounter<double> ratio_ghost;
  KFMCCounter<double> ratio_bg;
  KFMCCounter<double> ratio_clone;

  KFMCCounter<int> ghost;
  KFMCCounter<int> bg; // background
  KFMCCounter<int> clone; // background
};

#endif
