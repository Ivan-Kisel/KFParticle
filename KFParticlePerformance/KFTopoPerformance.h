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

#ifdef DO_TPCCATRACKER_EFF_PERFORMANCE

#ifndef KFTOPOPERFORMANCE_H
#define KFTOPOPERFORMANCE_H


#include "KFParticlePerformanceBase.h"

#include "KFMCVertex.h"
#include "KFMCTrack.h"
#include <cstdio>
#include <map>

#include "KFPartMatch.h"
#include "KFMCParticle.h"

class AliHLTTPCCAGBTracker;


class TObject;
class TParticle;
class KFParticleTopoReconstructor;
class KFPHistogram;
class TDirectory;
class TH1D;
class TH2D;
class TH3D;
class TProfile;

class TFile;

/**
 * @class KFTopoPerformance. Don't use w\o GlobalPerformance
 */
class KFTopoPerformance: public KFParticlePerformanceBase
{
 public:
  
  KFTopoPerformance();
  virtual ~KFTopoPerformance();
#ifdef KFPWITHTRACKER
  virtual void SetNewEvent(
    const AliHLTTPCCAGBTracker * const Tracker,
    AliHLTResizableArray<AliHLTTPCCAHitLabel> *hitLabels,
    AliHLTResizableArray<AliHLTTPCCAMCTrack> *mcTracks,
    AliHLTResizableArray<AliHLTTPCCALocalMCPoint> *localMCPoints);
#endif  
  void SetTopoReconstructor( const KFParticleTopoReconstructor * const TopoReconstructor ); // use together with SetNewEvent !!!
  const KFParticleTopoReconstructor * GetTopoReconstructor() const { return fTopoReconstructor; }
    
    /// Efficiency
    // Check if MC track is reconstructable. Calculate set of MC track. Etc.
  virtual void CheckMCTracks(); // fill mcData.
    // Find reco-MCTracks correspondence
  virtual void MatchTracks();   // fill recoData.
    // Calculate efficiencies
  virtual void EfficiencyPerformance(){}; // current don't use eff

  virtual void PrintEfficiencyStatistic(){}; // current don't use eff
  virtual void PrintEfficiency()         {};
  
    /// Histograms
    //     virtual void CreateHistos(string histoDir);
  virtual void FillHistos();
  void FillHistos(const KFPHistogram* histograms);
  void FillMCHistos();

  void AddV0Histos();
  
  void SetTrackMatch(const std::vector<int>& trackMatch) { fTrackMatch = trackMatch;}
  void SetMCTracks(const std::vector<KFMCTrack>& mcTracks)
  { 
    
    vMCTracks = mcTracks;
  }
  
  const KFPartEfficiencies GetEfficiency() const { return fParteff; }
  void SetPrintEffFrequency(int n) { fPrintEffFrequency = n;}

  const std::vector<KFMCVertex> GetPrimVertices() { return fPrimVertices; }
  const std::vector<KFMCParticle>& MCParticles() { return vMCParticles; }
  const std::vector<KFPartMatch>& ParticlesMatch() { return RtoMCParticleId; }
  const std::vector<KFPartMatch>& GetMCtoRPVId() { return MCtoRPVId; }
  const std::vector<KFPartMatch>& GetRtoMCPVId() { return RtoMCPVId; }
  const KFMCTrack& GetMCTrack(const int iRecoTrack)
  { 
    int iMCTrack = 0;
    if(RtoMCParticleId[iRecoTrack].IsMatched())
      iMCTrack = RtoMCParticleId[iRecoTrack].GetBestMatch();
    return vMCTracks[iMCTrack]; 
  }
  
  void SetCentralityBin(const int iBin) { fCentralityBin = iBin; }
  void SetCentralityWeight(const float weight) { fCentralityWeight = weight; }
  
 private:

  const KFTopoPerformance& operator = (const KFTopoPerformance&);
  KFTopoPerformance(const KFTopoPerformance&);
   
  void GetMCParticles();
  void MatchParticles();
  void MatchPV();
  void CalculateEfficiency();
  void CalculatePVEfficiency();
  void FindReconstructableMCParticles();
  void CheckMCParticleIsReconstructable(KFMCParticle &part);
  void FindReconstructableMCVertices();
  void FillParticleParameters(KFParticle& TempPart,
                              int iParticle,
                              int iP,
                              int iPV,
                              TH1F* histoParameters[4][KFPartEfficiencies::nParticles][nHistoPartParam],
                              TH2F* histoParameters2D[4][KFPartEfficiencies::nParticles][nHistoPartParam2D],
                              TH3F* histoParameters3D[1][KFPartEfficiencies::nParticles][nHistoPartParam3D],
                              TH1F* histoFit[KFPartEfficiencies::nParticles][nFitQA] = 0,
                              TH1F* histoFitDaughtersQA[KFPartEfficiencies::nParticles][nFitQA] = 0,
                              TH1F* histoDSToParticleQA[KFPartEfficiencies::nParticles][nDSToParticleQA] = 0,
                              std::vector<int>* multiplicities = 0);
  
  const KFParticleTopoReconstructor *fTopoReconstructor;

  std::vector<KFMCVertex> fPrimVertices; // primary vertex positions (currently only one vertex is implemented)
  std::vector<int> fMCTrackToMCPVMatch; // match between MC tracks and PV 
  std::vector<double> fPVPurity;
  std::vector<double> fPVTracksRate[4]; //0 - ghost tracks, 1 - from trigger PV, 2 - from pileup, 3 - from physics bg
  std::vector<int> fNCorrectPVTracks;

  std::vector<int> fTrackMatch;
  std::vector<KFMCTrack> vMCTracks;  // MC particles
  std::vector<KFMCParticle> vMCParticles;  // MC particles
  std::vector<int> fNeutralIndex;
  
  std::vector<KFPartMatch> MCtoRParticleId; // array for match MC and reco particles
  std::vector<KFPartMatch> RtoMCParticleId;

  std::vector<KFPartMatch> MCtoRPVId; // array for match MC and reco PV
  std::vector<KFPartMatch> RtoMCPVId;
  
  int fPrintEffFrequency;
  
  KFPartEfficiencies fPartInfo;
  
  int fCentralityBin;
  float fCentralityWeight;
};

#endif
#endif //DO_TPCCATRACKER_EFF_PERFORMANCE
