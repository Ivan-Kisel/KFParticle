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

#ifndef KFParticleTopoReconstructor_H
#define KFParticleTopoReconstructor_H

/*
 * Interface class for use KFParticle Topology Reconstructor
 */


#include "KFParticlePVReconstructor.h"
#include <vector>
#include <string>

#include "KFPTrackVector.h"
#include "KFParticleSIMD.h"

#ifdef USE_TIMERS
#include "Stopwatch.h"
#endif

#ifdef KFPWITHTRACKER
class AliHLTTPCCAGBTracker;
#endif

class KFParticleTopoReconstructor{
 public:
  KFParticleTopoReconstructor():fKFParticlePVReconstructor(0),fTracks(0),fNThreads(1)
#ifdef USE_TIMERS
  ,fTime(0.)
#endif
  {
#ifdef USE_TIMERS
    for ( int i = 0; i < fNTimers; i++ ) fStatTime[i] = 0;
#endif
    fKFParticlePVReconstructor = new KFParticlePVReconstructor;
  }
  ~KFParticleTopoReconstructor();

#ifdef KFPWITHTRACKER
  void Init(AliHLTTPCCAGBTracker* tracker, vector<int>* pdg=0); // init array of particles
#endif
  void Init(vector<KFParticle> &particles, vector<int>* pdg=0);
  void Init(const KFPTrackVector *particles, const std::vector<KFParticle>& pv, const std::vector<short int> clusterPV);
  void Init(KFPTrackVector &tracks, vector<int>* pdg);
  
  void DeInit() { fTracks = 0; }
  
  void ReconstructPrimVertex(bool isHeavySystem = 1); // find primary vertex
  void SortTracks(); //sort tracks according to the pdg hypothesis and pv index
  void ReconstructParticles(); //find short-lived particles 

   /// Accessors
  int NPrimaryVertices() const { return fKFParticlePVReconstructor->NPrimaryVertices(); }
  KFParticle &GetPrimVertex(int iPV=0) const { return fKFParticlePVReconstructor->GetPrimVertex(iPV); }
  vector<short int>& GetPVTrackIndexArray(int iPV=0) const { return fKFParticlePVReconstructor->GetPVTrackIndexArray(iPV); }
  
  vector<KFParticle> const &GetParticles() const { return fParticles; }
  const KFPTrackVector* GetTracks() const { return fTracks; }
  const kfvector_float* GetChiPrim() const { return fChiToPrimVtx; }
  
  void AddPV(const KFVertex &pv, const vector<short int> &tracks) { fKFParticlePVReconstructor->AddPV(pv,tracks);}
  void AddPV(const KFVertex &pv) { fKFParticlePVReconstructor->AddPV(pv);}

  void SetBeamLine(KFParticle& p) { fKFParticlePVReconstructor->SetBeamLine(p); }

  void SetField(double b);
  //speed measurements
#ifdef USE_TIMERS
  void SetTime(double d) { fTime = d; }
  double Time() const { return fTime; }
  double StatTime( int iTimer ) const { return fStatTime[iTimer]; }
  int NTimers() const { return fNTimers; }
#endif

  void SaveInputParticles(const std::string prefix = "KFPData", bool onlySecondary = 0);
  void SetNThreads(short int n) { fNThreads=n; }
  
 private:
  KFParticleTopoReconstructor &operator=(KFParticleTopoReconstructor &);
  KFParticleTopoReconstructor(KFParticleTopoReconstructor &);

  void GetChiToPrimVertex(KFParticleSIMD* pv, const int nPV);

  KFParticlePVReconstructor* fKFParticlePVReconstructor;

  KFPTrackVector *fTracks;
  kfvector_float fChiToPrimVtx[2];
  vector<KFParticle> fParticles;
  vector<KFParticleSIMD> fPV;
    
  short int fNThreads;
  
  //speed measurements
#ifdef USE_TIMERS
  double fTime; //* total time
  static const int fNTimers = 4;
  double fStatTime[fNTimers]; //* timers
  Stopwatch timer;
#endif // USE_TIMERS

}; // class KFParticleTopoReconstructor


  
#endif // KFParticleTopoReconstructor_H
  
