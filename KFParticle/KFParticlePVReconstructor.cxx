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

#include "KFParticlePVReconstructor.h"
#include "KFPTrackVector.h"
#include "KFParticle.h"

#include <iostream>

void KFParticlePVReconstructor::Init(KFPTrackVector *tracks, int nParticles)
{
    // copy tracks in particles.
  fNParticles = nParticles;
  fParticles.resize(fNParticles);
  fWeight.resize(fNParticles);
  
  KFPTrack track;
  for ( int iTr = 0; iTr < fNParticles; iTr++ ) {
    // {   // check cov matrix TODO
    //   bool ok = true;  
    //   int k = 0;
    //   for (int i = 0; i < 6; i++) {
    //     for (int j = 0; j <= i; j++, k++) {
    //       ok &= finite( KFPCov[k] );
    //     }
    //     ok &= ( KFPCov[k-1] > 0 );
    //   }
    //   if (!ok) continue;
    // }
    tracks->GetTrack(track,iTr);
    fParticles[iTr] = KFParticle( track, 211 ); // pi+ // TODO include PDG in KFPTrack?
    fParticles[iTr].AddDaughterId(track.Id());
  }
  
  fPrimVertices.clear();
  fClusters.clear();
  
  float pvEstimation[3] = {0.};

  float parTmp[8] = {0.};
  float covTmp[36] = {0.};
  
  for(int iIter=0; iIter<3; iIter++)
  {
    float C[3] = {100.,100.,100.};
    for ( int iTr = 0; iTr < fNParticles; iTr++ ) {

      fParticles[iTr].Transport( fParticles[iTr].GetDStoPoint(pvEstimation), parTmp, covTmp);

      float r2 = fParticles[iTr].X()*fParticles[iTr].X() + fParticles[iTr].Y()*fParticles[iTr].Y();
      if(r2 > 100 ) continue;  
      
      const float V[3] = {covTmp[0], covTmp[2], covTmp[5]}; 
      for(int iComp=0; iComp<3; iComp++)
      {
        float K = C[iComp]/(C[iComp]+V[iComp]);
        if(C[iComp] > 16*V[iComp])
          K = 1.f - V[iComp]/C[iComp];
        const float dzeta = parTmp[iComp]-pvEstimation[iComp];
        pvEstimation[iComp] += K*dzeta;
        C[iComp] -= K*C[iComp];
      }
    }
  }

  {
    int nPrimCand = fParticles.size();

    const KFParticle **pParticles = new const KFParticle*[nPrimCand+1]; // tmp array
    bool *vFlags = new bool[nPrimCand];  // flags returned by the vertex finder

    KFVertex primVtx;
    primVtx.SetVtxGuess(pvEstimation[0], pvEstimation[1], pvEstimation[2]);
    primVtx.X() = pvEstimation[0];
    primVtx.Y() = pvEstimation[1];
    primVtx.Z() = pvEstimation[2];
    
    nPrimCand = 0;
    for ( unsigned int iP = 0; iP < fParticles.size(); iP++ )
    {
      const KFParticle &p = fParticles[iP];
      float chi = p.GetDeviationFromVertex( primVtx );      
      if( chi >= 3.f )
        continue;
      
      pParticles[nPrimCand] = &fParticles[iP];
      vFlags[nPrimCand] = 1;
      nPrimCand++;
    }
  
    primVtx.SetConstructMethod(0);
    primVtx.ConstructPrimaryVertex( pParticles, nPrimCand, vFlags, 3.f );
      
    delete [] pParticles;
    delete [] vFlags;
    pvEstimation[0] = primVtx.GetX();
    pvEstimation[1] = primVtx.GetY();
    pvEstimation[2] = primVtx.GetZ();
  }

  for(int iP=0; iP<fNParticles; iP++)
  {
    fParticles[iP].TransportToPoint(pvEstimation);

    fWeight[iP] = 1.f/sqrt(fParticles[iP].CovarianceMatrix()[0]
      + fParticles[iP].CovarianceMatrix()[2] + fParticles[iP].CovarianceMatrix()[5]);
  }

} // void KFParticlePVReconstructor::Init

void KFParticlePVReconstructor::FindPrimaryClusters( int cutNDF )
{
  // The function finds a set of clusters of tracks.
  // Tracks are assumed to be transported to the beam line.
  // If a beam line is set - it will be used for a reconstruction but will not be added to a daughter cluster

  if( IsBeamLine() )
    cutNDF += 2;
  
  vector<unsigned short int> notUsedTracks(fNParticles);
  vector<unsigned short int> *notUsedTracksPtr = &notUsedTracks;
  int nNotUsedTracks = fNParticles;

  vector<unsigned short int> notUsedTracksNew(fNParticles);
  vector<unsigned short int> *notUsedTracksNewPtr = &notUsedTracksNew;
  int nNotUsedTracksNew = 0;

  for(int iTr=0; iTr<fNParticles; iTr++)
    (*notUsedTracksPtr)[iTr] = iTr;
    
  while(nNotUsedTracks>0)
  {
    short int bestTrack = 0;
    float bestWeight = -1.f;

    for(unsigned short int iTr = 0; iTr < nNotUsedTracks; iTr++)
    {
      unsigned short int &curTrack = (*notUsedTracksPtr)[iTr];
      
      if (fWeight[curTrack] > bestWeight)
      {
        bestWeight = fWeight[curTrack];
        bestTrack = curTrack;
      }
    }

    KFParticleCluster cluster;
    cluster.fTracks.reserve(nNotUsedTracks);

    const float *rBest = fParticles[bestTrack].Parameters();
    const float *covBest = fParticles[bestTrack].CovarianceMatrix();

    float rVertex[3] = {0.f};
    float covVertex[6] = {0.f};
    float weightVertex = 0.f;

    for(unsigned short int iTr = 0; iTr < nNotUsedTracks; iTr++)
    {
      unsigned short int &curTrack = (*notUsedTracksPtr)[iTr];

      float dr[3] = {rBest[0] - fParticles[curTrack].X(),
                     rBest[1] - fParticles[curTrack].Y(),
                     rBest[2] - fParticles[curTrack].Z() };

      if(fabs(dr[2]) > 5.)
      {
        (*notUsedTracksNewPtr)[nNotUsedTracksNew] = curTrack;
        nNotUsedTracksNew++;
        continue;
      }
      float cov[6] = {covBest[0] + fParticles[curTrack].CovarianceMatrix()[0],
                      covBest[1] + fParticles[curTrack].CovarianceMatrix()[1],
                      covBest[2] + fParticles[curTrack].CovarianceMatrix()[2],
                      covBest[3] + fParticles[curTrack].CovarianceMatrix()[3],
                      covBest[4] + fParticles[curTrack].CovarianceMatrix()[4],
                      covBest[5] + fParticles[curTrack].CovarianceMatrix()[5] };

      float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      float drError2 = dr[0]* (cov[0]* dr[0] + cov[1]* dr[1] + cov[3]* dr[2]) +
                       dr[1]* (cov[1]* dr[0] + cov[2]* dr[1] + cov[4]* dr[2]) + 
                       dr[2]* (cov[3]* dr[0] + cov[4]* dr[1] + cov[5]* dr[2]);
      if(!(dr2==dr2) || !(drError2==drError2))
        continue;
  
      if( dr2*dr2 <= (25 * drError2) )
      {
        for(int iP=0; iP<3; iP++)
          rVertex[iP] += fWeight[curTrack] * fParticles[curTrack].Parameters()[iP];
        
        float weight2 = fWeight[curTrack] * fWeight[curTrack];
        for(int iC=0; iC<6; iC++)
          covVertex[iC] += weight2 *fParticles[curTrack].CovarianceMatrix()[iC];

        weightVertex += fWeight[curTrack];
        cluster.fTracks.push_back(curTrack);
      }
      else
      {
        (*notUsedTracksNewPtr)[nNotUsedTracksNew] = curTrack;
        nNotUsedTracksNew++;
      }
    }
    
    vector<unsigned short int> *notUsedTracksPtrSave = notUsedTracksPtr;
    notUsedTracksPtr = notUsedTracksNewPtr;
    notUsedTracksNewPtr = notUsedTracksPtrSave;
    
    nNotUsedTracks = nNotUsedTracksNew;
    nNotUsedTracksNew = 0;

    if(cluster.fTracks.size() < 2) continue;
    if((cluster.fTracks.size() < 3) && (fNParticles>3)) continue;

    for(int iP=0; iP<3; iP++)
      cluster.fP[iP] = rVertex[iP]/weightVertex;

    for(int iC=0; iC<6; iC++)
      cluster.fC[iC] = covVertex[iC]/(weightVertex*weightVertex);


    int nPrimCand = cluster.fTracks.size(); // 1 is reserved for a beam line
    int nTracks = cluster.fTracks.size();
    
    const KFParticle **pParticles = new const KFParticle*[nPrimCand+1]; // tmp array

    for ( int iP = 0; iP < nPrimCand; iP++ )
      pParticles[iP] = &fParticles[ cluster.fTracks[iP] ];

    if( IsBeamLine() )
    {
      pParticles[nPrimCand] = &fBeamLine;
      nPrimCand++;
    }
      
      
      // find prim vertex
    KFVertex primVtx;

    if(fNParticles>1)
    {
      // construct PV candidate from a cluster
      bool *vFlags = new bool[nPrimCand];  // flags returned by the vertex finder
      for(int iFl=0; iFl<nPrimCand; iFl++)
        vFlags[iFl] = true;
      primVtx.SetVtxGuess(cluster.fP[0], cluster.fP[1], cluster.fP[2]);
      primVtx.SetConstructMethod(0);
      primVtx.ConstructPrimaryVertex( pParticles, nPrimCand, vFlags, 3.f );

      // clean cluster
      vector<short int> clearClusterInd;
      clearClusterInd.reserve(cluster.fTracks.size());
      for ( int iP = 0; iP < nTracks; iP++ ){
        if(cluster.fTracks[iP] == bestTrack) {
          clearClusterInd.push_back(cluster.fTracks[iP]);
          continue;
        }
        if(vFlags[iP])
          clearClusterInd.push_back(cluster.fTracks[iP]);
        else
        {
          (*notUsedTracksPtr)[nNotUsedTracks] = cluster.fTracks[iP];
          nNotUsedTracks++;
        }
      }
      cluster.fTracks = clearClusterInd;

      // save PV
      if( primVtx.GetNDF() >= cutNDF ) //at least 2 particles
      {
//         std::cout << primVtx.X() << " " << primVtx.Y() << " " << primVtx.Z() << " " << cluster.fTracks.size() << std::endl;
        fPrimVertices.push_back(primVtx);
        fClusters.push_back(cluster);
      }
      
      if(vFlags) delete [] vFlags;
    }
    if(pParticles) delete [] pParticles;
  }
//   if(fClusters.size()>1)
//   {
//     float dx = fClusters[0].fP[0] - fClusters[1].fP[0];
//     float dy = fClusters[0].fP[1] - fClusters[1].fP[1];
//     float dz = fClusters[0].fP[2] - fClusters[1].fP[2];
// 
//     float dr[3] = {dx, dy, dz};
//     float cov[6] = {fClusters[0].fC[0] + fClusters[1].fC[0],
//                     fClusters[0].fC[1] + fClusters[1].fC[1],
//                     fClusters[0].fC[2] + fClusters[1].fC[2],
//                     fClusters[0].fC[3] + fClusters[1].fC[3],
//                     fClusters[0].fC[4] + fClusters[1].fC[4],
//                     fClusters[0].fC[5] + fClusters[1].fC[5] };
//       float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
//       float drError2 = dr[0]* (cov[0]* dr[0] + cov[1]* dr[1] + cov[3]* dr[2]) +
//                        dr[1]* (cov[1]* dr[0] + cov[2]* dr[1] + cov[4]* dr[2]) + 
//                        dr[2]* (cov[3]* dr[0] + cov[4]* dr[1] + cov[5]* dr[2]);
//     drError2 /= dr2;
//   }

//   static int nVert[10]={0.};
//   if(fClusters.size()<10)
//     nVert[fClusters.size()]++;
//   std::cout << "N Vert     ";
//   for(int i=0; i<10; i++)
//     std::cout << i << ": " << nVert[i] << "     ";
//   std::cout << std::endl;
// 
//   static int nPart[10]={0.};
//   if(fNParticles<10)
//     nPart[fNParticles]++;
//   std::cout << "N Part     ";
//   for(int i=0; i<10; i++)
//     std::cout << i << ": " << nPart[i] << "     ";
//   std::cout << std::endl;
}

void KFParticlePVReconstructor::ReconstructPrimVertex()
{

  FindPrimaryClusters();

  if ( fPrimVertices.size() == 0 ) { // fill prim vertex by dummy values

    float X=0,Y=0,Z=0;

    KFPVertex primVtx_tmp;
    primVtx_tmp.SetXYZ(X, Y, Z);
    primVtx_tmp.SetCovarianceMatrix( 0, 0, 0, 0, 0, 0 );
    primVtx_tmp.SetNContributors(0);
    primVtx_tmp.SetChi2(-100);

    fPrimVertices.push_back(primVtx_tmp);
    
    KFParticleCluster cluster;
    fClusters.push_back(cluster);
  }
} // void KFParticlePVReconstructor::Run()

void KFParticlePVReconstructor::AddPV(const KFVertex &pv, const vector<short int> &tracks)
{
  fPrimVertices.push_back(pv);
  KFParticleCluster cluster;
  cluster.fTracks = tracks;
  fClusters.push_back(cluster);
}

void KFParticlePVReconstructor::AddPV(const KFVertex &pv)
{
  fPrimVertices.push_back(pv);
  KFParticleCluster cluster;
  fClusters.push_back(cluster);
}