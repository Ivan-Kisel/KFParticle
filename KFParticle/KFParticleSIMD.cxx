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



#include "KFParticleSIMD.h"
#include "KFParticle.h"
#include "KFParticleDatabase.h"

#ifdef NonhomogeneousField
#include "CbmKFTrack.h"
#endif

#ifdef HomogeneousField
float_v KFParticleSIMD::fgBz = -5.f;  //* Bz compoment of the magnetic field
#endif

static const float_v Zero = 0.f;
static const float_v One = 1.f;

KFParticleSIMD::KFParticleSIMD( const KFParticleSIMD &d1, const KFParticleSIMD &d2, Bool_t gamma )
{
  if (!gamma) {
    KFParticleSIMD mother;
    mother+= d1;
    mother+= d2;
    *this = mother;
  } else
    ConstructGamma(d1, d2);
}

void KFParticleSIMD::Create( const float_v Param[], const float_v Cov[], float_v Charge, float_v mass /*Int_t PID*/ )
{
  // Constructor from "cartesian" track, PID hypothesis should be provided
  //
  // Param[6] = { X, Y, Z, Px, Py, Pz } - position and momentum
  // Cov [21] = lower-triangular part of the covariance matrix:
  //
  //                (  0  .  .  .  .  . )
  //                (  1  2  .  .  .  . )
  //  Cov. matrix = (  3  4  5  .  .  . ) - numbering of covariance elements in Cov[]
  //                (  6  7  8  9  .  . )
  //                ( 10 11 12 13 14  . )
  //                ( 15 16 17 18 19 20 )
  float_v C[21];
  for( int i=0; i<21; i++ ) C[i] = Cov[i];
  
//  TParticlePDG* particlePDG = TDatabasePDG::Instance()->GetParticle(PID);
//  float_v mass = (particlePDG) ? particlePDG->Mass() :0.13957;
  
  KFParticleBaseSIMD::Initialize( Param, C, Charge, mass );
}

KFParticleSIMD::KFParticleSIMD( const KFPTrack *track, Int_t PID )
{
  Double_t r[3];
  Double_t C[21];

  for(Int_t iPart = 0; iPart<float_vLen; iPart++)
  {
    track[iPart].XvYvZv(r);
    for(Int_t i=0; i<3; i++)
      fP[i][iPart] = r[i];
    track[iPart].PxPyPz(r);
    for(Int_t i=0; i<3; i++)
      fP[i+3][iPart] = r[i];
    fQ[iPart] = track[iPart].Charge();
    track[iPart].GetCovarianceXYZPxPyPz( C );
    for(Int_t i=0; i<21; i++)
      fC[i][iPart] = C[i];
  }

  float_v mass = KFParticleDatabase::Instance()->GetMass(PID);
  Create(fP,fC,fQ,mass);

  for(Int_t iPart = 0; iPart<float_vLen; iPart++)
  {
    fChi2[iPart] = track[iPart].GetChi2();
    fNDF[iPart] = track[iPart].GetNDF();
  }
}

KFParticleSIMD::KFParticleSIMD(KFPTrack &Track, Int_t *qHypo, const Int_t *pdg)
{
  Double_t r[3];
  Double_t C[21];

  Track.XvYvZv(r);
  for(Int_t i=0; i<3; i++)
    fP[i] = r[i];
  Track.PxPyPz(r);
  for(Int_t i=0; i<3; i++)
    fP[i+3] = r[i];
  fQ = Track.Charge();
  Track.GetCovarianceXYZPxPyPz( C );
  for(Int_t i=0; i<21; i++)
    fC[i] = C[i];

  float_v mass = KFParticleDatabase::Instance()->GetMass(*pdg);
  Create(fP,fC,fQ,mass);

  fChi2 = Track.GetChi2();
  fNDF = Track.GetNDF();
}

KFParticleSIMD::KFParticleSIMD(KFPTrackVector &track, int n, Int_t *qHypo, const Int_t *pdg)
{
  for(int i=0; i<6; i++)
    fP[i] = track.Parameter(i)[n];
  for(int i=0; i<21; i++)
    fC[i] = track.Covariance(i)[n];
#ifdef NonhomogeneousField
  for(int i=0; i<10; i++)
    fField.fField[i] = track.FieldCoefficient(i)[n];
#endif
  fQ = track.Q()[n];

  float_v mass = KFParticleDatabase::Instance()->GetMass(*pdg);
  Create(fP,fC,fQ,mass);
}

  
KFParticleSIMD::KFParticleSIMD(KFPTrack* Track[], int NTracks, Int_t *qHypo, const Int_t *pdg)
{
  Create(Track, NTracks, qHypo, pdg);
}

void KFParticleSIMD::Create(KFPTrack* Track[], int NTracks, Int_t *qHypo, const Int_t *pdg)
{
  Double_t r[3];
  Double_t C[21];

  for(Int_t iPart = 0; iPart<float_vLen; iPart++)
  {
    Int_t iEntry = (iPart < NTracks) ? iPart : 0; 
    Track[iEntry]->XvYvZv(r);
    for(Int_t i=0; i<3; i++)
      fP[i][iEntry] = r[i];
    Track[iEntry]->PxPyPz(r);
    for(Int_t i=0; i<3; i++)
      fP[i+3][iEntry] = r[i];
    fQ[iEntry] = Track[iEntry]->Charge();
    Track[iEntry]->GetCovarianceXYZPxPyPz( C );
    for(Int_t i=0; i<21; i++)
      fC[i][iEntry] = C[i];
  }

  float_v mass = KFParticleDatabase::Instance()->GetMass(*pdg);
  Create(fP,fC,fQ,mass);

  for(Int_t iPart = 0; iPart<float_vLen; iPart++)
  {
    Int_t iEntry = (iPart < NTracks) ? iPart : 0; 
    fChi2[iEntry] = Track[iEntry]->GetChi2();
    fNDF[iEntry] = Track[iEntry]->GetNDF();
  }
}

KFParticleSIMD::KFParticleSIMD(KFPTrackVector &track, uint_v& index, const int_v& pdg)
{
  Create(track, index, pdg);
}

void KFParticleSIMD::Create(KFPTrackVector &track, uint_v& index, const int_v& pdg)
{
  for(int i=0; i<6; i++)
    fP[i].gather(&(track.Parameter(i)[0]), index);
  for(int i=0; i<21; i++)
    fC[i].gather(&(track.Covariance(i)[0]), index);
#ifdef NonhomogeneousField
  for(int i=0; i<10; i++)
    fField.fField[i].gather(&(track.FieldCoefficient(i)[0]), index);
#endif
  
  //   fPDG.gather(&(track.PDG()[0]), index);
  fQ = track.Q()[index[0]];

  float_v mass = KFParticleDatabase::Instance()->GetMass(pdg);
  Create(fP,fC,fQ,mass);
}

KFParticleSIMD::KFParticleSIMD( const KFPVertex &vertex )
{
  // Constructor from ALICE vertex

  Double_t r[3];
  Double_t C[21];

  vertex.GetXYZ( r );
  for(Int_t i=0; i<3; i++)
    fP[i] = r[i];
  vertex.GetCovarianceMatrix( C );
  for(Int_t i=0; i<21; i++)
    fC[i] = C[i];
  fChi2 = vertex.GetChi2();
  fNDF = 2*vertex.GetNContributors() - 3;
  fQ = Zero;
  fAtProductionVertex = 0;
  fIsLinearized = 0;
  fSFromDecay = 0;
}

void KFParticleSIMD::SetOneEntry(int iEntry, KFParticleSIMD& part, int iEntryPart)
{
  for( int i = 0; i < 7; ++i )
    fP[i][iEntry] = part.Parameters()[i][iEntryPart];
  for( int i = 0; i < 36; ++i )
    fC[i][iEntry] = part.CovarianceMatrix()[i][iEntryPart];
  
  fQ[iEntry] = part.Q()[iEntryPart];
  fNDF[iEntry] = part.NDF()[iEntryPart];
  fChi2[iEntry] = part.Chi2()[iEntryPart];
  
//   fSFromDecay[iEntry] = part.fSFromDecay[iEntryPart];
//   SumDaughterMass[iEntry] = part.SumDaughterMass[iEntryPart];
//   fMassHypo[iEntry] = part.fMassHypo[iEntryPart];

  fId[iEntry] = part.Id()[iEntryPart];

  fPDG[iEntry] = part.GetPDG()[iEntryPart];
  
  if(iEntry==0)
    fDaughterIds.resize( part.NDaughters(), int_v(-1) );
  
  for(int iD=0; iD<part.NDaughters(); iD++)
    fDaughterIds[iD][iEntry] = part.fDaughterIds[iD][iEntryPart];

#ifdef NonhomogeneousField
  fField.SetOneEntry( iEntry, part.fField, iEntryPart ); //CHECKME
#endif
}

KFParticleSIMD::KFParticleSIMD(KFParticle* parts[], const int nPart)
{
  { // check
    bool ok = 1;

    const int nD = (parts[0])->NDaughters();
    for ( int ie = 1; ie < nPart; ie++ ) {
      const KFParticle &part = *(parts[ie]);
      ok &= part.NDaughters() == nD;
    }
//    assert(ok);
    if (!ok) {
      std::cout << " void CbmKFParticle_simd::Create(CbmKFParticle *parts[], int N) " << std::endl;
      exit(1);
    }
  }
  fDaughterIds.resize( (parts[0])->NDaughters(), int_v(-1) );

  for ( int iPart = 0; iPart < float_vLen; iPart++ ) {
    Int_t iEntry = (iPart < nPart) ? iPart : 0; 
    KFParticle &part = *(parts[iEntry]);

    fId[iEntry] = part.Id();
    for(int iD=0; iD<part.NDaughters(); iD++)
      fDaughterIds[iD][iEntry] = part.DaughterIds()[iD];

    fPDG[iEntry] = part.GetPDG();
    
    for( int i = 0; i < 8; ++i )
      fP[i][iEntry] = part.Parameters()[i];
    for( int i = 0; i < 36; ++i )
      fC[i][iEntry] = part.CovarianceMatrix()[i];

    fNDF[iEntry] = part.GetNDF();
    fChi2[iEntry] = part.GetChi2();
    fQ[iEntry] = part.GetQ();
    fAtProductionVertex = part.GetAtProductionVertex(); // CHECKME
#ifdef NonhomogeneousField
    fField.SetOneEntry( part.GetFieldCoeff(), iEntry);
#endif
  }
}

KFParticleSIMD::KFParticleSIMD( KFParticle &part)
{

 fId = part.Id();
 fNDF = part.GetNDF();
 fChi2 = part.GetChi2();
 fQ = part.GetQ();
 fPDG = part.GetPDG();
 fAtProductionVertex = part.GetAtProductionVertex();
 fIsVtxGuess = 0;
 fIsVtxErrGuess = 0;

  SetNDaughters(part.NDaughters());
  for( int i = 0; i < part.NDaughters(); ++i ) {
    fDaughterIds.push_back( part.DaughterIds()[i] );
  }
  
  for( int i = 0; i < 8; ++i )
    fP[i] = part.Parameters()[i];
  for( int i = 0; i < 36; ++i )
    fC[i] = part.CovarianceMatrix()[i];

#ifdef NonhomogeneousField
  fField = KFParticleFieldRegion(part.GetFieldCoeff());
#endif
}

float_m KFParticleSIMD::GetDistanceFromVertexXY( const float_v vtx[], const float_v Cv[], float_v &val, float_v &err ) const
{
  //* Calculate DCA distance from vertex (transverse impact parameter) in XY
  //* v = [xy], Cv=[Cxx,Cxy,Cyy ]-covariance matrix
  
  float_v mP[8];
  float_v mC[36];
  
  Transport( GetDStoPoint(vtx), mP, mC );  

  float_v dx = mP[0] - vtx[0];
  float_v dy = mP[1] - vtx[1];
  float_v px = mP[3];
  float_v py = mP[4];
  float_v pt = sqrt(px*px + py*py);
  float_v ex(Vc::Zero), ey(Vc::Zero);
  float_m mask = ( pt < float_v(1.e-4) );

  pt(mask) = One;
  ex(!mask) = (px/pt);
  ey(!mask) = (py/pt);
  val(mask) = float_v(1.e4);
  val(!mask)= dy*ex - dx*ey;

  float_v h0 = -ey;
  float_v h1 = ex;
  float_v h3 = (dy*ey + dx*ex)*ey/pt;
  float_v h4 = -(dy*ey + dx*ex)*ex/pt;
  
  err = 
    h0*(h0*GetCovariance(0,0) + h1*GetCovariance(0,1) + h3*GetCovariance(0,3) + h4*GetCovariance(0,4) ) +
    h1*(h0*GetCovariance(1,0) + h1*GetCovariance(1,1) + h3*GetCovariance(1,3) + h4*GetCovariance(1,4) ) +
    h3*(h0*GetCovariance(3,0) + h1*GetCovariance(3,1) + h3*GetCovariance(3,3) + h4*GetCovariance(3,4) ) +
    h4*(h0*GetCovariance(4,0) + h1*GetCovariance(4,1) + h3*GetCovariance(4,3) + h4*GetCovariance(4,4) );

  if( Cv ){
    err+= h0*(h0*Cv[0] + h1*Cv[1] ) + h1*(h0*Cv[1] + h1*Cv[2] ); 
  }

  err = sqrt(abs(err));

  return mask;
}

float_m KFParticleSIMD::GetDistanceFromVertexXY( const float_v vtx[], float_v &val, float_v &err ) const
{
  return GetDistanceFromVertexXY( vtx, 0, val, err );
}


float_m KFParticleSIMD::GetDistanceFromVertexXY( const KFParticleSIMD &Vtx, float_v &val, float_v &err ) const 
{
  //* Calculate distance from vertex [cm] in XY-plane

  return GetDistanceFromVertexXY( Vtx.fP, Vtx.fC, val, err );
}

#ifdef HomogeneousField
float_m KFParticleSIMD::GetDistanceFromVertexXY( const KFPVertex &Vtx, float_v &val, float_v &err ) const 
{
  //* Calculate distance from vertex [cm] in XY-plane

  return GetDistanceFromVertexXY( KFParticleSIMD(Vtx), val, err );
}
#endif

float_v KFParticleSIMD::GetDistanceFromVertexXY( const float_v vtx[] ) const
{
  //* Calculate distance from vertex [cm] in XY-plane
  float_v val, err;
  GetDistanceFromVertexXY( vtx, 0, val, err );
  return val;
}

float_v KFParticleSIMD::GetDistanceFromVertexXY( const KFParticleSIMD &Vtx ) const 
{
  //* Calculate distance from vertex [cm] in XY-plane

  return GetDistanceFromVertexXY( Vtx.fP );
}

#ifdef HomogeneousField
float_v KFParticleSIMD::GetDistanceFromVertexXY( const KFPVertex &Vtx ) const 
{
  //* Calculate distance from vertex [cm] in XY-plane

  return GetDistanceFromVertexXY( KFParticleSIMD(Vtx).fP );
}
#endif

float_v KFParticleSIMD::GetDistanceFromParticleXY( const KFParticleSIMD &p ) const 
{
  //* Calculate distance to other particle [cm]

  float_v dS, dS1;
  GetDStoParticleXY( p, dS, dS1 );   
  float_v mP[8], mC[36], mP1[8], mC1[36];
  Transport( dS, mP, mC ); 
  p.Transport( dS1, mP1, mC1 ); 
  float_v dx = mP[0]-mP1[0]; 
  float_v dy = mP[1]-mP1[1]; 
  return sqrt(dx*dx+dy*dy);
}

float_v KFParticleSIMD::GetDeviationFromParticleXY( const KFParticleSIMD &p ) const 
{
  //* Calculate sqrt(Chi2/ndf) deviation from other particle

  float_v dS, dS1;
  GetDStoParticleXY( p, dS, dS1 );   
  float_v mP1[8], mC1[36];
  p.Transport( dS1, mP1, mC1 ); 

  float_v d[2]={ fP[0]-mP1[0], fP[1]-mP1[1] };

  float_v sigmaS = .1f+10.f*sqrt( (d[0]*d[0]+d[1]*d[1] )/
					(mP1[3]*mP1[3]+mP1[4]*mP1[4] )  );

  float_v h[2] = { mP1[3]*sigmaS, mP1[4]*sigmaS };       
  
  mC1[0] +=h[0]*h[0];
  mC1[1] +=h[1]*h[0]; 
  mC1[2] +=h[1]*h[1]; 

  return GetDeviationFromVertexXY( mP1, mC1 )*float_v(sqrt(2.f));
}


float_v KFParticleSIMD::GetDeviationFromVertexXY( const float_v vtx[], const float_v Cv[] ) const 
{
  //* Calculate sqrt(Chi2/ndf) deviation from vertex
  //* v = [xyz], Cv=[Cxx,Cxy,Cyy,Cxz,Cyz,Czz]-covariance matrix

  float_v val, err;
  float_m problem = GetDistanceFromVertexXY( vtx, Cv, val, err );
  float_v ret = float_v(1.e4f);
  float_m mask = (problem || (err<float_v(1.e-20)));
  ret(!mask) = val/err;
  return ret;
}


float_v KFParticleSIMD::GetDeviationFromVertexXY( const KFParticleSIMD &Vtx ) const  
{
  //* Calculate sqrt(Chi2/ndf) deviation from vertex
  //* v = [xyz], Cv=[Cxx,Cxy,Cyy,Cxz,Cyz,Czz]-covariance matrix

  return GetDeviationFromVertexXY( Vtx.fP, Vtx.fC );
}

#ifdef HomogeneousField
float_v KFParticleSIMD::GetDeviationFromVertexXY( const KFPVertex &Vtx ) const 
{
  //* Calculate sqrt(Chi2/ndf) deviation from vertex
  //* v = [xyz], Cv=[Cxx,Cxy,Cyy,Cxz,Cyz,Czz]-covariance matrix

  KFParticleSIMD v(Vtx);
  return GetDeviationFromVertexXY( v.fP, v.fC );
}
#endif

float_v KFParticleSIMD::GetAngle  ( const KFParticleSIMD &p ) const 
{
  //* Calculate the opening angle between two particles

  float_v dS, dS1;
  GetDStoParticle( p, dS, dS1 );   
  float_v mP[8], mC[36], mP1[8], mC1[36];
  Transport( dS, mP, mC ); 
  p.Transport( dS1, mP1, mC1 ); 
  float_v n = sqrt( mP[3]*mP[3] + mP[4]*mP[4] + mP[5]*mP[5] );
  float_v n1= sqrt( mP1[3]*mP1[3] + mP1[4]*mP1[4] + mP1[5]*mP1[5] );
  n*=n1;
  float_v a = Zero;
  float_m mask = (n>(1.e-8f));
  a(mask) = ( mP[3]*mP1[3] + mP[4]*mP1[4] + mP[5]*mP1[5] )/n;
  mask = ( abs(a) < One);
  float_m aPos = (a>=Zero);
  a(mask) = KFPMath::ACos(a);
  a((!mask) && aPos) = Zero;
  a((!mask) && (!aPos)) = 3.1415926535f;
  return a;
}

float_v KFParticleSIMD::GetAngleXY( const KFParticleSIMD &p ) const 
{
  //* Calculate the opening angle between two particles in XY plane

  float_v dS, dS1;
  GetDStoParticleXY( p, dS, dS1 );   
  float_v mP[8], mC[36], mP1[8], mC1[36];
  Transport( dS, mP, mC ); 
  p.Transport( dS1, mP1, mC1 ); 
  float_v n = sqrt( mP[3]*mP[3] + mP[4]*mP[4] );
  float_v n1= sqrt( mP1[3]*mP1[3] + mP1[4]*mP1[4] );
  n*=n1;
  float_v a = Zero;
  float_m mask = (n>(1.e-8f));
  a = ( mP[3]*mP1[3] + mP[4]*mP1[4] )/n;
  a(!mask) = 0.f;
  mask = ( abs(a) < One);
  float_m aPos = (a>=Zero);
  a(mask) = KFPMath::ACos(a);
  a((!mask) && aPos) = Zero;
  a((!mask) && (!aPos)) = 3.1415926535f;
  return a;
}

float_v KFParticleSIMD::GetAngleRZ( const KFParticleSIMD &p ) const 
{
  //* Calculate the opening angle between two particles in RZ plane

  float_v dS, dS1;
  GetDStoParticle( p, dS, dS1 );   
  float_v mP[8], mC[36], mP1[8], mC1[36];
  Transport( dS, mP, mC ); 
  p.Transport( dS1, mP1, mC1 ); 
  float_v nr = sqrt( mP[3]*mP[3] + mP[4]*mP[4] );
  float_v n1r= sqrt( mP1[3]*mP1[3] + mP1[4]*mP1[4]  );
  float_v n = sqrt( nr*nr + mP[5]*mP[5] );
  float_v n1= sqrt( n1r*n1r + mP1[5]*mP1[5] );
  n*=n1;
  float_v a = Zero;
  float_m mask = (n>(1.e-8f));
  a(mask) = ( nr*n1r +mP[5]*mP1[5])/n;
  mask = ( abs(a) < One);
  float_m aPos = (a>=Zero);
  a(mask) = KFPMath::ACos(a);
  a((!mask) && aPos) = Zero;
  a((!mask) && (!aPos)) = 3.1415926535f;
  return a;
}


/*

#include "AliExternalTrackParam.h"

void KFParticleSIMD::GetDStoParticleALICE( const KFParticleSIMDBaseSIMD &p, 
					   float_v &DS, float_v &DS1 ) 
  const
{ 
  DS = DS1 = 0;   
  float_v x1, a1, x2, a2;
  float_v par1[5], par2[5], cov[15];
  for(int i=0; i<15; i++) cov[i] = 0;
  cov[0] = cov[2] = cov[5] = cov[9] = cov[14] = .001;

  GetExternalTrackParam( *this, x1, a1, par1 );
  GetExternalTrackParam( p, x2, a2, par2 );

  AliExternalTrackParam t1(x1,a1, par1, cov);
  AliExternalTrackParam t2(x2,a2, par2, cov);

  float_v xe1=0, xe2=0;
  t1.GetDCA( &t2, -GetFieldAlice(), xe1, xe2 );
  t1.PropagateTo( xe1, -GetFieldAlice() );
  t2.PropagateTo( xe2, -GetFieldAlice() );

  float_v xyz1[3], xyz2[3];
  t1.GetXYZ( xyz1 );
  t2.GetXYZ( xyz2 );
  
  DS = GetDStoPoint( xyz1 );
  DS1 = p.GetDStoPoint( xyz2 );

  return;
}
*/

  // * Pseudo Proper Time of decay = (r*pt) / |pt| * M/|pt|
float_v KFParticleSIMD::GetPseudoProperDecayTime( const KFParticleSIMD &pV, const float_v& mass, float_v* timeErr2 ) const
{ // TODO optimize with respect to time and stability
  const float_v ipt2 = 1/( Px()*Px() + Py()*Py() );
  const float_v mipt2 = mass*ipt2;
  const float_v dx = X() - pV.X();
  const float_v dy = Y() - pV.Y();

  if ( timeErr2 ) {
      // -- calculate error = sigma(f(r)) = f'Cf'
      // r = {x,y,px,py,x_pV,y_pV}
      // df/dr = { px*m/pt^2,
      //     py*m/pt^2,
      //    ( x - x_pV )*m*(1/pt^2 - 2(px/pt^2)^2),
      //    ( y - y_pV )*m*(1/pt^2 - 2(py/pt^2)^2),
      //     -px*m/pt^2,
      //     -py*m/pt^2 }
    const float_v f0 = Px()*mipt2;
    const float_v f1 = Py()*mipt2;
    const float_v mipt2derivative = mipt2*(1-2*Px()*Px()*ipt2);
    const float_v f2 = dx*mipt2derivative;
    const float_v f3 = -dy*mipt2derivative;
    const float_v f4 = -f0;
    const float_v f5 = -f1;

    const float_v& mC00 =    GetCovariance(0,0);
    const float_v& mC10 =    GetCovariance(0,1);
    const float_v& mC11 =    GetCovariance(1,1);
    const float_v& mC20 =    GetCovariance(3,0);
    const float_v& mC21 =    GetCovariance(3,1);
    const float_v& mC22 =    GetCovariance(3,3);
    const float_v& mC30 =    GetCovariance(4,0);
    const float_v& mC31 =    GetCovariance(4,1);
    const float_v& mC32 =    GetCovariance(4,3);
    const float_v& mC33 =    GetCovariance(4,4);
    const float_v& mC44 = pV.GetCovariance(0,0);
    const float_v& mC54 = pV.GetCovariance(1,0);
    const float_v& mC55 = pV.GetCovariance(1,1);

    *timeErr2 =
      f5*mC55*f5 +
      f5*mC54*f4 +
      f4*mC44*f4 +
      f3*mC33*f3 +
      f3*mC32*f2 +
      f3*mC31*f1 +
      f3*mC30*f0 +
      f2*mC22*f2 +
      f2*mC21*f1 +
      f2*mC20*f0 +
      f1*mC11*f1 +
      f1*mC10*f0 +
      f0*mC00*f0;
  }
  return ( dx*Px() + dy*Py() )*mipt2;
}

void KFParticleSIMD::GetKFParticle(KFParticle &Part, int iPart)
{
  Part.SetId(static_cast<int>(Id()[iPart]));

  Part.CleanDaughtersId();
  Part.SetNDaughters(DaughterIds().size());
  for( unsigned int i = 0; i < DaughterIds().size(); i++ )
    Part.AddDaughter(static_cast<int>(DaughterIds()[i][iPart]));

  Part.SetPDG( static_cast<int>(GetPDG()[iPart]) );

  for(int iP=0; iP<8; iP++)
    Part.Parameters()[iP] = Parameters()[iP][iPart];
  for(int iC=0; iC<36; iC++)
    Part.CovarianceMatrix()[iC] = CovarianceMatrix()[iC][iPart];

  Part.NDF() = static_cast<int>(GetNDF()[iPart]);
  Part.Chi2() = GetChi2()[iPart];
  Part.Q()    = GetQ()[iPart];
  Part.SetAtProductionVertex(fAtProductionVertex);
#ifdef NonhomogeneousField
  for(int iF=0; iF<10; iF++)
    Part.SetFieldCoeff(fField.fField[iF][iPart], iF);
#endif
}

void KFParticleSIMD::GetKFParticle(KFParticle* Part, int nPart)
{
  for(int i=0; i<nPart; i++)
    GetKFParticle(Part[i],i);
}
