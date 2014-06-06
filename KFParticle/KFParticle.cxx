//----------------------------------------------------------------------------
// Implementation of the KFParticle class
// .
// @author  S.Gorbunov, I.Kisel, I.Kulakov, M.Zyzak
// @version 1.0
// @since   13.05.07
// 
// Class to reconstruct and store the decayed particle parameters.
// The method is described in CBM-SOFT note 2007-003, 
// ``Reconstruction of decayed particles based on the Kalman filter'', 
// http://www.gsi.de/documents/DOC-2007-May-14-1.pdf
//
// This class is ALICE interface to general mathematics in KFParticleCore
// 
//  -= Copyright &copy ALICE HLT and CBM L1 Groups =-
//____________________________________________________________________________


#include "KFParticle.h"
#include "KFParticleDatabase.h"
#ifdef HomogeneousField
#include "KFPTrack.h"
#include "KFPVertex.h"
#endif

#ifdef NonhomogeneousField
#include "CbmKFTrack.h"
#endif

#ifdef HomogeneousField
float KFParticle::fgBz = -5.;  //* Bz compoment of the magnetic field
#endif

KFParticle::KFParticle( const KFParticle &d1, const KFParticle &d2, Bool_t gamma )
{
  if (!gamma) {
    KFParticle mother;
    mother+= d1;
    mother+= d2;
    *this = mother;
  } else
    ConstructGamma(d1, d2);
}

void KFParticle::Create( const float Param[], const float Cov[], Int_t Charge, float mass /*Int_t PID*/ )
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
  float C[21];
  for( int i=0; i<21; i++ ) C[i] = Cov[i];
  
//  TParticlePDG* particlePDG = TDatabasePDG::Instance()->GetParticle(PID);
//  float mass = (particlePDG) ? particlePDG->Mass() :0.13957;
  
  KFParticleBase::Initialize( Param, C, Charge, mass );
}

#ifdef NonhomogeneousField
KFParticle::KFParticle( CbmKFTrackInterface* Track, float *z0, Int_t  *qHypo, Int_t *PID)
{
  CbmKFTrack Tr( *Track );

  float *m = Tr.GetTrack();
  float *V = Tr.GetCovMatrix();

  //* Fit of vertex track slopes and momentum (a,b,qp) to xyz0 vertex

//  if( z0 ) Tr.Extrapolate( *z0 );

  float a = m[2], b = m[3], qp = m[4];
  
  //* convert the track to (x,y,px,py,pz,e,t/p) parameterization

  double c2 = 1./(1. + a*a + b*b);
  double pq = 1./qp;
  if(qHypo) pq *= *qHypo;
  double p2 = pq*pq;
  double pz = sqrt(p2*c2);
  double px = a*pz;
  double py = b*pz;
  double H[3] = { -px*c2, -py*c2, -pz*pq };

  float r[6];

  r[0] = m[0];
  r[1] = m[1];
  r[2] = m[5];
  r[3] = px;
  r[4] = py;
  r[5] = pz;

  double cxpz = H[0]*V[ 3] + H[1]*V[ 6] + H[2]*V[10];
  double cypz = H[0]*V[ 4] + H[1]*V[ 7] + H[2]*V[11];
  double capz = H[0]*V[ 5] + H[1]*V[ 8] + H[2]*V[12];
  double cbpz = H[0]*V[ 8] + H[1]*V[ 9] + H[2]*V[13];
//  double cqpz = H[0]*V[12] + H[1]*V[13] + H[2]*V[14];
  double cpzpz = H[0]*H[0]*V[5] +H[1]*H[1]*V[9] + H[2]*H[2]*V[14] 
    + 2*( H[0]*H[1]*V[8] +H[0]*H[2]*V[12] +H[1]*H[2]*V[13]);

  float C[21];

  C[ 0] = V[0];
  C[ 1] = V[1];
  C[ 2] = V[2];
  C[ 3] = 0;
  C[ 4] = 0;
  C[ 5] = 0;
  C[ 6] = V[3]*pz + a*cxpz; 
  C[ 7] = V[4]*pz + a*cypz; 
  C[ 8] = 0;
  C[ 9] = V[5]*pz*pz + 2*a*pz*capz + a*a*cpzpz;
  C[10] = V[6]*pz+b*cxpz; 
  C[11] = V[7]*pz+b*cypz; 
  C[12] = 0;
  C[13] = V[8]*pz*pz + a*pz*cbpz + b*pz*capz + a*b*cpzpz;
  C[14] = V[9]*pz*pz + 2*b*pz*cbpz + b*b*cpzpz;
  C[15] = cxpz; 
  C[16] = cypz; 
  C[17] = 0;
  C[18] = capz*pz + a*cpzpz;
  C[19] = cbpz*pz + b*cpzpz;
  C[20] = cpzpz;

  float Charge = (qp>0.) ?1 :( (qp<0) ?-1 :0);
  if(qHypo) Charge *= *qHypo;

  float mass = Tr.GetMass();

  if(PID)
    mass=KFParticleDatabase::Instance()->GetMass(*PID);

  Create(r,C,Charge,mass);
}

// KFParticle::KFParticle( CbmKFVertexInterface &vertex )
// {
//   // Constructor from CBM KF vertex
//   fP[0] = vertex.GetRefX();
//   fP[1] = vertex.GetRefY();
//   fP[2] = vertex.GetRefZ();
//   for( Int_t i=0; i<6; i++)
//     fC[i] = vertex.GetCovMatrix( )[i];
//   fChi2 = vertex.GetRefChi2();
//   fNDF = vertex.GetRefNDF();//2*vertex.GetNContributors() - 3;
//   fQ = 0;
//   fAtProductionVertex = 0;
//   fIsLinearized = 0;
//   fSFromDecay = 0;
// }

#endif

#ifdef HomogeneousField
KFParticle::KFParticle( const KFPTrack &track, const int PID )
{
  // Constructor from ALICE track, PID hypothesis should be provided

  track.XvYvZv(fP);
  track.PxPyPz(fP+3);
  fQ = track.Charge();
  track.GetCovarianceXYZPxPyPz( fC );

  float mass = KFParticleDatabase::Instance()->GetMass(PID);

  Create(fP,fC,fQ,mass);
  fChi2 = track.GetChi2();
  fNDF = track.GetNDF();
}

KFParticle::KFParticle( const KFPVertex &vertex )
{
  // Constructor from ALICE vertex

  vertex.GetXYZ( fP );
  vertex.GetCovarianceMatrix( fC );  
  fChi2 = vertex.GetChi2();
  fNDF = 2*vertex.GetNContributors() - 3;
  fQ = 0;
  fAtProductionVertex = 0;
  fIsLinearized = 0;
  fSFromDecay = 0;
}

void KFParticle::GetExternalTrackParam( const KFParticleBase &p, float &X, float &Alpha, float P[5] ) 
{
  // Conversion to AliExternalTrackParam parameterization

  float cosA = p.GetPx(), sinA = p.GetPy(); 
  float pt = sqrt(cosA*cosA + sinA*sinA);
  float pti = 0;
  if( pt<1.e-4 ){
    cosA = 1;
    sinA = 0;
  } else {
    pti = 1./pt;
    cosA*=pti;
    sinA*=pti;
  }
  Alpha = atan2(sinA,cosA);  
  X   = p.GetX()*cosA + p.GetY()*sinA;   
  P[0]= p.GetY()*cosA - p.GetX()*sinA;
  P[1]= p.GetZ();
  P[2]= 0;
  P[3]= p.GetPz()*pti;
  P[4]= p.GetQ()*pti;
}

#endif


Bool_t KFParticle::GetDistanceFromVertexXY( const float vtx[], const float Cv[], float &val, float &err ) const
{
  //* Calculate DCA distance from vertex (transverse impact parameter) in XY
  //* v = [xy], Cv=[Cxx,Cxy,Cyy ]-covariance matrix

  Bool_t ret = 0;
  
  float mP[8];
  float mC[36];
  
  Transport( GetDStoPoint(vtx), mP, mC );  

  float dx = mP[0] - vtx[0];
  float dy = mP[1] - vtx[1];
  float px = mP[3];
  float py = mP[4];
  float pt = sqrt(px*px + py*py);
  float ex=0, ey=0;
  if( pt<1.e-4 ){
    ret = 1;
    pt = 1.;
    val = 1.e4;
  } else{
    ex = px/pt;
    ey = py/pt;
    val = dy*ex - dx*ey;
  }

  float h0 = -ey;
  float h1 = ex;
  float h3 = (dy*ey + dx*ex)*ey/pt;
  float h4 = -(dy*ey + dx*ex)*ex/pt;
  
  err = 
    h0*(h0*GetCovariance(0,0) + h1*GetCovariance(0,1) + h3*GetCovariance(0,3) + h4*GetCovariance(0,4) ) +
    h1*(h0*GetCovariance(1,0) + h1*GetCovariance(1,1) + h3*GetCovariance(1,3) + h4*GetCovariance(1,4) ) +
    h3*(h0*GetCovariance(3,0) + h1*GetCovariance(3,1) + h3*GetCovariance(3,3) + h4*GetCovariance(3,4) ) +
    h4*(h0*GetCovariance(4,0) + h1*GetCovariance(4,1) + h3*GetCovariance(4,3) + h4*GetCovariance(4,4) );

  if( Cv ){
    err+= h0*(h0*Cv[0] + h1*Cv[1] ) + h1*(h0*Cv[1] + h1*Cv[2] ); 
  }

  err = sqrt(fabs(err));

  return ret;
}

Bool_t KFParticle::GetDistanceFromVertexXY( const float vtx[], float &val, float &err ) const
{
  return GetDistanceFromVertexXY( vtx, 0, val, err );
}


Bool_t KFParticle::GetDistanceFromVertexXY( const KFParticle &Vtx, float &val, float &err ) const 
{
  //* Calculate distance from vertex [cm] in XY-plane

  return GetDistanceFromVertexXY( Vtx.fP, Vtx.fC, val, err );
}

#ifdef HomogeneousField
Bool_t KFParticle::GetDistanceFromVertexXY( const KFPVertex &Vtx, float &val, float &err ) const 
{
  //* Calculate distance from vertex [cm] in XY-plane

  return GetDistanceFromVertexXY( KFParticle(Vtx), val, err );
}
#endif

float KFParticle::GetDistanceFromVertexXY( const float vtx[] ) const
{
  //* Calculate distance from vertex [cm] in XY-plane
  float val, err;
  GetDistanceFromVertexXY( vtx, 0, val, err );
  return val;
}

float KFParticle::GetDistanceFromVertexXY( const KFParticle &Vtx ) const 
{
  //* Calculate distance from vertex [cm] in XY-plane

  return GetDistanceFromVertexXY( Vtx.fP );
}

#ifdef HomogeneousField
float KFParticle::GetDistanceFromVertexXY( const KFPVertex &Vtx ) const 
{
  //* Calculate distance from vertex [cm] in XY-plane

  return GetDistanceFromVertexXY( KFParticle(Vtx).fP );
}
#endif

float KFParticle::GetDistanceFromParticleXY( const KFParticle &p ) const 
{
  //* Calculate distance to other particle [cm]

  float dS, dS1;
  GetDStoParticleXY( p, dS, dS1 );   
  float mP[8], mC[36], mP1[8], mC1[36];
  Transport( dS, mP, mC ); 
  p.Transport( dS1, mP1, mC1 ); 
  float dx = mP[0]-mP1[0]; 
  float dy = mP[1]-mP1[1]; 
  return sqrt(dx*dx+dy*dy);
}

float KFParticle::GetDeviationFromParticleXY( const KFParticle &p ) const 
{
  //* Calculate sqrt(Chi2/ndf) deviation from other particle

  float dS, dS1;
  GetDStoParticleXY( p, dS, dS1 );   
  float mP1[8], mC1[36];
  p.Transport( dS1, mP1, mC1 ); 

  float d[2]={ fP[0]-mP1[0], fP[1]-mP1[1] };

  float sigmaS = .1+10.*sqrt( (d[0]*d[0]+d[1]*d[1] )/
					(mP1[3]*mP1[3]+mP1[4]*mP1[4] )  );

  float h[2] = { mP1[3]*sigmaS, mP1[4]*sigmaS };       
  
  mC1[0] +=h[0]*h[0];
  mC1[1] +=h[1]*h[0]; 
  mC1[2] +=h[1]*h[1]; 

  return GetDeviationFromVertexXY( mP1, mC1 )*sqrt(2./1.);
}


float KFParticle::GetDeviationFromVertexXY( const float vtx[], const float Cv[] ) const 
{
  //* Calculate sqrt(Chi2/ndf) deviation from vertex
  //* v = [xyz], Cv=[Cxx,Cxy,Cyy,Cxz,Cyz,Czz]-covariance matrix

  float val, err;
  Bool_t problem = GetDistanceFromVertexXY( vtx, Cv, val, err );
  if( problem || err<1.e-20 ) return 1.e4;
  else return val/err;
}


float KFParticle::GetDeviationFromVertexXY( const KFParticle &Vtx ) const  
{
  //* Calculate sqrt(Chi2/ndf) deviation from vertex
  //* v = [xyz], Cv=[Cxx,Cxy,Cyy,Cxz,Cyz,Czz]-covariance matrix

  return GetDeviationFromVertexXY( Vtx.fP, Vtx.fC );
}

#ifdef HomogeneousField
float KFParticle::GetDeviationFromVertexXY( const KFPVertex &Vtx ) const 
{
  //* Calculate sqrt(Chi2/ndf) deviation from vertex
  //* v = [xyz], Cv=[Cxx,Cxy,Cyy,Cxz,Cyz,Czz]-covariance matrix

  KFParticle v(Vtx);
  return GetDeviationFromVertexXY( v.fP, v.fC );
}
#endif

float KFParticle::GetAngle  ( const KFParticle &p ) const 
{
  //* Calculate the opening angle between two particles

  float dS, dS1;
  GetDStoParticle( p, dS, dS1 );   
  float mP[8], mC[36], mP1[8], mC1[36];
  Transport( dS, mP, mC ); 
  p.Transport( dS1, mP1, mC1 ); 
  float n = sqrt( mP[3]*mP[3] + mP[4]*mP[4] + mP[5]*mP[5] );
  float n1= sqrt( mP1[3]*mP1[3] + mP1[4]*mP1[4] + mP1[5]*mP1[5] );
  n*=n1;
  float a = 0;
  if( n>1.e-8 ) a = ( mP[3]*mP1[3] + mP[4]*mP1[4] + mP[5]*mP1[5] )/n;
  if (fabs(a)<1.) a = acos(a);
  else a = (a>=0) ?0 :3.14;
  return a;
}

float KFParticle::GetAngleXY( const KFParticle &p ) const 
{
  //* Calculate the opening angle between two particles in XY plane

  float dS, dS1;
  GetDStoParticleXY( p, dS, dS1 );   
  float mP[8], mC[36], mP1[8], mC1[36];
  Transport( dS, mP, mC ); 
  p.Transport( dS1, mP1, mC1 ); 
  float n = sqrt( mP[3]*mP[3] + mP[4]*mP[4] );
  float n1= sqrt( mP1[3]*mP1[3] + mP1[4]*mP1[4] );
  n*=n1;
  float a = 0;
  if( n>1.e-8 ) a = ( mP[3]*mP1[3] + mP[4]*mP1[4] )/n;
  if (fabs(a)<1.) a = acos(a);
  else a = (a>=0) ?0 :3.14;
  return a;
}

float KFParticle::GetAngleRZ( const KFParticle &p ) const 
{
  //* Calculate the opening angle between two particles in RZ plane

  float dS, dS1;
  GetDStoParticle( p, dS, dS1 );   
  float mP[8], mC[36], mP1[8], mC1[36];
  Transport( dS, mP, mC ); 
  p.Transport( dS1, mP1, mC1 ); 
  float nr = sqrt( mP[3]*mP[3] + mP[4]*mP[4] );
  float n1r= sqrt( mP1[3]*mP1[3] + mP1[4]*mP1[4]  );
  float n = sqrt( nr*nr + mP[5]*mP[5] );
  float n1= sqrt( n1r*n1r + mP1[5]*mP1[5] );
  n*=n1;
  float a = 0;
  if( n>1.e-8 ) a = ( nr*n1r +mP[5]*mP1[5])/n; 
  if (fabs(a)<1.) a = acos(a);
  else a = (a>=0) ?0 :3.14;
  return a;
}


/*

#include "AliExternalTrackParam.h"

void KFParticle::GetDStoParticleALICE( const KFParticleBase &p, 
					   float &DS, float &DS1 ) 
  const
{ 
  DS = DS1 = 0;   
  float x1, a1, x2, a2;
  float par1[5], par2[5], cov[15];
  for(int i=0; i<15; i++) cov[i] = 0;
  cov[0] = cov[2] = cov[5] = cov[9] = cov[14] = .001;

  GetExternalTrackParam( *this, x1, a1, par1 );
  GetExternalTrackParam( p, x2, a2, par2 );

  AliExternalTrackParam t1(x1,a1, par1, cov);
  AliExternalTrackParam t2(x2,a2, par2, cov);

  float xe1=0, xe2=0;
  t1.GetDCA( &t2, -GetFieldAlice(), xe1, xe2 );
  t1.PropagateTo( xe1, -GetFieldAlice() );
  t2.PropagateTo( xe2, -GetFieldAlice() );

  float xyz1[3], xyz2[3];
  t1.GetXYZ( xyz1 );
  t2.GetXYZ( xyz2 );
  
  DS = GetDStoPoint( xyz1 );
  DS1 = p.GetDStoPoint( xyz2 );

  return;
}
*/

  // * Pseudo Proper Time of decay = (r*pt) / |pt| * M/|pt|
float KFParticle::GetPseudoProperDecayTime( const KFParticle &pV, const float& mass, float* timeErr2 ) const
{ // TODO optimize with respect to time and stability
  const float ipt2 = 1/( Px()*Px() + Py()*Py() );
  const float mipt2 = mass*ipt2;
  const float dx = X() - pV.X();
  const float dy = Y() - pV.Y();

  if ( timeErr2 ) {
      // -- calculate error = sigma(f(r)) = f'Cf'
      // r = {x,y,px,py,x_pV,y_pV}
      // df/dr = { px*m/pt^2,
      //     py*m/pt^2,
      //    ( x - x_pV )*m*(1/pt^2 - 2(px/pt^2)^2),
      //    ( y - y_pV )*m*(1/pt^2 - 2(py/pt^2)^2),
      //     -px*m/pt^2,
      //     -py*m/pt^2 }
    const float f0 = Px()*mipt2;
    const float f1 = Py()*mipt2;
    const float mipt2derivative = mipt2*(1-2*Px()*Px()*ipt2);
    const float f2 = dx*mipt2derivative;
    const float f3 = -dy*mipt2derivative;
    const float f4 = -f0;
    const float f5 = -f1;

    const float& mC00 =    GetCovariance(0,0);
    const float& mC10 =    GetCovariance(0,1);
    const float& mC11 =    GetCovariance(1,1);
    const float& mC20 =    GetCovariance(3,0);
    const float& mC21 =    GetCovariance(3,1);
    const float& mC22 =    GetCovariance(3,3);
    const float& mC30 =    GetCovariance(4,0);
    const float& mC31 =    GetCovariance(4,1);
    const float& mC32 =    GetCovariance(4,3);
    const float& mC33 =    GetCovariance(4,4);
    const float& mC44 = pV.GetCovariance(0,0);
    const float& mC54 = pV.GetCovariance(1,0);
    const float& mC55 = pV.GetCovariance(1,1);

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
