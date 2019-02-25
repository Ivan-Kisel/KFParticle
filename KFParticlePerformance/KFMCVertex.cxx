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

#include "KFMCVertex.h"

KFMCVertex::KFMCVertex():fDaughterTracks(0),fIsReconstructable(0),fIsMCReconstructable(0),fIsReconstructed(0),fNReconstructedDaughters(0),fIsTriggerPV(0)
{
  for( int i = 0; i < 3; i++) fPar[i] = 0;
}

std::ostream& operator<<(std::ostream& out, const KFMCVertex &a)
{
  /** Operator to print coordinates of the MC vertex "a".
   ** \param[in] out - stream, where coordinates will be printed
   ** \param[in] a - vertex to be printed
   **/
  for (int i = 0; i < 3; i++) out << a.fPar[i] << std::endl;
  return out;
}

std::istream& operator>>(std::istream& in, KFMCVertex &a)
{
  /** Operator to read coordinates of the vertex from the input stream. 
   ** \param[in] in - input stream
   ** \param[in] a - vertex, where the coordinates will be read in
   **/
  for (int i = 0; i < 3; i++) in >> a.fPar[i];
  return in;
}

