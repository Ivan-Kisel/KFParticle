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

#include "KFPVertex.h"

KFPVertex::KFPVertex():fChi2(-1.f), fNContributors(0), fNDF(-1)  
{ 
  for(int iP=0; iP<3; iP++)
    fP[iP] = 0;
  for(int iC=0; iC<6; iC++)
    fC[iC] = 0;
}
