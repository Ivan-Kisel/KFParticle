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

#ifndef KFPartMatch_H
#define KFPartMatch_H

#include <vector>


/** @class KFPartMatch
 ** @brief A structure to store matching information between simulated Monte Carlo and reconstructed particles.
 ** @author  M.Zyzak, I.Kisel
 ** @date 05.02.2019
 ** @version 1.0
 **
 ** The class is used in both directions: to store links from Monte Carlo to reconstructed particles and vise versa.
 ** It contains two kind of links: 1) PDG hypothesis of the reconstructed and Monte Carlo particles should be the same;
 ** 2) PDG code differs. 
 **/

struct KFPartMatch // used for Reco to MC match as well as for MC to Reco
{
  KFPartMatch():ids(),idsMI() {}
  
  bool IsMatched() const { return ids.size() != 0 || idsMI.size() != 0; } ///< Returns true if at least one link exists independently of the PDG hypothesis.
  bool IsMatchedWithPdg() const { return ids.size() != 0; }               ///< Returns true is at least one link with the correct PDG exists.
  int  GetBestMatch() const { 
    if      (ids.size()   != 0) return ids[0];
    else if (idsMI.size() != 0) return idsMI[0];
    else return -1;
  } ///< Returns first link with correct PDG if exists, otherwise first link with incorrect PDG. If no link exists returns "-1".
  int  GetBestMatchWithPdg() const { 
    if      (ids.size()   != 0) return ids[0];
    else return -1;
  } ///< Returns first link with correct PDG if exists, otherwise returns "-1".
  std::vector<int> ids;   ///< Vector of links, PDG hypothesis of the reconstructed particle is required to be the same as the PDG code of the Monte Carlo particle.
  std::vector<int> idsMI; ///< Vector of links, PDG hypothesis of the reconstructed particle differs from the Monte Carlo particle.
};

#endif
