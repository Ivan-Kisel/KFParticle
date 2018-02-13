//----------------------------------------------------------------------------
// Structures to collect histograms online
// .
// @author  M.Zyzak
// @version 1.0
// @since   20.08.13
// 
// 
//  -= Copyright &copy ALICE HLT and CBM L1 Groups =-
//____________________________________________________________________________

#ifndef KFPHISTOGRAMSET
#define KFPHISTOGRAMSET

#include "KFPHistogram1D.h"
#include "KFParticle.h"

class KFPHistogramSet
{
 public:
  KFPHistogramSet(int iPart=0);
  ~KFPHistogramSet() {};
  
  void Fill(const KFParticle& particle);
  
  inline int GetNHisto1D() const { return NHisto1D; }
  inline int DataSize() const 
  {
    int dataSize = 0;
    for(int i=0; i<NHisto1D; i++)
      dataSize += fKFPHistogram1D[i].DataSize();
    return dataSize;
  }
  KFPHistogram1D GetHistogram1D(int iHistogram) const { return fKFPHistogram1D[iHistogram]; }
  
  inline void SetHisto1DBinContent(int iHisto, int iBin, int value)   { fKFPHistogram1D[iHisto].SetBinContent(iBin,value); }
  
  inline void operator += ( const KFPHistogramSet &h )
  {
    for(int i=0; i<NHisto1D; i++)
      fKFPHistogram1D[i] += h.fKFPHistogram1D[i];
  }
  
  void SetHistogramMemory(int* pointer);
  
 private:
  static const int NHisto1D = 17;
  KFPHistogram1D fKFPHistogram1D[NHisto1D];
};

#endif
