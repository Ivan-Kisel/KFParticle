//----------------------------------------------------------------------------
// Implementation of the KFParticle class
// .
// @author  M.Zyzak
// @version 1.0
// @since   09.02.16
// 
// 
//  -= Copyright &copy ALICE HLT and CBM L1 Groups =-
//____________________________________________________________________________

#ifndef KFPVector_H
#define KFPVector_H

#include <Vc/Vc>
#include <Vc/limits>
using ::Vc::float_v;

template <typename T>
class KFPVector
{
 public: 
  KFPVector(): fData(0), fSize(0) {};
  KFPVector(unsigned int n): fData(0), fSize(0) { resize(n); }
  KFPVector(unsigned int n, T defaultValue): fData(0), fSize(0) { resize(n, defaultValue); }
  ~KFPVector() { if(fData) _mm_free(fData); }

  KFPVector(const T& a)
  {
    if (this != &a)
    { 
      resize(a.size());
      for(int i=0; i<fSize; i++)
        fData[i] = a.fData[i];
    }
    return *this;
  }

  T& operator= (const T& a)
  {
    if (this != &a)
    { 
      resize(a.size());
      for(int i=0; i<fSize; i++)
        fData[i] = a.fData[i];
    }
    return *this;
  }
  
  T& operator[] (unsigned int i)             { return fData[i]; }
  const T& operator[] (unsigned int i) const { return fData[i]; }
  
  unsigned int size() const { return fSize; }
  
  void resize(unsigned int n)
  {
    if(fData) _mm_free(fData);
    fSize = n;
    if(fSize > 0)
      fData = (T*) _mm_malloc(fSize*sizeof(T), sizeof(Vc::float_v));
  }

  void resize(unsigned int n, T defaultValue)
  {
    if(fData) _mm_free(fData);
    fSize = n;
    if(fSize > 0)
      fData = (T*) _mm_malloc(fSize*sizeof(T), sizeof(Vc::float_v));
    memset(fData, defaultValue, fSize);
  }

  void *operator new(size_t size, void *ptr) { return ::operator new(size, ptr);}
  void *operator new[](size_t size, void *ptr) { return ::operator new(size, ptr);}
  void *operator new(size_t size) { return _mm_malloc(size, sizeof(Vc::float_v)); }
  void *operator new[](size_t size) { return _mm_malloc(size, sizeof(Vc::float_v)); }
  void operator delete(void *ptr, size_t) { _mm_free(ptr); }
  void operator delete[](void *ptr, size_t) { _mm_free(ptr); }
  
 private:
  T* fData;
  unsigned int fSize;
} __attribute__((aligned(sizeof(float_v))));

#endif
