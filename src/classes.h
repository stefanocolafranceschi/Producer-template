#ifndef DataFormats_SiPixelDigis_SiPixelDigisSoA_h
#define DataFormats_SiPixelDigis_SiPixelDigisSoA_h

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiPixelDigi/interface/SiPixelDigisSoA.h"  // Ensure this path is correct

// Forward declaration for SiPixelDigisLayout
template <int32_t size, bool aligned> class SiPixelDigisLayout;

typedef edmNew::DetSetVector<SiPixelDigisSoA> SiPixelDigisVector;

#endif  // DataFormats_SiPixelDigis_SiPixelDigisSoA_h
