#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiPixelDigi/interface/SiPixelDigisSoA.h"  // Ensure correct header inclusion

// Explicit template instantiations
template class edmNew::DetSetVector<SiPixelDigisSoA>;
template class edm::Wrapper<edmNew::DetSetVector<SiPixelDigisSoA>>;  // Wrapper instantiation

namespace {
  struct dictionary {
    SiPixelDigisSoA dummy1;  // Instantiate SiPixelDigisSoA
    edmNew::DetSetVector<SiPixelDigisSoA> dummy2;  // Instantiate DetSetVector<SiPixelDigisSoA>
    edm::Wrapper<edmNew::DetSetVector<SiPixelDigisSoA>> dummy3;  // Instantiate Wrapper
  };
}
