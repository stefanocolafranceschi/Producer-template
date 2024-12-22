#ifndef DataFormats_SiPixelClusterSoA_test_alpaka_Hits_test_h
#define DataFormats_SiPixelClusterSoA_test_alpaka_Hits_test_h

#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsSoA.h"
#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisSoA.h"
#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersSoA.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::Splitting {

  template <typename TrackerTraits>
  void runKernels(TrackingRecHitSoAView<TrackerTraits>& hits,
                  SiPixelDigisSoAView& digis,
                  SiPixelClustersSoAView& clusters,
                  Queue& queue);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::Splitting


#endif  // DataFormats_SiPixelClusterSoA_test_alpaka_Hits_test_h
