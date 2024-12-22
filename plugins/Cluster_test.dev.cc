#include <type_traits>

#include <alpaka/alpaka.hpp>

#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsDevice.h"
#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsSoA.h"
#include "DataFormats/TrackingRecHitSoA/interface/alpaka/TrackingRecHitsSoACollection.h"

#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisDevice.h"
#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisSoA.h"
#include "DataFormats/SiPixelDigiSoA/interface/alpaka/SiPixelDigisSoACollection.h"

#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersDevice.h"
#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersSoA.h"
#include "DataFormats/SiPixelClusterSoA/interface/alpaka/SiPixelClustersSoACollection.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

using namespace alpaka;

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;
  namespace Splitting {

    template <typename TrackerTraits>
    struct Printout {
      template <typename TAcc, typename = std::enable_if_t<isAccelerator<TAcc>>>
      ALPAKA_FN_ACC void operator()(TAcc const& acc, 
                                    TrackingRecHitSoAConstView<TrackerTraits> hitView, 
                                    SiPixelDigisSoAConstView digiView,
                                    SiPixelClustersSoAConstView clustersView) const {
        
        // Output information from TrackingRecHitSoAConstView
        if (cms::alpakatools::once_per_grid(acc)) {
          printf("TrackingRecHits Info:\n");
          printf("nbins = %d\n", hitView.phiBinner().nbins());
          printf("offsetBPIX = %d\n", hitView.offsetBPIX2());
          printf("nHits = %d\n", hitView.metadata().size());
        }

        // Print debug info for hits
        for (uint32_t i : cms::alpakatools::uniform_elements(acc, 10)) {
          printf("Hit iPhi %d -> %d\n", i, hitView[i].iphi());
        }

        // Output information from SiPixelDigisSoAConstView
        if (cms::alpakatools::once_per_grid(acc)) {
          printf("SiPixelDigis Info:\n");
          printf("nDigis = %d\n", digiView.metadata().size());
        }

        // Print debug info for digis
        for (uint32_t j : cms::alpakatools::uniform_elements(acc, 10)) {
          printf("Digi Module %d -> ADC %d\n", j, digiView[j].adc());
        }

        // Output information from SiPixelClustersSoAConstView
        if (cms::alpakatools::once_per_grid(acc)) {
          printf("SiPixelClusters Info:\n");
          printf("nClusters = %d\n", clustersView.metadata().size());
        }

        // Print debug info for clusters
        for (uint32_t k : cms::alpakatools::uniform_elements(acc, 10)) {
          //printf("Cluster %d -> Size %d\n", k, clustersView[k].size());
        }
      }
    };




    template <typename TrackerTraits>
    void runKernels(TrackingRecHitSoAView<TrackerTraits>& hitView,
                    SiPixelDigisSoAView& digiView,
                    SiPixelClustersSoAView& clustersView,
                    Queue& queue) {
        uint32_t items = 64;
        uint32_t groupsHits = divide_up_by(hitView.metadata().size(), items);
        uint32_t groupsDigis = divide_up_by(digiView.metadata().size(), items);
        uint32_t groupsClusters = divide_up_by(clustersView.metadata().size(), items);
        
        uint32_t groups = std::max({groupsHits, groupsDigis, groupsClusters});  // Ensure work division covers all three views

        auto workDiv = make_workdiv<Acc1D>(groups, items);

        // Kernel execution for hits, digis, and clusters
        alpaka::exec<Acc1D>(queue, workDiv, Printout<TrackerTraits>{}, hitView, digiView, clustersView);
    }

    template void runKernels<pixelTopology::Phase1>(TrackingRecHitSoAView<pixelTopology::Phase1>& hitView,
                                                    SiPixelDigisSoAView& digiView,
                                                    SiPixelClustersSoAView& clustersView,
                                                    Queue& queue);

    template void runKernels<pixelTopology::Phase2>(TrackingRecHitSoAView<pixelTopology::Phase2>& hitView,
                                                    SiPixelDigisSoAView& digiView,
                                                    SiPixelClustersSoAView& clustersView,
                                                    Queue& queue);
  }  // namespace Splitting
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE