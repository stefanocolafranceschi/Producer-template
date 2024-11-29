#include <type_traits>

#include <alpaka/alpaka.hpp>

#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsDevice.h"
#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsSoA.h"
#include "DataFormats/TrackingRecHitSoA/interface/alpaka/TrackingRecHitsSoACollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "Hits_test.h"

using namespace alpaka;

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;
  namespace testTrackingRecHitSoA {

    template <typename TrackerTraits>    
    struct TestFillKernel {
      template <typename TAcc, typename = std::enable_if_t<isAccelerator<TAcc>>>
      ALPAKA_FN_ACC void operator()(TAcc const& acc, TrackingRecHitSoAView<TrackerTraits> soa) const {

        const uint32_t i(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        const uint32_t j(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);

        //printf("Thread Index [i] in the Grid: %u\n", i);
        //printf("Thread Index [j] in the Block: %u\n", j);

        // Initialization of global/shared data
        // once_per_grid does what it says
        // (put a random number in here and all threads will have the same random number :-()
        if (cms::alpakatools::once_per_grid(acc)) {
          soa.offsetBPIX2() = 22;
        }

        // Using Linear Congruential Generator
        // i and j of each thread used as seed 
        uint32_t seed = i + j;
        float randomValue = (seed * 1664525u + 1013904223u) % 10000 / 1000.0f;  // Linear congruential generator
        soa[i].xLocal() = randomValue;

        seed = i + j + 1000 ;
        randomValue = (seed * 1664525u + 1013904223u) % 10000 / 1000.0f;  // Linear congruential generator
        soa[i].yLocal() = randomValue;

        soa[i].iphi() = i % 10;
        soa.hitsLayerStart()[j] = j;

        // Print out the initialized values for debugging from each thread!
        //printf("Thread [i = %u, j = %u] Initialized:\n", i, j);
        //printf("  xLocal = %f\n", soa[i].xLocal());
        //printf("  yLocal = %f\n", soa[i].yLocal());
        //printf("  iphi = %d\n", soa[i].iphi());
        //printf("  hitsLayerStart[%u] = %d\n", j, soa.hitsLayerStart()[j]);

      }
    };

    template <typename TrackerTraits>
    struct ShowKernel {
      template <typename TAcc, typename = std::enable_if_t<isAccelerator<TAcc>>>
      ALPAKA_FN_ACC void operator()(TAcc const& acc, TrackingRecHitSoAConstView<TrackerTraits> soa) const {

        if (cms::alpakatools::once_per_grid(acc)) {
          printf("nbins = %d\n", soa.phiBinner().nbins());
          printf("offsetBPIX = %d\n", soa.offsetBPIX2());
          printf("nHits = %d\n", soa.metadata().size());
          //printf("hitsModuleStart[28] = %d\n", soa[28].hitsModuleStart());
        }

        // can be increased to soa.nHits() for debugging
        for (uint32_t i : cms::alpakatools::uniform_elements(acc, 10)) {
          printf("iPhi %d -> %d\n", i, soa[i].iphi());
        }
      }
    };


    template <typename TrackerTraits>
    struct CheckXLocalKernel {
      template <typename TAcc, typename = std::enable_if_t<isAccelerator<TAcc>>>
      ALPAKA_FN_ACC void operator()(TAcc const& acc, TrackingRecHitSoAConstView<TrackerTraits> soa) const {

        // Get the thread index
        const uint32_t i(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        const uint32_t j(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);

        // Define the range for xLocal (this is just an example range)
        const float minValidXLocal = 0.0f;
        const float maxValidXLocal = 1.0f;

        // Check if the xLocal value is within the specified range
        float xLocalValue = soa[i].xLocal();
        if (xLocalValue < minValidXLocal || xLocalValue > maxValidXLocal) {
          // Print a warning if xLocal is outside the valid range
          printf("Warning: Thread [i = %u, j = %u] has invalid xLocal value: %f\n", i, j, xLocalValue);
        }
      }
    };


    template <typename TrackerTraits>
    void runKernels(TrackingRecHitSoAView<TrackerTraits>& view, Queue& queue) {
      uint32_t items = 64;
      uint32_t groups = divide_up_by(view.metadata().size(), items);
      auto workDiv = make_workdiv<Acc1D>(groups, items);
      alpaka::exec<Acc1D>(queue, workDiv, TestFillKernel<TrackerTraits>{}, view);
      alpaka::exec<Acc1D>(queue, workDiv, ShowKernel<TrackerTraits>{}, view);
      alpaka::exec<Acc1D>(queue, workDiv, CheckXLocalKernel<TrackerTraits>{}, view);
    }

    template void runKernels<pixelTopology::Phase1>(TrackingRecHitSoAView<pixelTopology::Phase1>& view, Queue& queue);
    template void runKernels<pixelTopology::Phase2>(TrackingRecHitSoAView<pixelTopology::Phase2>& view, Queue& queue);

  }  // namespace testTrackingRecHitSoA
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
