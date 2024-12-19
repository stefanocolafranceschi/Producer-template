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

  namespace testSiPixelClustersSoA {

    struct ShowKernel {
        template <typename TAcc, typename = std::enable_if_t<isAccelerator<TAcc>>>
        ALPAKA_FN_ACC void operator()(TAcc const& acc, SiPixelClustersSoAConstView soa) const {
            if (cms::alpakatools::once_per_grid(acc)) {
                printf("well this works!\n");
            }
/*
            printf("Inspecting available members:\n");
            printf("moduleStart[0]: %u\n", soa.moduleStart(0));
            printf("clusInModule[0]: %u\n", soa.clusInModule(0));
            printf("moduleId[0]: %u\n", soa.moduleId(0));
            printf("clusModuleStart[0]: %u\n", soa.clusModuleStart(0));
*/
for (size_t i = 0; i < 5; ++i) {

            printf("moduleStart[i]: %u\n", soa.moduleStart(i));
            printf("clusInModule[i]: %u\n", soa.clusInModule(i));
            printf("moduleId[i]: %u\n", soa.moduleId(i));
            printf("clusModuleStart[i]: %u\n", soa.clusModuleStart(i));

}

        }
    };

    // Reintroduce template specialization
    void runKernels(SiPixelClustersSoAView& clusters, Queue& queue) {
        uint32_t items = 64;
        uint32_t groups = divide_up_by(2, items);  // Use actual cluster count
        auto workDiv = make_workdiv<Acc1D>(groups, items);

        // Launch kernels
        alpaka::exec<Acc1D>(queue, workDiv, ShowKernel{}, clusters);
    }

  }  // namespace testSiPixelClustersSoA

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
