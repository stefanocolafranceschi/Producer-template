#include <type_traits>
#include <alpaka/alpaka.hpp>

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
                printf("Inspecting available members:\n");
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
