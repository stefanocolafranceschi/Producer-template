#include <cstdlib>
#include <unistd.h>

#include <alpaka/alpaka.hpp>

#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsDevice.h"
#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsHost.h"
#include "DataFormats/TrackingRecHitSoA/interface/alpaka/TrackingRecHitsSoACollection.h"
#include "FWCore/Utilities/interface/stringize.h"
#include "Geometry/CommonTopologies/interface/SimplePixelTopology.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/devices.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "Hits_test.h"

using namespace ALPAKA_ACCELERATOR_NAMESPACE;

class HitsTestRunner {
public:
    HitsTestRunner(const std::vector<Device>& devices, uint32_t nHits, int32_t offset)
        : nHits_(nHits), offset_(offset) {
        devices_ = cms::alpakatools::devices<Platform>();
        if (devices_.empty()) {
            std::cerr << "No devices available for the " EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE) " backend, "
                      << "the test will be skipped.\n";
            exit(EXIT_FAILURE);
        }

        printDeviceInfo();
    }

    void run() {
        for (const auto& device : devices_) {
            Queue queue(device);
            runOnDevice(queue);
        }
    }

private:
    void printDeviceInfo() const {
        std::cout << "Found " << devices_.size() << " device(s)." << std::endl;
        for (size_t i = 0; i < devices_.size(); ++i) {
            const auto& device = devices_[i];
            std::cout << "Device " << i << " type: " << typeid(device).name() << std::endl;
        }
    }

    void runOnDevice(Queue& queue) {
        // Allocate and initialize memory
        auto moduleStartH =
            cms::alpakatools::make_host_buffer<uint32_t[]>(queue, pixelTopology::Phase1::numberOfModules + 1);
        for (size_t i = 0; i < pixelTopology::Phase1::numberOfModules + 1; ++i) {
            moduleStartH[i] = i * 2;
        }
        auto moduleStartD =
            cms::alpakatools::make_device_buffer<uint32_t[]>(queue, pixelTopology::Phase1::numberOfModules + 1);
        alpaka::memcpy(queue, moduleStartD, moduleStartH);

        TrackingRecHitsSoACollection<pixelTopology::Phase1> tkhit(queue, nHits_, offset_, moduleStartD.data());

        // Execute kernels
        testTrackingRecHitSoA::runKernels<pixelTopology::Phase1>(tkhit.view(), queue);

        // Synchronize and check results
        tkhit.updateFromDevice(queue);
#if defined ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED or defined ALPAKA_ACC_CPU_B_TBB_T_SEQ_ENABLED
        TrackingRecHitHost<pixelTopology::Phase1> const& host_collection = tkhit;
#else
        TrackingRecHitHost<pixelTopology::Phase1> host_collection =
            cms::alpakatools::CopyToHost<TrackingRecHitDevice<pixelTopology::Phase1, Device> >::copyAsync(queue, tkhit);
#endif
        alpaka::wait(queue);

        // Verify results
        assert(tkhit.nHits() == nHits_);
        assert(tkhit.offsetBPIX2() == 22);  // Set in the kernel
        assert(tkhit.nHits() == host_collection.nHits());
        assert(tkhit.offsetBPIX2() == host_collection.offsetBPIX2());
    }

    std::vector<Device> devices_;
    uint32_t nHits_;
    int32_t offset_;
};