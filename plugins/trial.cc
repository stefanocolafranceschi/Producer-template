#ifdef ALPAKA_ACC_GPU_CUDA_ENABLED  // Ensure this file is only used for CUDA

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

// system include files
#include "TFile.h"
#include "TString.h"


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

class trial : public edm::stream::EDProducer<> {
public:
  explicit trial(const edm::ParameterSet&);
  ~trial() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  std::string configString_;
  uint32_t nHits_;
  int32_t offset_;
  TFile* rootFile_;
  std::vector<Device> devices_;
};

trial::trial(const edm::ParameterSet& iConfig)
    : configString_(iConfig.getParameter<std::string>("configString")),
      nHits_(iConfig.getParameter<uint32_t>("nHits")),
      offset_(iConfig.getParameter<int32_t>("offset")),
      rootFile_(nullptr) {
  // Initialize ROOT file
  rootFile_ = new TFile("config_output.root", "RECREATE");

  // Register the product that this module will produce
  produces<std::vector<int>>("outputHits");
}

trial::~trial() {
  if (rootFile_) {
    rootFile_->Close();
    delete rootFile_;
  }
}

void trial::beginStream(edm::StreamID) {
  // Initialize devices
  devices_ = cms::alpakatools::devices<Platform>();
  if (devices_.empty()) {
    edm::LogError("trial") << "No devices available for the backend. The test will be skipped.";
    return;
  }
  edm::LogInfo("trial") << "Found " << devices_.size() << " device(s).";
}

void trial::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if (devices_.empty()) {
    edm::LogWarning("trial") << "Skipping event because no devices are available.";
    return;
  }

  // Extract data from the event
  int32_t eventOffset = iEvent.id().event();

  for (const auto& device : devices_) {
    Queue queue(device);

    // Define moduleStart data
    auto moduleStartH =
        cms::alpakatools::make_host_buffer<uint32_t[]>(queue, pixelTopology::Phase1::numberOfModules + 1);
    for (size_t i = 0; i < pixelTopology::Phase1::numberOfModules + 1; ++i) {
      moduleStartH[i] = i * 2;
    }
    auto moduleStartD =
        cms::alpakatools::make_device_buffer<uint32_t[]>(queue, pixelTopology::Phase1::numberOfModules + 1);
    alpaka::memcpy(queue, moduleStartD, moduleStartH);

    // Pass event-based offset to kernel
    TrackingRecHitsSoACollection<pixelTopology::Phase1> tkhit(queue, nHits_, eventOffset, moduleStartD.data());

    // Execute kernels
    testTrackingRecHitSoA::runKernels<pixelTopology::Phase1>(tkhit.view(), queue);

    // Update data from device back to host
    tkhit.updateFromDevice(queue);

    TrackingRecHitHost<pixelTopology::Phase1> host_collection =
        cms::alpakatools::CopyToHost<TrackingRecHitDevice<pixelTopology::Phase1, Device>>::copyAsync(queue, tkhit);

    // Wait for kernel completion
    alpaka::wait(queue);

    // Perform assertions
    assert(tkhit.nHits() == nHits_);
    assert(tkhit.offsetBPIX2() == eventOffset);
    assert(tkhit.nHits() == host_collection.nHits());
    assert(tkhit.offsetBPIX2() == host_collection.offsetBPIX2());
  }
}

void trial::endStream() {
  edm::LogInfo("trial") << "Processing completed.";
}

void trial::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("configString", "default")->setComment("Configuration string to store in ROOT file");
  desc.add<uint32_t>("nHits", 100)->setComment("Number of hits for the test");
  desc.add<int32_t>("offset", 0)->setComment("Offset for hits");
  descriptions.add("trial", desc);
}

DEFINE_FWK_MODULE(trial);

#endif  // ALPAKA_ACC_GPU_CUDA_ENABLED
