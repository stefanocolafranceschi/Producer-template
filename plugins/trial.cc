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

#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersDevice.h"
#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersHost.h"
#include "DataFormats/SiPixelClusterSoA/interface/alpaka/SiPixelClustersSoACollection.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "FWCore/Utilities/interface/stringize.h"
#include "Geometry/CommonTopologies/interface/SimplePixelTopology.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/devices.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "Cluster_test.h"

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

  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> clusterToken_;

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
  clusterToken_ = consumes<edmNew::DetSetVector<SiPixelCluster>>(edm::InputTag("siPixelClusters"));

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

  // Handle to access SiPixelCluster
  edm::Handle<edmNew::DetSetVector<SiPixelCluster>> clustersHandle;

  // Retrieve clusters from the event
  iEvent.getByToken(clusterToken_, clustersHandle);

  // Check if the handle is valid
  if (!clustersHandle.isValid()) {
    edm::LogError("trial") << "Could not retrieve siPixelClusters.";
    return;
  }

  uint32_t nClusters = 0;
  for (const auto& detSet : *clustersHandle) {
    nClusters += detSet.size();  // Correctly count clusters for edmNew::DetSetVector
  }

  // Print total clusters
  std::cout << "Total clusters in this event: " << nClusters << std::endl;

  // Use event ID as the offset
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

    // Create device buffers for SiPixelClustersSoA
    SiPixelClustersDevice<Device> clustersSoA(nClusters, queue);

    // Debugging output
    std::cout << "Type of clustersSoA: " << typeid(clustersSoA).name() << std::endl;
    std::cout << "Type of clustersSoA.view(): " << typeid(clustersSoA.view()).name() << std::endl;

    // Run the kernel on the clusters
    testSiPixelClustersSoA::runKernels(clustersSoA.view(), queue);

    // Update data back from device to host (if needed)
    //clustersSoA.updateFromDevice(queue);

    // Wait for kernel completion
    alpaka::wait(queue);

    // Verify results (optional)
    // assert(clustersSoA.nClusters() == nClusters);
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
