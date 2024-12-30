#ifdef ALPAKA_ACC_GPU_CUDA_ENABLED  // Ensure this file is only used for CUDA

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/stringize.h"

// system include files
#include "TFile.h"
#include "TString.h"

#include <cstdlib>
#include <unistd.h>

#include <alpaka/alpaka.hpp>

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsDevice.h"
#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsHost.h"
#include "DataFormats/TrackingRecHitSoA/interface/alpaka/TrackingRecHitsSoACollection.h"

#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisDevice.h"
#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisHost.h"
#include "DataFormats/SiPixelDigiSoA/interface/alpaka/SiPixelDigisSoACollection.h"

#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersDevice.h"
#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersHost.h"
#include "DataFormats/SiPixelClusterSoA/interface/alpaka/SiPixelClustersSoACollection.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/VertexSoA/interface/ZVertexSoA.h"
#include "DataFormats/VertexSoA/interface/ZVertexHost.h"
#include "DataFormats/VertexSoA/interface/ZVertexDevice.h"
#include "DataFormats/VertexSoA/interface/alpaka/ZVertexSoACollection.h"

#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/SimplePixelTopology.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

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

  edm::EDGetTokenT<ALPAKA_ACCELERATOR_NAMESPACE::SiPixelClustersSoACollection> clusterToken_;
  edm::EDGetTokenT<ALPAKA_ACCELERATOR_NAMESPACE::SiPixelDigisSoACollection> digisToken_;
  edm::EDGetTokenT<ALPAKA_ACCELERATOR_NAMESPACE::TrackingRecHitsSoACollection<pixelTopology::Phase1>> recHitsToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> candidateToken_;
  edm::EDGetTokenT<ALPAKA_ACCELERATOR_NAMESPACE::ZVertexSoACollection> zVertexToken_;
  edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> const tTrackingGeom_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> const tTrackerTopo_;

  const double ptMin_;
};


trial::trial(const edm::ParameterSet& iConfig)
    : configString_(iConfig.getParameter<std::string>("configString")),
      nHits_(iConfig.getParameter<uint32_t>("nHits")),
      offset_(iConfig.getParameter<int32_t>("offset")),
      rootFile_(nullptr),
      clusterToken_(consumes<ALPAKA_ACCELERATOR_NAMESPACE::SiPixelClustersSoACollection>(edm::InputTag("siPixelClusters"))),
      digisToken_(consumes<ALPAKA_ACCELERATOR_NAMESPACE::SiPixelDigisSoACollection>(edm::InputTag("SiPixelDigisSoA"))),
      recHitsToken_(consumes<ALPAKA_ACCELERATOR_NAMESPACE::TrackingRecHitsSoACollection<pixelTopology::Phase1>>(edm::InputTag("TrackingRecHitsSoA"))),
      candidateToken_(consumes<edm::View<reco::Candidate>>(edm::InputTag("candidateInput"))),
      zVertexToken_(consumes<ALPAKA_ACCELERATOR_NAMESPACE::ZVertexSoACollection>(edm::InputTag("ZVertex"))),
      tTrackingGeom_(esConsumes<GlobalTrackingGeometry, GlobalTrackingGeometryRecord>()),
      tTrackerTopo_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
      ptMin_(iConfig.getParameter<double>("ptMin"))
        {
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

    edm::Handle<ALPAKA_ACCELERATOR_NAMESPACE::SiPixelClustersSoACollection> clustersHandle;
    iEvent.getByToken(clusterToken_, clustersHandle);
    if (!clustersHandle.isValid()) {
        edm::LogError("trial") << "Could not retrieve SiPixelClusters.";
        return;
    }

    edm::Handle<ALPAKA_ACCELERATOR_NAMESPACE::SiPixelDigisSoACollection> digisHandle;
    iEvent.getByToken(digisToken_, digisHandle);
    if (!digisHandle.isValid()) {
        edm::LogError("trial") << "Could not retrieve SiPixelDigisSoA.";
        return;
    }

    edm::Handle<ALPAKA_ACCELERATOR_NAMESPACE::TrackingRecHitsSoACollection<pixelTopology::Phase1>> recHitsHandle;
    iEvent.getByToken(recHitsToken_, recHitsHandle);
    if (!recHitsHandle.isValid()) {
        edm::LogError("trial") << "Could not retrieve TrackingRecHitsSoA.";
        return;
    }        

    edm::Handle<edm::View<reco::Candidate>> candidatesHandle;
    iEvent.getByToken(candidateToken_, candidatesHandle);
    if (!candidatesHandle.isValid()) {
        edm::LogError("trial") << "Could not retrieve Candidate collection.";
        return;
    }

    std::vector<CandidateGPUData> gpuCandidates;
    gpuCandidates.reserve(candidatesHandle->size());


    // Iterate over the candidates in the collection, filter them based on a minimum transverse momentum (ptMin_),
    // and prepare their data in a GPU-compatible format (CandidateGPUData) for further processing on the device.    
    for (const auto& candidate : *candidatesHandle) {
        if (candidate.pt() > ptMin_) {
            CandidateGPUData data = {
                static_cast<float>(candidate.px()),  // Cast to float
                static_cast<float>(candidate.py()),  // Cast to float
                static_cast<float>(candidate.pz()),  // Cast to float
                static_cast<float>(candidate.pt()),  // Cast to float
                static_cast<float>(candidate.eta()), // Cast to float
                static_cast<float>(candidate.phi())  // Cast to float
            };
            gpuCandidates.push_back(data);
        }
    }
    auto const& candidates = *candidatesHandle;

    edm::Handle<ALPAKA_ACCELERATOR_NAMESPACE::ZVertexSoACollection> zVertexHandle;
    iEvent.getByToken(zVertexToken_, zVertexHandle);
    if (!zVertexHandle.isValid()) {
        edm::LogError("trial") << "Could not retrieve zVertexHandle";
        return;
    }
    
    // Process TrackingRecHitsSoACollection
    const auto& recHits = *recHitsHandle;
    size_t nHits = recHits.nHits();
    std::cout << "Number of hits: " << recHits.nHits() << std::endl;

    // Process SiPixelDigisSoACollection
    const auto& digis = *digisHandle;
    size_t nDigis = digis.nDigis();
    std::cout << "Number of digis: " << nDigis << std::endl;

    // Process SiPixelClustersSoACollection
    const auto& clusters = *clustersHandle;
    uint32_t nClusters = clusters.nClusters();
    std::cout << "Total clusters in this event: " << nClusters << std::endl;

    // Process ZVertexSoACollection
    const auto& zVertices = *zVertexHandle;
    uint32_t nVertices = zVertices.view().nvFinal();
    std::cout << "Number of Vertices: " << nVertices << std::endl;


    // Use event ID as the offset
    int32_t eventOffset = iEvent.id().event();
    std::cout << "Event offset: " << eventOffset << std::endl;

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
        alpaka::wait(queue);            // Ensure the data copy is complete


        // ------------- CREATE DEVICE BUFFERS -------------------------------
        /* RecHits
           the TrackingRecHitsSoACollection is an alias for: TrackingRecHitDevice (gpu) 
                                                            TrackingRecHitHost (cpu)  */
        TrackingRecHitsSoACollection<pixelTopology::Phase1> tkhit(queue, nHits, eventOffset, moduleStartD.data());
        //- - - - - - - - - - - - - - - - - - -


        /* Digis 
        the SiPixelDigisSoACollection is an alias for: SiPixelDigisDevice (gpu) or 
                                                          SiPixelDigisHost (cpu)
        but it's not templated so <pixelTopology> won't work
        I could also: SiPixelDigisDevice<Device> digisDevice(nDigis, queue); */
        SiPixelDigisSoACollection tkdigi(nDigis, queue);
        tkdigi.setNModules(pixelTopology::Phase1::numberOfModules);         // Set additional metadata
        //- - - - - - - - - - - - - - - - - - -


        /* Clusters
           the SiPixelClustersSoACollection is an alias for: SiPixelClustersDevice (gpu) 
                                                            SiPixelClustersHost (cpu)  */
        SiPixelClustersSoACollection tkclusters(nClusters, queue);
        // It seems the above class has no topology and no Modules.. not sure why


        /* Candidates
           I create Alpaka host buffer and device buffer (for GPU)*/
        auto gpuCandidatesHost = cms::alpakatools::make_host_buffer<CandidateGPUData[]>(queue, gpuCandidates.size());
        auto gpuCandidatesDevice = cms::alpakatools::make_device_buffer<CandidateGPUData[]>(queue, gpuCandidates.size());


        ZVertexSoACollection tkvertices(queue);

        //- - - - - - - - - - - - - - - - - - -



        // ------------- COPY FROM HOST TO DEVICE BUFFERS -------------------------------
        alpaka::memcpy(queue, tkhit.buffer(), recHitsHandle->buffer());
        alpaka::memcpy(queue, tkdigi.buffer(), digisHandle->buffer());
        alpaka::memcpy(queue, tkclusters.buffer(), clustersHandle->buffer());
        alpaka::memcpy(queue, tkvertices.buffer(), zVertexHandle->buffer());  // Transfer ZVertex data from host to device
        std::copy(gpuCandidates.begin(), gpuCandidates.end(), gpuCandidatesHost.data());  // Copy data from std::vector to Alpaka host buffer
        alpaka::memcpy(queue, gpuCandidatesDevice, gpuCandidatesHost);          // Transfer data from host to device
        alpaka::wait(queue);  // Ensure the transfer is complete


        // Run the kernel with GPU candidates
        Splitting::runKernels<pixelTopology::Phase1>(
            tkhit.view(), tkdigi.view(), tkclusters.view(), tkvertices.view(), gpuCandidatesDevice.data(), gpuCandidates.size(), queue
        );


        // Update from device to host (RecHits and Digis)
        tkhit.updateFromDevice(queue);
        TrackingRecHitHost<pixelTopology::Phase1> hostRecHits = cms::alpakatools::CopyToHost<TrackingRecHitDevice<pixelTopology::Phase1, Device>>::copyAsync(queue, tkhit);
        SiPixelDigisHost digisHost = cms::alpakatools::CopyToHost<SiPixelDigisDevice<Device>>::copyAsync(queue, tkdigi);
        SiPixelClustersHost clustersHost = cms::alpakatools::CopyToHost<SiPixelClustersDevice<Device>>::copyAsync(queue, tkclusters);

        alpaka::wait(queue);
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
