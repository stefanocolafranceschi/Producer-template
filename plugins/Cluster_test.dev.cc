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

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/JetReco/interface/Jet.h"

#include "Cluster_test.h"

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
                                SiPixelClustersSoAConstView clustersView,
                                const CandidateGPUData* candidates,  // Updated to CandidateGPUData
                                size_t nCandidates) const {         // Added number of candidates
 

        // Print debug info for RecHits -----------------------------------
        if (cms::alpakatools::once_per_grid(acc)) {
          printf("TrackingRecHits Info:\n");
          printf("nbins = %d\n", hitView.phiBinner().nbins());
          printf("offsetBPIX = %d\n", hitView.offsetBPIX2());
          printf("nHits = %d\n", hitView.metadata().size());
        }
        for (uint32_t i : cms::alpakatools::uniform_elements(acc, 10)) {
          printf("Hit %d -> xLocal: %.2f, yLocal: %.2f, xerrLocal: %.2f, yerrLocal: %.2f, "
                 "xGlobal: %.2f, yGlobal: %.2f, zGlobal: %.2f, rGlobal: %.2f, iPhi: %d, "
                 "charge: %d, isBigX: %d, isOneX: %d, isBigY: %d, isOneY: %d, qBin: %d, "
                 "clusterSizeX: %d, clusterSizeY: %d, detectorIndex: %d\n",
                 i,
                 hitView[i].xLocal(),
                 hitView[i].yLocal(),
                 hitView[i].xerrLocal(),
                 hitView[i].yerrLocal(),
                 hitView[i].xGlobal(),
                 hitView[i].yGlobal(),
                 hitView[i].zGlobal(),
                 hitView[i].rGlobal(),
                 hitView[i].iphi(),
                 hitView[i].chargeAndStatus().charge,
                 hitView[i].chargeAndStatus().status.isBigX,
                 hitView[i].chargeAndStatus().status.isOneX,
                 hitView[i].chargeAndStatus().status.isBigY,
                 hitView[i].chargeAndStatus().status.isOneY,
                 hitView[i].chargeAndStatus().status.qBin,
                 hitView[i].clusterSizeX(),
                 hitView[i].clusterSizeY(),
                 hitView[i].detectorIndex());
        }


        // Print debug info for digis -----------------------------------
        if (cms::alpakatools::once_per_grid(acc)) {
          printf("SiPixelDigis Info:\n");
          printf("nDigis = %d\n", digiView.metadata().size());
        }
        for (uint32_t j : cms::alpakatools::uniform_elements(acc, 10)) {
          uint16_t x = digiView[j].xx();
          uint16_t y = digiView[j].yy();
          uint16_t adc = digiView[j].adc();
          printf("Digi %d -> x: %d, y: %d, ADC: %d\n", j, x, y, adc);
        }


        // Print debug info for Clusters -----------------------------------
        if (cms::alpakatools::once_per_grid(acc)) {
          printf("SiPixelClusters Info:\n");
          printf("nClusters = %d\n", clustersView.metadata().size());
        }
        for (uint32_t k : cms::alpakatools::uniform_elements(acc, 10)) {
            printf("Cluster %d -> moduleStart: %d, clusInModule: %d, moduleId: %d, clusModuleStart: %d\n",
                   k,
                   clustersView[k].moduleStart(),
                   clustersView[k].clusInModule(),
                   clustersView[k].moduleId(),
                   clustersView[k].clusModuleStart());
        }

        // Iterate over all clusters (assuming clusters are indexed from 0 to nClusters-1)
        for (uint32_t clusterIdx : cms::alpakatools::uniform_elements(acc, clustersView.metadata().size())) {
            // Temporary storage for cluster properties
            uint32_t minX = std::numeric_limits<uint32_t>::max();
            uint32_t maxX = 0;
            uint32_t minY = std::numeric_limits<uint32_t>::max();
            uint32_t maxY = 0;
            uint32_t totalADC = 0;
            uint32_t numPixels = 0;

            // Iterate over all digis to find those belonging to the current cluster
            for (uint32_t j : cms::alpakatools::uniform_elements(acc, digiView.metadata().size())) {
                if (static_cast<uint32_t>(digiView[j].clus()) == clusterIdx) { // Fixed comparison
                    uint16_t x = digiView[j].xx();
                    uint16_t y = digiView[j].yy();
                    uint16_t adc = digiView[j].adc();

                    // Update cluster properties
                    minX = std::min(minX, (uint32_t)x);
                    maxX = std::max(maxX, (uint32_t)x);
                    minY = std::min(minY, (uint32_t)y);
                    maxY = std::max(maxY, (uint32_t)y);
                    totalADC += adc;
                    numPixels++;
                }
            }

            // Print cluster properties
            if (numPixels > 0) { // Only print clusters that contain pixels
                printf("Cluster %d -> Pixels: %d, Total ADC: %d, Bounds: x[%d-%d], y[%d-%d]\n",
                       clusterIdx, numPixels, totalADC, minX, maxX, minY, maxY);
            }
        }
      }
    };


/*
template <typename TrackerTraits>
struct JetSplit {
  // Define ClusterProperties to store cluster-related information
  struct ClusterProperties {
    uint32_t minX = std::numeric_limits<uint32_t>::max();
    uint32_t minY = std::numeric_limits<uint32_t>::max();
    uint32_t maxX = 0;
    uint32_t maxY = 0;
    uint32_t charge = 0;

    uint32_t sizeX() const { return maxX - minX + 1; }
    uint32_t sizeY() const { return maxY - minY + 1; }
  };

  // Main operator function
  template <typename TAcc, typename = std::enable_if_t<isAccelerator<TAcc>>>
  ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                TrackingRecHitSoAConstView<TrackerTraits> hitView,
                                SiPixelDigisSoAConstView digiView,
                                SiPixelClustersSoAConstView clustersView,
                                const std::vector<reco::Candidate>& cores, // Jet collection
                                const GlobalPoint& ppv, // Primary vertex position
                                const edm::ESHandle<GeometricSearchTracker>& det,
                                float ptMin_, float deltaR_) const {

      // Iterate over clusters
      for (uint32_t clusterIdx : cms::alpakatools::uniform_elements(acc, clustersView.metadata().size())) {
        // Compute properties of the cluster from digis
        auto clusterProps = computeClusterProperties(acc, digiView, clusterIdx);

        // Placeholder: Get cluster position (this might come from the cluster data)
        GlobalPoint cPos = getClusterPosition(clusterIdx, digiView);  // Assume this is defined

        // Compute the cluster direction (vector from cluster to ppv)
        GlobalVector clusterDir = cPos - ppv;

        // Iterate over jets (cores) to check deltaR with the cluster


        for (unsigned int ji = 0; ji < cores.size(); ji++) {
          if (ji < cores.size()) {
            const reco::Candidate& jet = cores.at(ji);
            GlobalVector jetDir(static_cast<double>(jet.px()), 
                                static_cast<double>(jet.py()), 
                                static_cast<double>(jet.pz()));

            // Calculate deltaR between jet and cluster
            if (Geom::deltaR(jetDir, clusterDir) < deltaR_) {
              // If deltaR condition met, check if the cluster should be split
              if (shouldSplit(clusterProps, jetDir)) {
                // Split the cluster into sub-clusters
                auto subClusters = splitCluster(acc, clusterProps, jetDir, digiView, clusterIdx, det);

                // Debug output for sub-clusters
                for (const auto& subCluster : subClusters) {
                  printf("Sub-cluster -> minX: %d, minY: %d, maxX: %d, maxY: %d, charge: %d\n",
                         subCluster.minX, subCluster.minY, subCluster.maxX, subCluster.maxY, subCluster.charge);
                }
              }
            }
          }
        }
      }

  // Function to compute cluster properties
  template <typename TAcc>
  ALPAKA_FN_ACC ClusterProperties computeClusterProperties(
      TAcc const& acc,
      const SiPixelDigisSoAConstView& digiView,
      uint32_t clusterIdx) const {
    ClusterProperties props;

    for (uint32_t j : cms::alpakatools::uniform_elements(acc, digiView.metadata().size())) {
      if (static_cast<uint32_t>(digiView[j].clus()) == clusterIdx) {
        uint16_t x = digiView[j].xx();
        uint16_t y = digiView[j].yy();
        uint16_t adc = digiView[j].adc();

        props.minX = std::min(props.minX, static_cast<uint32_t>(x));
        props.minY = std::min(props.minY, static_cast<uint32_t>(y));
        props.maxX = std::max(props.maxX, static_cast<uint32_t>(x));
        props.maxY = std::max(props.maxY, static_cast<uint32_t>(y));
        props.charge += adc;
      }
    }
    return props;
  }

  // Function to decide whether a cluster should be split
  ALPAKA_FN_ACC bool shouldSplit(const ClusterProperties& cluster, const GlobalVector& jetDir) const {
    // Placeholder logic: split if charge or size exceeds threshold
    return cluster.charge > 500 && (cluster.sizeX() > 2 || cluster.sizeY() > 2);
  }

  // Function to split a cluster into sub-clusters
  template <typename TAcc>
  ALPAKA_FN_ACC std::vector<ClusterProperties> splitCluster(
      TAcc const& acc,
      const ClusterProperties& clusterProps,
      const GlobalVector& jetDir,
      const SiPixelDigisSoAConstView& digiView,
      uint32_t clusterIdx,
      const edm::ESHandle<GeometricSearchTracker>& det) const {
        std::vector<ClusterProperties> subClusters;

        for (uint32_t j : cms::alpakatools::uniform_elements(acc, digiView.metadata().size())) {
          if (static_cast<uint32_t>(digiView[j].clus()) == clusterIdx) {
            uint16_t x = digiView[j].xx();
            uint16_t y = digiView[j].yy();
            uint16_t adc = digiView[j].adc();

            // Calculate distance from jet direction
            float dx = x - jetDir.x();
            float dy = y - jetDir.y();
            float distance = std::sqrt(dx * dx + dy * dy);

            if (distance < 2.0f) { // Example threshold
              // Check if this pixel belongs to an existing sub-cluster
              bool found = false;
              for (auto& subCluster : subClusters) {
                if (isCloseToSubCluster(subCluster, x, y)) { // Placeholder function
                  subCluster.minX = std::min(subCluster.minX, static_cast<uint32_t>(x));
                  subCluster.minY = std::min(subCluster.minY, static_cast<uint32_t>(y));
                  subCluster.maxX = std::max(subCluster.maxX, static_cast<uint32_t>(x));
                  subCluster.maxY = std::max(subCluster.maxY, static_cast<uint32_t>(y));
                  subCluster.charge += adc;
                  found = true;
                  break;
                }
              }
              // If no sub-cluster matches, create a new one
              if (!found) {
                ClusterProperties newSubCluster;
                newSubCluster.minX = x;
                newSubCluster.minY = y;
                newSubCluster.maxX = x;
                newSubCluster.maxY = y;
                newSubCluster.charge = adc;
                subClusters.push_back(newSubCluster);
              }
            }
          }
        }

        return subClusters;
      }


  // Helper function to determine if a pixel is close to a sub-cluster
  ALPAKA_FN_ACC bool isCloseToSubCluster(const ClusterProperties& subCluster, uint16_t x, uint16_t y) const {
    return std::abs(static_cast<int>(subCluster.minX) - x) <= 1 &&
           std::abs(static_cast<int>(subCluster.minY) - y) <= 1;
  }

  // Helper function to get the cluster position (implement this based on your data)
  GlobalPoint getClusterPosition(uint32_t clusterIdx, const SiPixelDigisSoAConstView& digiView) const {
    // Placeholder: Assuming that cluster position is derived from the cluster data.
    return GlobalPoint(0.0, 0.0, 0.0);  // Replace with actual cluster position calculation
  }


};

};
*/

template <typename TrackerTraits>
void runKernels(TrackingRecHitSoAView<TrackerTraits>& hitView,
                SiPixelDigisSoAView& digiView,
                SiPixelClustersSoAView& clustersView,
                const CandidateGPUData* candidates,  // Updated to use CandidateGPUData
                size_t nCandidates,                  // Added number of candidates
                Queue& queue) {
    uint32_t items = 64;
    uint32_t groupsHits = divide_up_by(hitView.metadata().size(), items);
    uint32_t groupsDigis = divide_up_by(digiView.metadata().size(), items);
    uint32_t groupsClusters = divide_up_by(clustersView.metadata().size(), items);

    uint32_t groups = std::max({groupsHits, groupsDigis, groupsClusters});  // Ensure work division covers all three views

    auto workDiv = make_workdiv<Acc1D>(groups, items);

    // Kernel execution for hits, digis, clusters, and candidate data
    alpaka::exec<Acc1D>(queue, workDiv, Printout<TrackerTraits>{}, hitView, digiView, clustersView, candidates, nCandidates);
}

// Explicit template instantiation for Phase 1
template void runKernels<pixelTopology::Phase1>(TrackingRecHitSoAView<pixelTopology::Phase1>& hitView,
                                                SiPixelDigisSoAView& digiView,
                                                SiPixelClustersSoAView& clustersView,
                                                const CandidateGPUData* candidates,  // Updated to CandidateGPUData
                                                size_t nCandidates,                  // Added number of candidates
                                                Queue& queue);

// Explicit template instantiation for Phase 2
template void runKernels<pixelTopology::Phase2>(TrackingRecHitSoAView<pixelTopology::Phase2>& hitView,
                                                SiPixelDigisSoAView& digiView,
                                                SiPixelClustersSoAView& clustersView,
                                                const CandidateGPUData* candidates,  // Updated to CandidateGPUData
                                                size_t nCandidates,                  // Added number of candidates
                                                Queue& queue);
  }  // namespace Splitting
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
