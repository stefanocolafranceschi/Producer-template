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

#include "DataFormats/VertexSoA/interface/ZVertexSoA.h"
#include "DataFormats/VertexSoA/interface/ZVertexHost.h"
#include "DataFormats/VertexSoA/interface/ZVertexDevice.h"
#include "DataFormats/VertexSoA/interface/alpaka/ZVertexSoACollection.h"

#include "DataFormats/GeometryVector/interface/Basic3DVector.h"
#include "DataFormats/Math/interface/SSEVec.h"
#include "DataFormats/Math/interface/ExtVec.h"

#include "DataFormats/ClusterGeometrySoA/interface/ClusterGeometryLayout.h"
#include "DataFormats/ClusterGeometrySoA/interface/alpaka/ClusterGeometrySoACollection.h"

#include "DataFormats/CandidateSoA/interface/CandidateLayout.h"
#include "DataFormats/CandidateSoA/interface/alpaka/CandidateSoACollection.h"

#include "Cluster_test.h"

using namespace alpaka;
using namespace reco;

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
                                ZVertexSoAView vertexView,
                                CandidateSoAView candidates,
                                ClusterGeometrySoAView geoclusters) const {         
 

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


        // Print debug info for Candidates -----------------------------------
        if (cms::alpakatools::once_per_grid(acc)) {
            printf("Candidate Info:\n");
            printf("nCandidates = %d\n", candidates.metadata().size());
        }
        // Iterate over the candidates (assuming candidates are indexed from 0 to nCandidates-1)
        for (uint32_t c : cms::alpakatools::uniform_elements(acc, candidates.metadata().size())) {
            printf("Candidate %d -> px: %.2f, py: %.2f, pz: %.2f, pt: %.2f, eta: %.2f, phi: %.2f\n",
                   c,
                   candidates[c].px(),
                   candidates[c].py(),
                   candidates[c].pz(),
                   candidates[c].pt(),
                   candidates[c].eta(),
                   candidates[c].phi());
        }


        // Print debug info for Vertices -----------------------------------
        if (cms::alpakatools::once_per_grid(acc)) {
          printf("Vertex Info:\n");
          printf("nVertices = %d\n", vertexView.metadata().size());
        }
        for (uint32_t v : cms::alpakatools::uniform_elements(acc, 10)) {
          printf("Vertex %d -> z: %.2f, w: %.2f, chi2: %.2f, pt^2: %.2f, sortedIndex: %d\n",
                 v,
                 vertexView[v].zv(),
                 vertexView[v].wv(),
                 vertexView[v].chi2(),
                 vertexView[v].ptv2(),
                 vertexView[v].sortInd());
        }


        // Print debug info for ClusterGeometry -----------------------------------
        if (cms::alpakatools::once_per_grid(acc)) {
          printf("ClusterGeometry Info:\n");
          printf("nClusterGeometries = %d\n", geoclusters.metadata().size());
        }
        for (uint32_t g : cms::alpakatools::uniform_elements(acc, 10)) {
          printf("geoclusters %d -> clusterId: %d, pitchX: %.2f, pitchY: %.2f, thickness: %.2f, tanLorentzAngle: %.2f\n",
                 g,
                 geoclusters[g].clusterIds(),
                 geoclusters[g].pitchX(),
                 geoclusters[g].pitchY(),
                 geoclusters[g].thickness(),
                 geoclusters[g].tanLorentzAngles());
        }



      }
    };





    template <typename TrackerTraits>
    struct JetSplit {


        // Main operator function
        template <typename TAcc, typename = std::enable_if_t<isAccelerator<TAcc>>>
        ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                      TrackingRecHitSoAConstView<TrackerTraits> hitView,
                                      SiPixelDigisSoAConstView digiView,
                                      SiPixelClustersSoAConstView clustersView,
                                      ZVertexSoAView vertexView,
                                      CandidateSoAView candidates,
                                      ClusterGeometrySoAView geoclusters,
                                      double ptMin_,
                                      double deltaR_,
                                      double chargeFracMin_,
                                      float expSizeXAtLorentzAngleIncidence_,
                                      float expSizeXDeltaPerTanAlpha_,
                                      float expSizeYAtNormalIncidence_,
                                      double centralMIPCharge_) const {
            // Iterate over clusters
            for (uint32_t clusterIdx : cms::alpakatools::uniform_elements(acc, clustersView.metadata().size())) {
                ClusterProperties clusterProps = computeClusterProperties(acc, digiView, clusterIdx);

                // Fetch the cluster's position and geometry
                auto clusterPosition = getClusterPosition(clusterIdx, hitView, clustersView);
                float pitchX = geoclusters[clusterIdx].pitchX();
                float pitchY = geoclusters[clusterIdx].pitchY();
                float thickness = geoclusters[clusterIdx].thickness();

                // Loop through all candidates (jets)
                for (uint32_t candIdx : cms::alpakatools::uniform_elements(acc, candidates.metadata().size())) {
                    const auto& jet = candidates[candIdx];

                    // Skip low-pt jets
                    if (jet.pt() < ptMin_)
                        continue;

                    // Extract jet direction components
                    float jetPx = jet.px();
                    float jetPy = jet.py();
                    float jetPz = jet.pz();

                    // Calculate deltaR directly with scalar values
                    float deltaEta = clusterPosition.y - jetPy;
                    float deltaPhi = clusterPosition.x - jetPx;
                    float deltaR = sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);

                    // Check deltaR condition and split clusters if applicable
                    if (deltaR < deltaR_ && shouldSplit(clusterProps, jetPx, jetPy, jetPz, pitchX, pitchY, thickness,
                                                        geoclusters[clusterIdx].tanLorentzAngles(),
                                                        chargeFracMin_,
                                                        expSizeXAtLorentzAngleIncidence_,
                                                        expSizeXDeltaPerTanAlpha_,
                                                        expSizeYAtNormalIncidence_,
                                                        centralMIPCharge_)) {
                        // Using a fixed-size array for sub-clusters
                        const int maxSubClusters = 10;          // Max limit
                        ClusterProperties subClusters[maxSubClusters];

                        int subClusterCount = 0;

                        // Split the cluster into sub-clusters
                        splitCluster(acc, clusterProps, candidates, digiView, clusterIdx, subClusters, subClusterCount);

                        // Print debug info for sub-clusters
                        for (int i = 0; i < subClusterCount; ++i) {
                            const auto& subCluster = subClusters[i];
                            printf("Sub-cluster -> minX: %d, minY: %d, maxX: %d, maxY: %d, charge: %d\n",
                                   subCluster.minX, subCluster.minY, subCluster.maxX, subCluster.maxY, subCluster.charge);
                        }
                    }
                }
            }
        }

        // Computes cluster properties
        template <typename TAcc>
        ALPAKA_FN_ACC ClusterProperties computeClusterProperties(TAcc const& acc,
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

        // Determines if a cluster should be split
        ALPAKA_FN_ACC bool shouldSplit(const ClusterProperties& cluster,
                                       float jetPx,
                                       float jetPy,
                                       float jetPz,
                                       float pitchX,
                                       float pitchY,
                                       float thickness,
                                       float tanLorentzAngles,
                                       double chargeFracMin_,
                                       float expSizeXAtLorentzAngleIncidence_,
                                       float expSizeXDeltaPerTanAlpha_,
                                       float expSizeYAtNormalIncidence_,
                                       double centralMIPCharge_) const {

            float jetTanAlpha = jetPx / jetPz;
            float jetTanBeta = jetPy / jetPz;
            float jetZOverRho = std::sqrt(jetTanAlpha * jetTanAlpha + jetTanBeta * jetTanBeta);

            float expSizeX = expSizeXAtLorentzAngleIncidence_ +
                             std::abs(expSizeXDeltaPerTanAlpha_ * (jetTanAlpha - tanLorentzAngles));
            float expSizeY = std::sqrt((expSizeYAtNormalIncidence_ * expSizeYAtNormalIncidence_) +
                                       thickness * thickness / (pitchY * pitchY) * jetTanBeta * jetTanBeta);

            if (expSizeX < 1.f)
                expSizeX = 1.f;
            if (expSizeY < 1.f)
                expSizeY = 1.f;

            float expCharge = std::sqrt(1.08f + jetZOverRho * jetZOverRho) * centralMIPCharge_;
            return cluster.charge > expCharge * chargeFracMin_ &&
                   (cluster.sizeX() > expSizeX + 1 || cluster.sizeY() > expSizeY + 1);
        }

        // Splits clusters into sub-clusters and stores them in a fixed-size array
        template <typename TAcc>
        ALPAKA_FN_ACC void splitCluster(TAcc const& acc,
                                         const ClusterProperties& clusterProps,
                                         const CandidateSoAView& candidates,
                                         const SiPixelDigisSoAConstView& digiView,
                                         uint32_t clusterIdx,
                                         ClusterProperties* subClusters,
                                         int& subClusterCount) const {
            for (uint32_t j : cms::alpakatools::uniform_elements(acc, digiView.metadata().size())) {
                if (static_cast<uint32_t>(digiView[j].clus()) == clusterIdx) {
                    uint16_t x = digiView[j].xx();
                    uint16_t y = digiView[j].yy();
                    uint16_t adc = digiView[j].adc();

                    for (int candIdx = 0; candIdx < static_cast<int>(candidates.metadata().size()); ++candIdx) {
                        const auto jetDir = candidates[candIdx];

                        float dx = x - static_cast<float>(jetDir.px());
                        float dy = y - static_cast<float>(jetDir.py());
                        float distance = std::sqrt(dx * dx + dy * dy);

                        if (distance < 2.0f) {
                            bool found = false;
                            for (int i = 0; i < subClusterCount; ++i) {
                                if (isCloseToSubCluster(subClusters[i], x, y)) {
                                    subClusters[i].charge += adc;
                                    subClusters[i].minX = std::min(subClusters[i].minX, static_cast<uint32_t>(x));
                                    subClusters[i].maxX = std::max(subClusters[i].maxX, static_cast<uint32_t>(x));
                                    subClusters[i].minY = std::min(subClusters[i].minY, static_cast<uint32_t>(y));
                                    subClusters[i].maxY = std::max(subClusters[i].maxY, static_cast<uint32_t>(y));
                                    found = true;
                                    break;
                                }
                            }

                            if (!found && subClusterCount < 10) { // Ensure we don't exceed the max sub-clusters
                                subClusters[subClusterCount].minX = x;
                                subClusters[subClusterCount].maxX = x;
                                subClusters[subClusterCount].minY = y;
                                subClusters[subClusterCount].maxY = y;
                                subClusters[subClusterCount].charge = adc;
                                ++subClusterCount;
                            }
                        }
                    }
                }
            }
        }

        // Checks if a pixel is close to an existing sub-cluster
        ALPAKA_FN_ACC bool isCloseToSubCluster(const ClusterProperties& subCluster, uint16_t x, uint16_t y) const {
            return (x >= subCluster.minX - 1 && x <= subCluster.maxX + 1 &&
                    y >= subCluster.minY - 1 && y <= subCluster.maxY + 1);
        }

        struct Position {
            float x, y, z;
        };

        ALPAKA_FN_ACC Position getClusterPosition(uint32_t clusterIdx,
                                                  TrackingRecHitSoAConstView<TrackerTraits> hitView,
                                                  SiPixelClustersSoAConstView clustersView) const {
            Position pos = {0.0f, 0.0f, 0.0f};      // Default to (0,0,0)

            if (clusterIdx < static_cast<uint32_t>(hitView.metadata().size())) {
                pos.x = hitView[clusterIdx].xGlobal();
                pos.y = hitView[clusterIdx].yGlobal();
                pos.z = hitView[clusterIdx].zGlobal();
            }
            return pos;
        }
    };



    template <typename TrackerTraits>
    void runKernels(TrackingRecHitSoAView<TrackerTraits>& hitView,
                    SiPixelDigisSoAView& digiView,
                    SiPixelClustersSoAView& clustersView,
                    ZVertexSoAView& vertexView,
                    CandidateSoAView& candidates,
                    ClusterGeometrySoAView& geoclusters,
                    double ptMin_,
                    double deltaR_,
                    double chargeFracMin_,
                    float expSizeXAtLorentzAngleIncidence_,
                    float expSizeXDeltaPerTanAlpha_,
                    float expSizeYAtNormalIncidence_,
                    double centralMIPCharge_,
                    Queue& queue) {
        uint32_t items = 64;
        uint32_t groupsHits = divide_up_by(hitView.metadata().size(), items);
        uint32_t groupsDigis = divide_up_by(digiView.metadata().size(), items);
        uint32_t groupsClusters = divide_up_by(clustersView.metadata().size(), items);

        uint32_t groups = std::max({groupsHits, groupsDigis, groupsClusters});

        auto workDiv = make_workdiv<Acc1D>(groups, items);

        // Kernel execution
        alpaka::exec<Acc1D>(queue, workDiv, Printout<TrackerTraits>{}, hitView, digiView, clustersView, vertexView, candidates, geoclusters);
        alpaka::exec<Acc1D>(queue, 
                            workDiv, 
                            JetSplit<TrackerTraits>{}, 
                            hitView, 
                            digiView, 
                            clustersView, 
                            vertexView, 
                            candidates, 
                            geoclusters,
                            ptMin_,
                            deltaR_,
                            chargeFracMin_,
                            expSizeXAtLorentzAngleIncidence_,
                            expSizeXDeltaPerTanAlpha_,
                            expSizeYAtNormalIncidence_,
                            centralMIPCharge_);
    }

    // Explicit template instantiation for Phase 1
    template void runKernels<pixelTopology::Phase1>(TrackingRecHitSoAView<pixelTopology::Phase1>& hitView,
                                                    SiPixelDigisSoAView& digiView,
                                                    SiPixelClustersSoAView& clustersView,
                                                    ZVertexSoAView& vertexView,
                                                    CandidateSoAView& candidates,
                                                    ClusterGeometrySoAView& geoclusters,
                                                    double ptMin_,
                                                    double deltaR_,
                                                    double chargeFracMin_,
                                                    float expSizeXAtLorentzAngleIncidence_,
                                                    float expSizeXDeltaPerTanAlpha_,
                                                    float expSizeYAtNormalIncidence_,
                                                    double centralMIPCharge_,
                                                    Queue& queue);

    // Explicit template instantiation for Phase 2
    template void runKernels<pixelTopology::Phase2>(TrackingRecHitSoAView<pixelTopology::Phase2>& hitView,
                                                    SiPixelDigisSoAView& digiView,
                                                    SiPixelClustersSoAView& clustersView,
                                                    ZVertexSoAView& vertexView,
                                                    CandidateSoAView& candidates,
                                                    ClusterGeometrySoAView& geoclusters,
                                                    double ptMin_,
                                                    double deltaR_,
                                                    double chargeFracMin_,
                                                    float expSizeXAtLorentzAngleIncidence_,
                                                    float expSizeXDeltaPerTanAlpha_,
                                                    float expSizeYAtNormalIncidence_,
                                                    double centralMIPCharge_,                                                    
                                                    Queue& queue);
  }  // namespace Splitting
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
