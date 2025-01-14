#ifndef DataFormats_SiPixelClusterSoA_test_alpaka_Hits_test_h
#define DataFormats_SiPixelClusterSoA_test_alpaka_Hits_test_h

//#include <alpaka/alpaka.hpp>

#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsSoA.h"
#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisSoA.h"
#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersSoA.h"
#include "DataFormats/VertexSoA/interface/ZVertexSoA.h"

#include "DataFormats/ClusterGeometrySoA/interface/ClusterGeometryLayout.h"
#include "DataFormats/ClusterGeometrySoA/interface/alpaka/ClusterGeometrySoACollection.h"

#include "DataFormats/CandidateSoA/interface/CandidateLayout.h"
#include "DataFormats/CandidateSoA/interface/alpaka/CandidateSoACollection.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include <alpaka/alpaka.hpp>

constexpr int maxPixels = 200;
constexpr int maxSubClusters = 200;

// This represent a per-cluster data needed in the Splitting algorithm
struct clusterProperties {
    float clx[maxSubClusters];
    float cly[maxSubClusters];
    float cls[maxSubClusters];
    float oldclx[maxSubClusters];
    float oldcly[maxSubClusters];
    uint32_t pixelCounter;                // how many pixels in this Cluster

    float distanceMap[maxPixels][maxSubClusters];
    int scoresIndices[maxPixels];
    float scoresValues[maxPixels];

    float clusterForPixel[maxPixels];
    float pixels[maxPixels];              // Storing the index of the pixel
    float weightOfPixel[maxPixels];
};

using namespace reco;

namespace ALPAKA_ACCELERATOR_NAMESPACE::Splitting {

  template <typename TrackerTraits>
  void runKernels(TrackingRecHitSoAView<TrackerTraits>& hits,
                  SiPixelDigisSoAView& digis,
                  SiPixelClustersSoAView& clusters,
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
                  double chargePerUnit_,
                  double fractionalWidth_,
                  SiPixelDigisSoAView& outputDigis,
                  SiPixelClustersSoAView& outputClusters,
                  clusterProperties* clusterPropertiesDevice,
                  Queue& queue);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::Splitting

#endif  // DataFormats_SiPixelClusterSoA_test_alpaka_Hits_test_h
