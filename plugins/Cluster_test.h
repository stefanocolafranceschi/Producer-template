#ifndef DataFormats_SiPixelClusterSoA_test_alpaka_Hits_test_h
#define DataFormats_SiPixelClusterSoA_test_alpaka_Hits_test_h

#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsSoA.h"
#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisSoA.h"
#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersSoA.h"
#include "DataFormats/VertexSoA/interface/ZVertexSoA.h"

#include "DataFormats/ClusterGeometrySoA/interface/ClusterGeometryLayout.h"
#include "DataFormats/ClusterGeometrySoA/interface/alpaka/ClusterGeometrySoACollection.h"

#include "DataFormats/CandidateSoA/interface/CandidateLayout.h"
#include "DataFormats/CandidateSoA/interface/alpaka/CandidateSoACollection.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

using namespace reco;

        // Cluster properties structure
        struct ClusterProperties {
            uint32_t minX = std::numeric_limits<uint32_t>::max();
            uint32_t minY = std::numeric_limits<uint32_t>::max();
            uint32_t maxX = 0;
            uint32_t maxY = 0;
            uint32_t charge = 0;

            // Fields to store error information after splitting
            float splitClusterErrorX = 0.0f;
            float splitClusterErrorY = 0.0f;

            // Device-compatible sizeX() and sizeY() methods
            ALPAKA_FN_ACC uint32_t sizeX() const {
                return maxX - minX + 1;
            }

            ALPAKA_FN_ACC uint32_t sizeY() const {
                return maxY - minY + 1;
            }

            ALPAKA_FN_ACC float getChargeAsFloat() const {
                return static_cast<float>(charge);
            }
        };

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
                  Queue& queue);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::Splitting

#endif  // DataFormats_SiPixelClusterSoA_test_alpaka_Hits_test_h
