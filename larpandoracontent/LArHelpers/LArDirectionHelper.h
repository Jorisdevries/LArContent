/**
 *  @file   larpandoracontent/LArHelpers/LArDirectionHelper.h
 *
 *  @brief  Header file for the file helper class
 *
 *  $Log: $
 */
#ifndef LAR_DIRECTION_HELPER_H
#define LAR_DIRECTION_HELPER_H 1

#include <string>

#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/Vertex.h"
#include "Objects/MCParticle.h"

#include "larpandoracontent/LArDirection/TrackDirectionTool.h"

namespace lar_content
{

/**
 *  @brief  LArDirectionHelper class
 */
class LArDirectionHelper
{
public:

    static float GetAngleWithVector(const pandora::Pandora &pandora, const pandora::ParticleFlowObject* pPfo, pandora::CartesianVector &axisVector);

    static std::vector<pandora::CartesianVector> GetLowHighYPoints(const pandora::Pandora &pandora, const pandora::ParticleFlowObject* pPfo);

    static std::vector<pandora::CartesianVector> GetLowHighZPoints(const pandora::Cluster* pCluster);

    static bool IntersectsYFace(const pandora::Pandora &pandora, const pandora::ParticleFlowObject* pPfo);

    static bool IntersectsYFace(const pandora::CartesianVector &lowYVector, const pandora::CartesianVector &highYVector);

    static bool IntersectsYFace(TrackDirectionTool::DirectionFitObject &fitResult);

    static bool HasFiducialLowY(const pandora::Pandora &pandora, const pandora::ParticleFlowObject* pPfo);

    static bool HasHighTopY(const pandora::Pandora &pandora, const pandora::ParticleFlowObject* pPfo, float threshold);

    static float CalculateCosmicProbability(TrackDirectionTool::DirectionFitObject &directionFit);

    static bool IsStoppingTopFaceMCParticle(const pandora::MCParticle* pMCParticle);

    static bool IsInFiducialVolume(pandora::CartesianVector positionVector);

    static pandora::CartesianVector GetApproximateNeutrinoMomentum(const pandora::Pandora &pandora, pandora::PfoList pfoList, const pandora::ParticleFlowObject* pLongestPfo);

    static pandora::CartesianVector GetApproximatePfoMomentum(const pandora::Pandora &pandora, const pandora::ParticleFlowObject* pPfo, const float &particleMass);

    static float GetPfoCharge(const pandora::ParticleFlowObject* pPfo); 

    static void GetClusterHits(pandora::PfoList &pfoList, pandora::CaloHitList &caloHitList); 
};

} // namespace lar_content

#endif // #ifndef LAR_DIRECTION_HELPER_H
