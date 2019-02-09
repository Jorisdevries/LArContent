/**
 *  @file   larpandoracontent/LArHelpers/LArDirectionHelper.cc
 *
 *  @brief  Implementation of the file helper class.
 *
 *  $Log: $
 */

#include "Helpers/XmlHelper.h"

#include "larpandoracontent/LArHelpers/LArDirectionHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include <cstdlib>
#include <sys/stat.h>

using namespace pandora;

namespace lar_content

//------------------------------------------------------------------------------------------------------------------------------------------
{

float LArDirectionHelper::GetAngleWithVector(const Pandora &pandora, const pandora::ParticleFlowObject* pPfo, CartesianVector &axisVector)
{
    LArTrackStateVector trackStateVector; 
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, LArPfoHelper::GetVertex(pPfo), 20, LArGeometryHelper::GetWireZPitch(pandora), trackStateVector);

    if (trackStateVector.size() == 0)
        return -1.f;

    //ascending order: front is low y
    std::sort(trackStateVector.begin(), trackStateVector.end(), [](LArTrackState &state1, LArTrackState &state2){return state1.GetPosition().GetY() < state2.GetPosition().GetY();});

    pandora::CartesianVector vertexDirection(trackStateVector.back().GetPosition() - trackStateVector.front().GetPosition());

    //if (trackStateVector.back().GetZ() > trackStateVector.front().GetZ())
    //    vertexDirection = (trackStateVector.back().GetPosition() - trackStateVector.front().GetPosition());

    //if (trackStateVector.back().GetPosition().GetZ() < trackStateVector.front().GetPosition().GetZ())
    //    vertexDirection = (trackStateVector.front().GetPosition() - trackStateVector.back().GetPosition());

    return axisVector.GetOpeningAngle(vertexDirection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<CartesianVector> LArDirectionHelper::GetLowHighYPoints(const Pandora &pandora, const pandora::ParticleFlowObject* pPfo)
{
    LArTrackStateVector trackStateVector; 
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, LArPfoHelper::GetVertex(pPfo), 20, LArGeometryHelper::GetWireZPitch(pandora), trackStateVector);

    std::vector<CartesianVector> positions;

    if (trackStateVector.size() == 0)
        return positions;

    //ascending order: front is low y
    std::sort(trackStateVector.begin(), trackStateVector.end(), [](LArTrackState &state1, LArTrackState &state2){return state1.GetPosition().GetY() < state2.GetPosition().GetY();});

    pandora::CartesianVector initialPosition(trackStateVector.front().GetPosition());
    pandora::CartesianVector endPosition(trackStateVector.back().GetPosition());

    pandora::CartesianVector lowYVector(initialPosition.GetY() < endPosition.GetY() ? initialPosition : endPosition);
    pandora::CartesianVector highYVector(initialPosition.GetY() > endPosition.GetY() ? initialPosition : endPosition);

    positions.push_back(lowYVector);
    positions.push_back(highYVector);

    return positions; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<CartesianVector> LArDirectionHelper::GetLowHighZPoints(const pandora::Cluster* pCluster)
{
    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);
    CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());

    std::vector<CartesianVector> positions;

    if (caloHitVector.size() == 0)
        return positions;

    //ascending order: front is low z
    std::sort(caloHitVector.begin(), caloHitVector.end(), [](const pandora::CaloHit* const pCaloHit1, const pandora::CaloHit* const pCaloHit2){return pCaloHit1->GetPositionVector().GetZ() < pCaloHit2->GetPositionVector().GetZ();});

    pandora::CartesianVector initialPosition(caloHitVector.front()->GetPositionVector());
    pandora::CartesianVector endPosition(caloHitVector.back()->GetPositionVector());

    pandora::CartesianVector lowZVector(initialPosition.GetZ() < endPosition.GetZ() ? initialPosition : endPosition);
    pandora::CartesianVector highZVector(initialPosition.GetZ() > endPosition.GetZ() ? initialPosition : endPosition);

    positions.push_back(lowZVector);
    positions.push_back(highZVector);

    return positions; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArDirectionHelper::IntersectsYFace(const Pandora &pandora, const pandora::ParticleFlowObject* pPfo)
{
    std::vector<CartesianVector> positions(GetLowHighYPoints(pandora, pPfo));

    if (positions.size() == 0)
        return false;

    pandora::CartesianVector lowYVector(positions.front());
    pandora::CartesianVector highYVector(positions.back());

    float xExtent(highYVector.GetX() - lowYVector.GetX()), yExtent(highYVector.GetY() - lowYVector.GetY()), zExtent(highYVector.GetZ() - lowYVector.GetZ());
    float xSlope(xExtent/yExtent), zSlope(zExtent/yExtent);
    float yDistanceToTravel(116.5 - highYVector.GetY());
    CartesianVector yFaceIntersection(highYVector.GetX() + xSlope*yDistanceToTravel, 116.5, highYVector.GetZ() + zSlope*yDistanceToTravel);

    if (yFaceIntersection.GetX() > 0.0 && yFaceIntersection.GetX() < 256.35 && yFaceIntersection.GetZ() > 0.0 && yFaceIntersection.GetZ() < 1036.8)
        return true;
    else 
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArDirectionHelper::IntersectsYFace(const pandora::CartesianVector &lowYVector, const pandora::CartesianVector &highYVector)
{
    float xExtent(highYVector.GetX() - lowYVector.GetX()), yExtent(highYVector.GetY() - lowYVector.GetY()), zExtent(highYVector.GetZ() - lowYVector.GetZ());
    float xSlope(xExtent/yExtent), zSlope(zExtent/yExtent);
    float yDistanceToTravel(116.5 - highYVector.GetY());
    CartesianVector yFaceIntersection(highYVector.GetX() + xSlope*yDistanceToTravel, 116.5, highYVector.GetZ() + zSlope*yDistanceToTravel);

    if (yFaceIntersection.GetX() > 0.0 && yFaceIntersection.GetX() < 256.35 && yFaceIntersection.GetZ() > 0.0 && yFaceIntersection.GetZ() < 1036.8)
        return true;
    else 
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArDirectionHelper::IntersectsYFace(TrackDirectionTool::DirectionFitObject &fitResult) 
{
    const pandora::CartesianVector initialPosition(fitResult.GetBeginpoint());
    const pandora::CartesianVector endPosition(fitResult.GetEndpoint());

    const pandora::CartesianVector lowYVector(initialPosition.GetY() < endPosition.GetY() ? initialPosition : endPosition);
    const pandora::CartesianVector highYVector(initialPosition.GetY() > endPosition.GetY() ? initialPosition : endPosition);

    float xExtent(highYVector.GetX() - lowYVector.GetX()), yExtent(highYVector.GetY() - lowYVector.GetY()), zExtent(highYVector.GetZ() - lowYVector.GetZ());
    float xSlope(xExtent/yExtent), zSlope(zExtent/yExtent);
    float yDistanceToTravel(116.5 - highYVector.GetY());
    CartesianVector yFaceIntersection(highYVector.GetX() + xSlope*yDistanceToTravel, 116.5, highYVector.GetZ() + zSlope*yDistanceToTravel);

    if (yFaceIntersection.GetX() > 0.0 && yFaceIntersection.GetX() < 256.35 && yFaceIntersection.GetZ() > 0.0 && yFaceIntersection.GetZ() < 1036.8)
        return true;
    else 
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArDirectionHelper::HasFiducialLowY(const Pandora &pandora, const pandora::ParticleFlowObject* pPfo)
{
    std::vector<CartesianVector> positions(GetLowHighYPoints(pandora, pPfo));

    if (positions.size() == 0)
        return false;

    pandora::CartesianVector lowYVector(positions.front());
    pandora::CartesianVector highYVector(positions.back());

    if (IsInFiducialVolume(lowYVector))
        return true;
    else 
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArDirectionHelper::HasHighTopY(const Pandora &pandora, const pandora::ParticleFlowObject* pPfo, float threshold)
{
    std::vector<CartesianVector> positions(GetLowHighYPoints(pandora, pPfo));

    if (positions.size() == 0)
        return false;

    pandora::CartesianVector lowYVector(positions.front());
    pandora::CartesianVector highYVector(positions.back());

    if (highYVector.GetY() > threshold)
        return true;
    else 
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArDirectionHelper::CalculateCosmicProbability(TrackDirectionTool::DirectionFitObject &directionFit)
{
    float upDownDeltaChiSquaredPerHit(std::abs(directionFit.GetUpDownDeltaChiSquaredPerHit())), minChiSquaredPerHit(directionFit.GetMinChiSquaredPerHit()); 
    int numberHits(directionFit.GetNHits());

    float p_max_par1(9.57688e-01), p_max_par2(2.07682e-04), p_max_par3(-4.73225e-02);
    float alpha_par1(4.89049e-01), alpha_par2(4.91588e-03), alpha_par3(-3.68713e-01);
    float beta_par1(1.62781e-02), beta_par2(-5.80356e-05), beta_par3(1.82657e-03);

    float p_max(std::max(p_max_par1 + p_max_par3 * minChiSquaredPerHit + p_max_par2 * numberHits, 0.5f)), alpha(std::max(alpha_par1 + alpha_par3 * minChiSquaredPerHit + alpha_par2 * numberHits, 0.1f)), beta(std::max(beta_par1 + beta_par3 * minChiSquaredPerHit + beta_par2 * numberHits, 0.001f)); 

    float p0(0.5 + ((1/(2*alpha)) * ( (2*alpha*p_max + 2*beta*p_max - alpha - beta) * std::pow(((alpha + beta)/beta), (beta/alpha)))));

    if (p0 > 1.0) 
        p0 = 1.0; 

    if (alpha > 1.6) 
        alpha = 1.6; 

    float probability(0.5 + (p0 - 0.5) * (1 - exp(-alpha*upDownDeltaChiSquaredPerHit)) * exp(-beta*upDownDeltaChiSquaredPerHit));

    if (!IntersectsYFace(directionFit))
        probability = 0.0; 

    return probability;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArDirectionHelper::IsStoppingTopFaceMCParticle(const pandora::MCParticle* pMCParticle)
{
    const pandora::CartesianVector initialPosition(pMCParticle->GetVertex());
    const pandora::CartesianVector endPosition(pMCParticle->GetEndpoint());

    pandora::CartesianVector lowYVector(initialPosition.GetY() < endPosition.GetY() ? initialPosition : endPosition);
    pandora::CartesianVector highYVector(initialPosition.GetY() > endPosition.GetY() ? initialPosition : endPosition);

    float xExtent(highYVector.GetX() - lowYVector.GetX()), yExtent(highYVector.GetY() - lowYVector.GetY()), zExtent(highYVector.GetZ() - lowYVector.GetZ());
    float xSlope(xExtent/yExtent), zSlope(zExtent/yExtent);
    float yDistanceToTravel(116.5 - lowYVector.GetY());
    CartesianVector yFaceIntersection(lowYVector.GetX() + xSlope*yDistanceToTravel, 116.5, lowYVector.GetZ() + zSlope*yDistanceToTravel);

    if (IsInFiducialVolume(lowYVector) && highYVector.GetZ() > 116.5 && yFaceIntersection.GetX() > 0.0 && yFaceIntersection.GetX() < 256.35 && yFaceIntersection.GetZ() > 0.0 && yFaceIntersection.GetZ() < 1036.8) 
        return true;
    else 
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArDirectionHelper::IsInFiducialVolume(pandora::CartesianVector positionVector)
{
    if ((positionVector.GetX() > 12.0 && positionVector.GetX() < (256.35 - 12.0)) && (positionVector.GetY() > (-116.5 + 35.0) && positionVector.GetY() < (116.5 - 35.0)) && ((positionVector.GetZ() > 25.0 && positionVector.GetZ() < 675.0) || (positionVector.GetZ() > 775.0 && positionVector.GetZ() < (1036.8 - 85.0))))
        return true;
    else 
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
