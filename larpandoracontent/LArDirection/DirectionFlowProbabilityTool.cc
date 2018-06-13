/**
 *  @file   larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.cc
 *
 *  @brief  Implementation of the Svm vertex selection algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArDirection/DirectionFlowProbabilityTool.h"

using namespace pandora;

namespace lar_content
{
DirectionFlowProbabilityTool::DirectionFlowProbabilityTool() :
    m_slidingFitWindow(20),
    m_impactRadius(5.f),
    m_extrapolationNSteps(200),
    m_extrapolationStepSize(0.1f),
    m_minimumClusterLength(10.f),
    m_maxBestChiSquaredPerHit(3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DirectionFlowProbabilityTool::GetDirectionFlowProbability(std::function< TrackDirectionTool::DirectionFitObject(const pandora::Cluster* const pCluster) > &directionToolLambda, const pandora::CartesianVector &vertexPosition, const pandora::ClusterList &inputClusterList) const
{
    try
    {
        //Always interested only in W projection of vertex position
        const pandora::CartesianVector vertexProjection(LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition, TPC_VIEW_W));
        
        //Creates sliding fits and direction fits
        pandora::ClusterVector inputClusterVector;
        this->SelectClusters(directionToolLambda, inputClusterList, inputClusterVector);

        //Get primary clusters
        pandora::ClusterVector primaryClusters(this->GetPrimaryClusters(vertexProjection, inputClusterVector));

        //Check if there are any sizable clusters at all
        if (primaryClusters.size() == 0) //|| GetSmallestClusterLength(inputClusterVector) < 30.f
            return 0.5;
        
        float accumulatedProbability(1.f);
        GetDirectionFlowContribution(inputClusterVector, primaryClusters, vertexProjection, accumulatedProbability);
    
        return accumulatedProbability;
    }
    catch (...)
    {
        std::cout << "Something went wrong with direction flow." << std::endl;
        return 0.5f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int DirectionFlowProbabilityTool::GetNumberConsideredClusters(std::function< TrackDirectionTool::DirectionFitObject(const pandora::Cluster* const pCluster) > &directionToolLambda, const pandora::CartesianVector &vertexPosition, const pandora::ClusterList &inputClusterList) const
{
    try
    {
        //Always interested only in W projection of vertex position
        const pandora::CartesianVector vertexProjection(LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition, TPC_VIEW_W));
        
        //Creates sliding fits and direction fits
        pandora::ClusterVector inputClusterVector;
        this->SelectClusters(directionToolLambda, inputClusterList, inputClusterVector);

        //Get primary clusters
        pandora::ClusterVector primaryClusters(this->GetPrimaryClusters(vertexProjection, inputClusterVector));

        int nConsideredClusters(primaryClusters.size());

        for (const auto pPrimaryCluster : primaryClusters)
        {
            pandora::ClusterVector daughterClusters(this->GetOrderedDaughters(vertexProjection, pPrimaryCluster, inputClusterVector, primaryClusters));
            nConsideredClusters += daughterClusters.size();
        }
    
        return nConsideredClusters;
    }
    catch (...)
    {
        std::cout << "Something went wrong with direction flow." << std::endl;
        return 0;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFlowProbabilityTool::GetDirectionFlowContribution(pandora::ClusterVector &inputClusterVector, pandora::ClusterVector &primaryClusters, const pandora::CartesianVector &vertexProjection, float &accumulatedProbability) const
{
    for (const auto pPrimaryCluster : primaryClusters)
    {
        TrackDirectionTool::DirectionFitObject fitResult(this->GetCachedDirectionFit(pPrimaryCluster));

        //DRAW
        //std::cout << "Min chi squared per hit: " << fitResult.GetMinChiSquaredPerHit() << std::endl;
        //fitResult.DrawFit();
        //PANDORA_MONITORING_API(Pause(this->GetPandora()));

        if (fitResult.GetMinChiSquaredPerHit() < m_maxBestChiSquaredPerHit && (vertexProjection - fitResult.GetBeginpoint()).GetMagnitude() < (vertexProjection - fitResult.GetEndpoint()).GetMagnitude())
        {
            float fractionalCloseness(1.0 - ((fitResult.GetBeginpoint() - vertexProjection).GetMagnitude())/((fitResult.GetBeginpoint() - fitResult.GetEndpoint()).GetMagnitude()));
            accumulatedProbability *= fractionalCloseness * fitResult.GetProbability();
        }
        else if (fitResult.GetMinChiSquaredPerHit() < m_maxBestChiSquaredPerHit && (vertexProjection - fitResult.GetBeginpoint()).GetMagnitude() > (vertexProjection - fitResult.GetEndpoint()).GetMagnitude())
        {
            float fractionalCloseness(1.0 - ((fitResult.GetEndpoint() - vertexProjection).GetMagnitude())/((fitResult.GetBeginpoint() - fitResult.GetEndpoint()).GetMagnitude()));
            accumulatedProbability *= fractionalCloseness * (1.0 - fitResult.GetProbability());
        }
        else
        {
            accumulatedProbability *= 0.5;
        }
        
        pandora::ClusterVector daughterClusters(this->GetOrderedDaughters(vertexProjection, pPrimaryCluster, inputClusterVector, primaryClusters));
        DaughtersDirectionFlowContribution(daughterClusters, pPrimaryCluster, vertexProjection, accumulatedProbability);
        
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFlowProbabilityTool::DaughtersDirectionFlowContribution(pandora::ClusterVector &daughterClusters, const pandora::Cluster* const pPrimaryCluster, const pandora::CartesianVector &vertexProjection, float &accumulatedProbability) const
{
    TwoDSlidingFitResult parentSlidingFit(this->GetCachedSlidingFit(pPrimaryCluster));
    pandora::CartesianVector parentOrigin(vertexProjection);
    
    for (const auto pDaughterCluster : daughterClusters)
    {
        pandora::CartesianVector currentOrigin((parentSlidingFit.GetGlobalMinLayerPosition() - parentOrigin).GetMagnitude() > (parentSlidingFit.GetGlobalMaxLayerPosition() - parentOrigin).GetMagnitude() ? parentSlidingFit.GetGlobalMinLayerPosition() : parentSlidingFit.GetGlobalMaxLayerPosition()); //whichever sliding fit endpoint is furthest away from the previous origin
        
        TrackDirectionTool::DirectionFitObject daughterFitResult(this->GetCachedDirectionFit(pDaughterCluster));

        //DRAW
        //std::cout << "Min chi squared per hit: " << daughterFitResult.GetMinChiSquaredPerHit() << std::endl;
        //daughterFitResult.DrawFit();
        //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        
        if (daughterFitResult.GetMinChiSquaredPerHit() < m_maxBestChiSquaredPerHit && (currentOrigin - daughterFitResult.GetBeginpoint()).GetMagnitude() < (currentOrigin - daughterFitResult.GetEndpoint()).GetMagnitude())
        {
            float fractionalCloseness(1.0-((daughterFitResult.GetBeginpoint() - currentOrigin).GetMagnitude())/((daughterFitResult.GetBeginpoint() - daughterFitResult.GetEndpoint()).GetMagnitude()));
            accumulatedProbability *= fractionalCloseness * daughterFitResult.GetProbability();
        }
        else if (daughterFitResult.GetMinChiSquaredPerHit() < m_maxBestChiSquaredPerHit && (currentOrigin - daughterFitResult.GetBeginpoint()).GetMagnitude() > (currentOrigin - daughterFitResult.GetEndpoint()).GetMagnitude())
        {
            float fractionalCloseness(1.0-((daughterFitResult.GetEndpoint() - currentOrigin).GetMagnitude())/((daughterFitResult.GetBeginpoint() - daughterFitResult.GetEndpoint()).GetMagnitude()));
            accumulatedProbability *= fractionalCloseness * (1.0 - daughterFitResult.GetProbability());
        }
        else
        {
            accumulatedProbability *= 0.5;
        }
        
        parentOrigin = currentOrigin;
        parentSlidingFit = this->GetCachedSlidingFit(pDaughterCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DirectionFlowProbabilityTool::GetSmallestClusterLength(pandora::ClusterVector &inputClusterVector) const
{
    float longestClusterLength(0.f);

    for (const auto pCluster : inputClusterVector)
    {
        if (LArClusterHelper::GetLength(pCluster) > longestClusterLength)
            longestClusterLength = LArClusterHelper::GetLength(pCluster);
    }

    return longestClusterLength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFlowProbabilityTool::SelectClusters(std::function< TrackDirectionTool::DirectionFitObject(const pandora::Cluster* const pCluster) > &directionToolLambda, const pandora::ClusterList &clusterList, pandora::ClusterVector &selectedClusterVector) const
{
    ClusterVector sortedClusters(clusterList.begin(), clusterList.end());
    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);
    
    for (const Cluster *const pCluster : sortedClusters)
    {
        if (LArClusterHelper::GetLength(pCluster) < m_minimumClusterLength)
            continue;
        
        try
        {
            this->AddToSlidingFitCache(directionToolLambda, pCluster);
            if (this->GetCachedDirectionFit(pCluster).GetMinChiSquaredPerHit() < 4.0)
                selectedClusterVector.push_back(pCluster);
        }
        catch (StatusCodeException &statusCodeException)
        {
            std::cout << "Failure to cache." << std::endl;
            continue;
            //if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            //    throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::ClusterVector DirectionFlowProbabilityTool::GetPrimaryClusters(const pandora::CartesianVector &positionVector, const pandora::ClusterVector &inputClusterVector) const
{
    ClusterToSpacepointsMap clusterToSpacepointsMap(this->FillClusterToSpacepointsMap(inputClusterVector));
    ClusterVector primaryClusterVector;

    const pandora::ClusterList clusterList(inputClusterVector.begin(), inputClusterVector.end());
    const pandora::CartesianVector closestHitPosition(clusterList.size() != 0 ? LArClusterHelper::GetClosestPosition(positionVector, clusterList) : positionVector);
    
    for (const auto pCluster : inputClusterVector)
    {
        if (this->ClusterPointsToPosition(pCluster, positionVector, clusterToSpacepointsMap))
        {
            //Don't consider clusters that happen to point to vertex but are far away
            if ((this->GetCachedDirectionFit(pCluster).GetBeginpoint() - positionVector).GetMagnitude() <= 5 * m_impactRadius
                || (this->GetCachedDirectionFit(pCluster).GetEndpoint() - positionVector).GetMagnitude() <= 5 * m_impactRadius)
                primaryClusterVector.push_back(pCluster);
        }
    }
    
    return primaryClusterVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::ClusterVector DirectionFlowProbabilityTool::GetOrderedDaughters(const pandora::CartesianVector &positionVector, const pandora::Cluster* const pPrimaryCluster, const pandora::ClusterVector &inputClusterVector, pandora::ClusterVector &primaryClusters) const
{
    ClusterToSpacepointsMap clusterToSpacepointsMap(this->FillClusterToSpacepointsMap(inputClusterVector));
    
    const pandora::Cluster* pCurrentCluster(pPrimaryCluster);
    TwoDSlidingFitResult currentSlidingFit(this->GetCachedSlidingFit(pCurrentCluster));
    pandora::CartesianVector currentOrigin((currentSlidingFit.GetGlobalMinLayerPosition() - positionVector).GetMagnitude() > (currentSlidingFit.GetGlobalMaxLayerPosition() - positionVector).GetMagnitude() ? currentSlidingFit.GetGlobalMinLayerPosition() : currentSlidingFit.GetGlobalMaxLayerPosition()); //whichever sliding fit endpoint is furthest away from the previous origin
    
    pandora::ClusterVector orderedClusterVector;
    bool daughterFound(true);
    
    pandora::ClusterVector checkedClusters;
    checkedClusters.push_back(pPrimaryCluster);
    
    while (daughterFound)
    {
        daughterFound = false;
        
        for (const auto pCluster : inputClusterVector)
        {
            if (this->ClusterPointsToPosition(pCluster, currentOrigin, clusterToSpacepointsMap) && std::find(checkedClusters.begin(), checkedClusters.end(), pCluster) == checkedClusters.end())
            {
                //If daughter alredy considered or if cluster in primaries vector, disregard
                if (std::find(checkedClusters.begin(), checkedClusters.end(), pCluster) != checkedClusters.end()
                    || std::find(primaryClusters.begin(), primaryClusters.end(), pCluster) != primaryClusters.end())
                    continue;
                
                pCurrentCluster = pCluster;
                currentSlidingFit = this->GetCachedSlidingFit(pCluster);
                currentOrigin = ((currentSlidingFit.GetGlobalMinLayerPosition() - currentOrigin).GetMagnitude() > (currentSlidingFit.GetGlobalMaxLayerPosition() - currentOrigin).GetMagnitude() ? currentSlidingFit.GetGlobalMinLayerPosition() : currentSlidingFit.GetGlobalMaxLayerPosition()); //whichever sliding fit endpoint is furthest away from the previous origin
                
                orderedClusterVector.push_back(pCluster);
                checkedClusters.push_back(pCluster);
                
                daughterFound = true;
                break;
            }
        }
    }
    
    return orderedClusterVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionFlowProbabilityTool::ClusterPointsToPosition(const Cluster *const pCluster, const pandora::CartesianVector &positionVector, ClusterToSpacepointsMap &clusterToSpacepointsMap) const
{
    const auto clusterSpacePoints(clusterToSpacepointsMap.at(pCluster));
    
    for (const auto position : clusterSpacePoints)
    {
        if ((position - positionVector).GetMagnitude() <= m_impactRadius)
            return true;
    }
    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

DirectionFlowProbabilityTool::ClusterToSpacepointsMap DirectionFlowProbabilityTool::FillClusterToSpacepointsMap(const ClusterVector &clusterVector) const
{
    ClusterToSpacepointsMap clusterToSpacepointsMap;
    
    for (const Cluster *const pCluster : clusterVector)
    {
        ClusterToSpacepointsMap::iterator mapIter(clusterToSpacepointsMap.emplace(pCluster, CartesianPointVector()).first);
        this->GetSpacepoints(pCluster, mapIter->second);
    }
    
    return clusterToSpacepointsMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFlowProbabilityTool::GetSpacepoints(const Cluster *const pCluster, CartesianPointVector &spacepoints) const
{
    LArClusterHelper::GetCoordinateVector(pCluster, spacepoints);
    
    const TwoDSlidingFitResult &fitResult(this->GetCachedSlidingFit(pCluster));
    const float minLayerRL(fitResult.GetL(fitResult.GetMinLayer()));
    const float maxLayerRL(fitResult.GetL(fitResult.GetMaxLayer()));
    
    for (unsigned int iStep = 0; iStep < m_extrapolationNSteps; ++ iStep)
    {
        const float deltaRL(static_cast<float>(iStep) * m_extrapolationStepSize);
        
        CartesianVector positionPositive(0.f, 0.f, 0.f), positionNegative(0.f, 0.f, 0.f);
        fitResult.GetExtrapolatedPosition(maxLayerRL + deltaRL, positionPositive);
        fitResult.GetExtrapolatedPosition(minLayerRL - deltaRL, positionNegative);
        
        spacepoints.push_back(positionPositive);
        spacepoints.push_back(positionNegative);
    }
    
    std::sort(spacepoints.begin(), spacepoints.end(), LArClusterHelper::SortCoordinatesByPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFlowProbabilityTool::AddToSlidingFitCache(std::function< TrackDirectionTool::DirectionFitObject(const pandora::Cluster* const pCluster) > &directionToolLambda, const Cluster *const pCluster) const
{
    if (m_slidingFitResultMap.find(pCluster) != m_slidingFitResultMap.end())
        return;
    
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFit(pCluster, m_slidingFitWindow, slidingFitPitch);
    
    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, slidingFit)).second)
    {
        std::cout << "Sliding fit failure" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }
    
    try
    {
        TrackDirectionTool::DirectionFitObject directionFit(directionToolLambda(pCluster));
        m_directionFitMap.insert(std::make_pair(pCluster, directionFit));
    }
    catch (...)
    {
        std::cout << "Direction fit failure" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const TwoDSlidingFitResult &DirectionFlowProbabilityTool::GetCachedSlidingFit(const Cluster *const pCluster) const
{
    TwoDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);
    
    if (m_slidingFitResultMap.end() == iter)
    {
        std::cout << "Sliding fit retrieval failure" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
    
    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackDirectionTool::DirectionFitObject DirectionFlowProbabilityTool::GetCachedDirectionFit(const Cluster *const pCluster) const
{
    const auto iter = m_directionFitMap.find(pCluster);
    
    if (m_directionFitMap.end() == iter)
    {
        std::cout << "Direction fit retrieval failure" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
    
    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DirectionFlowProbabilityTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImpactRadius", m_impactRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ExtrapolationNumberSteps", m_extrapolationNSteps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ExtrapolationStepSize", m_extrapolationStepSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumClusterLength", m_minimumClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxBestChiSquaredPerHit", m_maxBestChiSquaredPerHit));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
