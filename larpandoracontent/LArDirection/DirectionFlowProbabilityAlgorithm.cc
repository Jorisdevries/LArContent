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

#include "larpandoracontent/LArDirection/DirectionFlowProbabilityAlgorithm.h"

using namespace pandora;

namespace lar_content
{
DirectionFlowProbabilityAlgorithm::DirectionFlowProbabilityAlgorithm() :
    m_slidingFitWindow(20),
    m_impactRadius(2.f),
    m_extrapolationNSteps(200),
    m_extrapolationStepSize(0.1f),
    m_minimumClusterLength(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFlowProbabilityAlgorithm::Test() 
{
    std::cout << "Using DirectionFlowProbabilityAlgorithm." << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DirectionFlowProbabilityAlgorithm::GetDirectionFlowProbability(const pandora::CartesianVector &positionVector, const pandora::ClusterList *pClusterList) 
{
    pandora::ClusterVector inputClusterVector;
    this->SelectClusters(pClusterList, inputClusterVector);

    if (inputClusterVector.size() == 0)
        throw STATUS_CODE_FAILURE;

    float accumulatedProbability(1.f);
    EmergingClusterVector finalDaughters(this->CreateEmergingClusters(positionVector, inputClusterVector));

    for (const auto pEmergingCluster : finalDaughters)
    {
        while (pEmergingCluster->ParentEmergingCluster() != nullptr)
        {
            const EmergingCluster* const pParentEmergingCluster(pEmergingCluster->ParentEmergingCluster());

            if (this->BeginpointIsParentOrigin(pEmergingCluster, pParentEmergingCluster))
                accumulatedProbability *= pEmergingCluster->ParentEmergingCluster()->DirectionFit()->GetProbability();
            else
                accumulatedProbability *= (1.f - pEmergingCluster->ParentEmergingCluster()->DirectionFit()->GetProbability());
        }
    }

    return accumulatedProbability;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFlowProbabilityAlgorithm::SelectClusters(const pandora::ClusterList *pClusterList, pandora::ClusterVector &selectedClusterVector) 
{
    ClusterVector sortedClusters(pClusterList->begin(), pClusterList->end());
    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : sortedClusters)
    {   
        if (LArClusterHelper::GetLength(pCluster) < m_minimumClusterLength)
            continue;

        try 
        {   
            this->AddToSlidingFitCache(pCluster);
            selectedClusterVector.push_back(pCluster);
        }   
        catch (StatusCodeException &statusCodeException)
        {   
            continue;
            //if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            //    throw statusCodeException;
        }  
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionFlowProbabilityAlgorithm::BeginpointIsParentOrigin(const EmergingCluster* const pEmergingCluster, const EmergingCluster* const pParentEmergingCluster) const
{
    if ((pParentEmergingCluster->DirectionFit()->GetBeginpoint() - *(pEmergingCluster->Origin())).GetMagnitude() < (pParentEmergingCluster->DirectionFit()->GetEndpoint() - *(pEmergingCluster->Origin())).GetMagnitude())
        return true;
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

DirectionFlowProbabilityAlgorithm::EmergingClusterVector DirectionFlowProbabilityAlgorithm::CreateEmergingClusters(const pandora::CartesianVector &vertexPosition, const pandora::ClusterVector &inputClusterVector) 
{
    //Creates all EmergingClusters in memory but returns only the last daughter for each primary (so we can walk 'up' each branch) 
    pandora::ClusterVector primaryClusters(this->GetPrimaryClusters(vertexPosition, inputClusterVector));
    EmergingClusterVector finalDaughters;

    for (const auto pPrimaryCluster : primaryClusters)
    {
        TrackDirectionTool::DirectionFitObject* pDirectionFit(this->GetCachedDirectionFit(pPrimaryCluster)); 
        EmergingCluster emergingPrimary(vertexPosition, nullptr, pDirectionFit);
        pandora::ClusterVector daughterClusters(this->GetOrderedDaughters(vertexPosition, pPrimaryCluster, inputClusterVector));

        EmergingClusterVector latestDaughter;
        EmergingCluster* pParentEmergingCluster(&emergingPrimary);
        pandora::CartesianVector parentOrigin(vertexPosition);

        for (const auto pDaughterCluster : daughterClusters)
        {
            TwoDSlidingFitResult currentSlidingFit(this->GetCachedSlidingFit(pDaughterCluster));
            pandora::CartesianVector currentOrigin((currentSlidingFit.GetGlobalMinLayerPosition() - parentOrigin).GetMagnitude() > (currentSlidingFit.GetGlobalMaxLayerPosition() - parentOrigin).GetMagnitude() ? currentSlidingFit.GetGlobalMinLayerPosition() : currentSlidingFit.GetGlobalMaxLayerPosition()); //whichever sliding fit endpoint is furthest away from the previous origin

            TrackDirectionTool::DirectionFitObject* pDaughterDirectionFit(this->GetCachedDirectionFit(pDaughterCluster)); 
            EmergingCluster currentEmergingCluster(currentOrigin, pParentEmergingCluster, pDaughterDirectionFit);

            latestDaughter.clear();
            latestDaughter.push_back(&currentEmergingCluster);

            parentOrigin = currentOrigin;
            pParentEmergingCluster = &currentEmergingCluster;
        }

        finalDaughters.push_back(latestDaughter.front());
    }

    return finalDaughters;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::ClusterVector DirectionFlowProbabilityAlgorithm::GetPrimaryClusters(const pandora::CartesianVector &positionVector, const pandora::ClusterVector &inputClusterVector) const
{
    ClusterToSpacepointsMap clusterToSpacepointsMap(this->FillClusterToSpacepointsMap(inputClusterVector));
    ClusterVector primaryClusterVector;

    for (const auto pCluster : inputClusterVector)
    {
        if (this->ClusterPointsToPosition(pCluster, positionVector, clusterToSpacepointsMap))
            primaryClusterVector.push_back(pCluster);
    }

    std::cout << "There are " << primaryClusterVector.size() << " primary clusters in this event." << std::endl;
    return primaryClusterVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::ClusterVector DirectionFlowProbabilityAlgorithm::GetOrderedDaughters(const pandora::CartesianVector &positionVector, const pandora::Cluster* const pPrimaryCluster, const pandora::ClusterVector &inputClusterVector) const
{
    ClusterToSpacepointsMap clusterToSpacepointsMap(this->FillClusterToSpacepointsMap(inputClusterVector));
    const pandora::Cluster* pCurrentCluster(pPrimaryCluster);
    TwoDSlidingFitResult currentSlidingFit(this->GetCachedSlidingFit(pCurrentCluster));
    pandora::CartesianVector currentOrigin((currentSlidingFit.GetGlobalMinLayerPosition() - positionVector).GetMagnitude() > (currentSlidingFit.GetGlobalMaxLayerPosition() - positionVector).GetMagnitude() ? currentSlidingFit.GetGlobalMinLayerPosition() : currentSlidingFit.GetGlobalMaxLayerPosition()); //whichever sliding fit endpoint is furthest away from the previous origin

    pandora::ClusterVector orderedClusterVector;
    bool daughterFound(true);

    while (daughterFound)
    {
        bool currentClusterHasDaughter(false);

        for (const auto pCluster : inputClusterVector)
        {
            if (this->ClusterPointsToPosition(pCluster, currentOrigin, clusterToSpacepointsMap))
            {
                pCurrentCluster = pCluster;
                currentSlidingFit = this->GetCachedSlidingFit(pCluster);
                currentOrigin = ((currentSlidingFit.GetGlobalMinLayerPosition() - currentOrigin).GetMagnitude() > (currentSlidingFit.GetGlobalMaxLayerPosition() - currentOrigin).GetMagnitude() ? currentSlidingFit.GetGlobalMinLayerPosition() : currentSlidingFit.GetGlobalMaxLayerPosition()); //whichever sliding fit endpoint is furthest away from the previous origin

                orderedClusterVector.push_back(pCluster);       

                currentClusterHasDaughter = true;
                break;
            }
        }

        daughterFound = currentClusterHasDaughter;
    }

    return orderedClusterVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionFlowProbabilityAlgorithm::ClusterPointsToPosition(const Cluster *const pCluster, const pandora::CartesianVector &positionVector, ClusterToSpacepointsMap &clusterToSpacepointsMap) const
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

DirectionFlowProbabilityAlgorithm::ClusterToSpacepointsMap DirectionFlowProbabilityAlgorithm::FillClusterToSpacepointsMap(const ClusterVector &clusterVector) const
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

void DirectionFlowProbabilityAlgorithm::GetSpacepoints(const Cluster *const pCluster, CartesianPointVector &spacepoints) const
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

void DirectionFlowProbabilityAlgorithm::AddToSlidingFitCache(const Cluster *const pCluster)
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
        TrackDirectionTool::DirectionFitObject directionFit(m_pTrackDirectionTool->GetClusterDirection(pCluster));
        m_directionFitMap.insert(std::make_pair(pCluster, &directionFit));
    } 
    catch (...)
    {
        std::cout << "Direction fit failure" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const TwoDSlidingFitResult &DirectionFlowProbabilityAlgorithm::GetCachedSlidingFit(const Cluster *const pCluster) const
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

TrackDirectionTool::DirectionFitObject* DirectionFlowProbabilityAlgorithm::GetCachedDirectionFit(const Cluster *const pCluster) 
{
    const auto iter = m_directionFitMap.find(pCluster);

    if (m_directionFitMap.end() == iter)
    {    
        std::cout << "Sliding fit retrieval failure" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }    

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

DirectionFlowProbabilityAlgorithm::EmergingCluster::EmergingCluster(const pandora::CartesianVector originPosition, const EmergingCluster* pParentEmergingCluster, TrackDirectionTool::DirectionFitObject* pDirectionFit) :
m_origin(&originPosition),
m_parent(pParentEmergingCluster),
m_directionfit(pDirectionFit) 
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DirectionFlowProbabilityAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImpactRadius", m_impactRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NExtrapolationSteps", m_extrapolationNSteps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumClusterLength", m_minimumClusterLength));

    AlgorithmTool *pAlgorithmTool(nullptr);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "TrackDirection", pAlgorithmTool));

    if (!(this->m_pTrackDirectionTool = dynamic_cast<TrackDirectionTool *>(pAlgorithmTool)))
        throw STATUS_CODE_FAILURE;

    return VertexSelectionBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
