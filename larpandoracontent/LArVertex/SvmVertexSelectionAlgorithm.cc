/**
 *  @file   larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.cc
 *
 *  @brief  Implementation of the Svm vertex selection algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArVertex/EnergyKickFeatureTool.h"
#include "larpandoracontent/LArVertex/LocalAsymmetryFeatureTool.h"
#include "larpandoracontent/LArVertex/GlobalAsymmetryFeatureTool.h"
#include "larpandoracontent/LArVertex/ShowerAsymmetryFeatureTool.h"
#include "larpandoracontent/LArVertex/RPhiFeatureTool.h"

#include "larpandoracontent/LArVertex/SvmVertexSelectionAlgorithm.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include "larpandoracontent/LArHelpers/LArSpaceChargeHelper.h"

#include <random>

using namespace pandora;

namespace lar_content
{
SvmVertexSelectionAlgorithm::SvmVertexSelectionAlgorithm() :
    VertexSelectionBaseAlgorithm(),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH"),
    m_trainingSetMode(false),
    m_allowClassifyDuringTraining(false),
    m_mcVertexXCorrection(0.f),
    m_minClusterCaloHits(12),
    m_slidingFitWindow(100),
    m_minShowerSpineLength(15.f),
    m_beamDeweightingConstant(0.4),
    m_localAsymmetryConstant(3.f),
    m_globalAsymmetryConstant(1.f),
    m_showerAsymmetryConstant(1.f),
    m_energyKickConstant(0.06),
    m_showerClusteringDistance(3.f),
    m_minShowerClusterHits(1),
    m_useShowerClusteringApproximation(false),
    m_regionRadius(5.f),
    m_rPhiFineTuningRadius(2.f),
    m_maxTrueVertexRadius(1.f),
    m_useRPhiFeatureForRegion(false),
    m_dropFailedRPhiFastScoreCandidates(true),
    m_enableIntersectionFeatures(true),
    m_enableDirectionFeatures(false), //whether to enable direction features
    m_directionScoreReweighting(false),
    m_onlyCorrectDirection(false),
    m_extrapolationNSteps(200),
    m_extrapolationStepSize(0.1f),
    m_visualiseEvent(false),
    m_fileIdentifier(1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::GetVertexScoreList(const VertexVector &vertexVector, const BeamConstants &beamConstants,
    HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const
{
    ClusterList clustersU, clustersV, clustersW;
    this->GetClusterLists(m_inputClusterListNames, clustersU, clustersV, clustersW);

    SlidingFitDataList slidingFitDataListU, slidingFitDataListV, slidingFitDataListW;
    this->CalculateClusterSlidingFits(clustersU, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListU);
    this->CalculateClusterSlidingFits(clustersV, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListV);
    this->CalculateClusterSlidingFits(clustersW, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListW);

    ShowerClusterList showerClusterListU, showerClusterListV, showerClusterListW;
    this->CalculateShowerClusterList(clustersU, showerClusterListU);
    this->CalculateShowerClusterList(clustersV, showerClusterListV);
    this->CalculateShowerClusterList(clustersW, showerClusterListW);

    // Create maps from hit types to objects for passing to feature tools.
    const ClusterListMap clusterListMap{{TPC_VIEW_U, clustersU},
                                        {TPC_VIEW_V, clustersV},
                                        {TPC_VIEW_W, clustersV}};

    const SlidingFitDataListMap slidingFitDataListMap{{TPC_VIEW_U, slidingFitDataListU},
                                                      {TPC_VIEW_V, slidingFitDataListV},
                                                      {TPC_VIEW_W, slidingFitDataListW}};

    const ShowerClusterListMap showerClusterListMap{{TPC_VIEW_U, showerClusterListU},
                                                    {TPC_VIEW_V, showerClusterListV},
                                                    {TPC_VIEW_W, showerClusterListW}};

    const KDTreeMap kdTreeMap{{TPC_VIEW_U, kdTreeU},
                              {TPC_VIEW_V, kdTreeV},
                              {TPC_VIEW_W, kdTreeW}};

    // Calculate the event feature list and the vertex feature map.
    EventFeatureInfo eventFeatureInfo(this->CalculateEventFeatures(clustersU, clustersV, clustersW, vertexVector));

    LArMvaHelper::MvaFeatureVector eventFeatureList;
    this->AddEventFeaturesToVector(eventFeatureInfo, eventFeatureList);

    float longestClusterDirectionProbability(0.5), minChiSquaredPerHit(100.0), longestClusterLength(1.f); 
    pandora::CartesianVector beginpointPosition(1.f, 1.f, 1.f);

    try
    {
        const pandora::Cluster* const pLongestCluster(m_pDirectionFlowProbabilityTool->GetLongestCluster(clustersW));
        TrackDirectionTool::DirectionFitObject directionFit = m_pTrackDirectionTool->GetClusterDirection(pLongestCluster); 

        longestClusterDirectionProbability = (directionFit.GetProbability());
        minChiSquaredPerHit = directionFit.GetMinChiSquaredPerHit();
        longestClusterLength = LArClusterHelper::GetLength(pLongestCluster); 
        beginpointPosition = (directionFit.GetBeginpoint());
    }
    catch (...)
    {
        std::cout << "Could not determine longest cluster direction probability." << std::endl;
    }

    if (m_visualiseEvent)
    {
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clustersU, "All Clusters U", BLACK));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clustersV, "All Clusters V", BLACK));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clustersW, "All Clusters W", BLACK));
    }

    VertexFeatureInfoMap vertexFeatureInfoMap;
    for (const Vertex *const pVertex : vertexVector)
    {
        if (m_visualiseEvent)
        {
            const pandora::CartesianVector vertexProjection(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjection, "Vertex Projection", BLUE, 1));
        }

        this->PopulateVertexFeatureInfoMap(beamConstants, longestClusterDirectionProbability, minChiSquaredPerHit, longestClusterLength, beginpointPosition, clusterListMap, slidingFitDataListMap, showerClusterListMap, kdTreeMap, pVertex,
            vertexFeatureInfoMap);
    }

    if (m_visualiseEvent)
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    // Use a simple score to get the list of vertices representing good regions.

    VertexScoreList initialScoreList;
    for (const Vertex *const pVertex : vertexVector)
        PopulateInitialScoreList(vertexFeatureInfoMap, pVertex, initialScoreList);

    VertexVector bestRegionVertices;
    this->GetBestRegionVertices(initialScoreList, bestRegionVertices);

    if (m_visualiseEvent)
    {
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clustersU, "All Clusters U", BLACK));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clustersV, "All Clusters V", BLACK));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clustersW, "All Clusters W", BLACK));

        for (const Vertex *const pVertex : bestRegionVertices)
        {
            const pandora::CartesianVector vertexProjection(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjection, "Region Vertex Projection", GREEN, 1));
        }

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    if (m_trainingSetMode)
        this->ProduceTrainingSets(vertexVector, bestRegionVertices, vertexFeatureInfoMap, eventFeatureList, kdTreeMap);

    if ((!m_trainingSetMode || m_allowClassifyDuringTraining) && !bestRegionVertices.empty())
    {
        // Use svm to choose the region.
        const Vertex* pBestRegionVertex;

        if (m_directionScoreReweighting)
        {
            try
            {
                pBestRegionVertex = this->ScoreVertices(bestRegionVertices, clustersW, vertexFeatureInfoMap, eventFeatureList, m_svMachineRegion,
                m_useRPhiFeatureForRegion, m_enableDirectionFeatures);
            }
            catch (...)
            {
                std::cout << "Something went wrong, reverting to normal vertex region selection." << std::endl;

                pBestRegionVertex = this->CompareVertices(bestRegionVertices, vertexFeatureInfoMap, eventFeatureList, m_svMachineRegion,
                m_useRPhiFeatureForRegion, m_enableDirectionFeatures);
            }
        }
        else
        {
            pBestRegionVertex = this->CompareVertices(bestRegionVertices, vertexFeatureInfoMap, eventFeatureList, m_svMachineRegion,
            m_useRPhiFeatureForRegion, m_enableDirectionFeatures);
        }

        // Get all the vertices in the best region.
        VertexVector regionalVertices{pBestRegionVertex};
        for (const Vertex *const pVertex : vertexVector)
        {
            if (pVertex == pBestRegionVertex)
                continue;

            if ((pBestRegionVertex->GetPosition() - pVertex->GetPosition()).GetMagnitude() < m_regionRadius)
                regionalVertices.push_back(pVertex);
        }

        this->CalculateRPhiScores(regionalVertices, vertexFeatureInfoMap, kdTreeMap);

        //if (m_enableIntersectionFeatures)
        //    this->CalculateIntersectionScores(regionalVertices, vertexFeatureInfoMap, clustersW);

        if (!regionalVertices.empty())
        {
            // Use svm to choose the vertex and then fine-tune using the RPhi score.
            const Vertex *const pBestVertex(this->CompareVertices(regionalVertices, vertexFeatureInfoMap, eventFeatureList, m_svMachineVertex, true, false));
            this->PopulateFinalVertexScoreList(vertexFeatureInfoMap, pBestVertex, vertexVector, vertexScoreList);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::CalculateShowerClusterList(const ClusterList &inputClusterList, ShowerClusterList &showerClusterList) const
{
    ClusterEndPointsMap clusterEndPointsMap;
    ClusterList showerLikeClusters;
    this->GetShowerLikeClusterEndPoints(inputClusterList, showerLikeClusters, clusterEndPointsMap);

    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    ClusterList availableShowerLikeClusters(showerLikeClusters.begin(), showerLikeClusters.end());

    HitKDTree2D kdTree;
    HitToClusterMap hitToClusterMap;

    if (!m_useShowerClusteringApproximation)
        this->PopulateKdTree(availableShowerLikeClusters, kdTree, hitToClusterMap);

    while (!availableShowerLikeClusters.empty())
    {
        ClusterList showerCluster;
        showerCluster.push_back(availableShowerLikeClusters.back());
        availableShowerLikeClusters.pop_back();

        bool addedCluster(true);
        while (addedCluster && !availableShowerLikeClusters.empty())
        {
            addedCluster = false;
            for (const Cluster *const pCluster : showerCluster)
            {
                if (!m_useShowerClusteringApproximation)
                {
                    addedCluster = this->AddClusterToShower(kdTree, hitToClusterMap, availableShowerLikeClusters, pCluster, showerCluster);
                }
                else
                {
                    addedCluster = this->AddClusterToShower(clusterEndPointsMap, availableShowerLikeClusters, pCluster, showerCluster);
                }

                if (addedCluster)
                    break;
            }
        }
        unsigned int totHits(0);
        for (const Cluster *const pCluster : showerCluster)
            totHits += pCluster->GetNCaloHits();

        if (totHits < m_minClusterCaloHits)
            continue;

        showerClusterList.emplace_back(showerCluster, slidingFitPitch, m_slidingFitWindow);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::GetShowerLikeClusterEndPoints(const ClusterList &clusterList, ClusterList &showerLikeClusters,
    ClusterEndPointsMap &clusterEndPointsMap) const
{
    for (const Cluster *const pCluster : clusterList)
    {
        if (pCluster->GetNCaloHits() < m_minShowerClusterHits)
            continue;

        if (this->IsClusterShowerLike(pCluster))
            showerLikeClusters.push_back(pCluster);

        CaloHitList clusterCaloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);

        CaloHitVector clusterCaloHitVector(clusterCaloHitList.begin(), clusterCaloHitList.end());
        std::sort(clusterCaloHitVector.begin(), clusterCaloHitVector.end(), LArClusterHelper::SortHitsByPosition);

        if (clusterCaloHitVector.empty())
            continue;

        ClusterEndPoints clusterEndPoints(clusterCaloHitVector.front()->GetPositionVector(), clusterCaloHitVector.back()->GetPositionVector());
        clusterEndPointsMap.emplace(pCluster, clusterEndPoints);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::PopulateKdTree(const ClusterList &clusterList, HitKDTree2D &kdTree, HitToClusterMap &hitToClusterMap) const
{
    CaloHitList allCaloHits;

    for (const Cluster *const pCluster : clusterList)
    {
        CaloHitList daughterHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(daughterHits);
        allCaloHits.insert(allCaloHits.end(), daughterHits.begin(), daughterHits.end());

        for (const CaloHit *const pCaloHit : daughterHits)
            (void) hitToClusterMap.insert(HitToClusterMap::value_type(pCaloHit, pCluster));
    }

    HitKDNode2DList hitKDNode2DList;
    KDTreeBox hitsBoundingRegion2D(fill_and_bound_2d_kd_tree(allCaloHits, hitKDNode2DList));
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SvmVertexSelectionAlgorithm::AddClusterToShower(const ClusterEndPointsMap &clusterEndPointsMap, ClusterList &availableShowerLikeClusters,
    const Cluster *const pCluster, ClusterList &showerCluster) const
{
    const auto existingEndPointsIter(clusterEndPointsMap.find(pCluster));
    if (existingEndPointsIter == clusterEndPointsMap.end())
        return false;

    const ClusterEndPoints &existingClusterEndPoints(existingEndPointsIter->second);

    for (auto iter = availableShowerLikeClusters.begin(); iter != availableShowerLikeClusters.end(); ++iter)
    {
        const Cluster *const pAvailableShowerLikeCluster(*iter);
        const auto &newEndPointsIter(clusterEndPointsMap.find(pAvailableShowerLikeCluster));

        if (newEndPointsIter == clusterEndPointsMap.end())
            continue;

        const ClusterEndPoints &newClusterEndPoints(newEndPointsIter->second);
        const float startStartDistance((newClusterEndPoints.first - existingClusterEndPoints.first).GetMagnitude());
        const float startEndDistance((newClusterEndPoints.first - existingClusterEndPoints.second).GetMagnitude());
        const float endStartDistance((newClusterEndPoints.second - existingClusterEndPoints.first).GetMagnitude());
        const float endEndDistance((newClusterEndPoints.second - existingClusterEndPoints.second).GetMagnitude());

        const float smallestDistance(std::min(std::min(startStartDistance, startEndDistance), std::min(endStartDistance, endEndDistance)));

        if (smallestDistance < m_showerClusteringDistance)
        {
            showerCluster.push_back(pAvailableShowerLikeCluster);
            availableShowerLikeClusters.erase(iter); // Now must return, after invalidating current iterator
            return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SvmVertexSelectionAlgorithm::AddClusterToShower(HitKDTree2D &kdTree, const HitToClusterMap &hitToClusterMap,
    ClusterList &availableShowerLikeClusters, const Cluster *const pCluster, ClusterList &showerCluster) const
{
    ClusterSet nearbyClusters;
    CaloHitList daughterHits;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(daughterHits);

    for (const CaloHit *const pCaloHit : daughterHits)
    {
        KDTreeBox searchRegionHits = build_2d_kd_search_region(pCaloHit, m_showerClusteringDistance, m_showerClusteringDistance);

        HitKDNode2DList found;
        kdTree.search(searchRegionHits, found);

        for (const auto &hit : found)
            (void) nearbyClusters.insert(hitToClusterMap.at(hit.data));
    }

    for (auto iter = availableShowerLikeClusters.begin(); iter != availableShowerLikeClusters.end(); ++iter)
    {
        const Cluster *const pAvailableShowerLikeCluster(*iter);

        if ((pAvailableShowerLikeCluster != pCluster) && nearbyClusters.count(pAvailableShowerLikeCluster))
        {
            showerCluster.push_back(pAvailableShowerLikeCluster);
            availableShowerLikeClusters.erase(iter); // Now must return, after invalidating current iterator
            return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

SvmVertexSelectionAlgorithm::EventFeatureInfo SvmVertexSelectionAlgorithm::CalculateEventFeatures(const ClusterList &clusterListU,
    const ClusterList &clusterListV, const ClusterList &clusterListW, const VertexVector &vertexVector) const
{
    float eventEnergy(0.f);
    unsigned int nShoweryHits(0), nHits(0);

    this->IncrementShoweryParameters(clusterListU, nShoweryHits, nHits, eventEnergy);
    this->IncrementShoweryParameters(clusterListV, nShoweryHits, nHits, eventEnergy);
    this->IncrementShoweryParameters(clusterListW, nShoweryHits, nHits, eventEnergy);

    const unsigned int nClusters(clusterListU.size() + clusterListV.size() + clusterListW.size());
    const float eventShoweryness((nHits > 0) ? static_cast<float>(nShoweryHits) / static_cast<float>(nHits) : 0.f);

    ClusterList allClusters(clusterListU);
    allClusters.insert(allClusters.end(), clusterListV.begin(), clusterListV.end());
    allClusters.insert(allClusters.end(), clusterListW.begin(), clusterListW.end());

    float eventVolume(0.f), longitudinality(0.f);
    this->GetEventShapeFeatures(allClusters, eventVolume, longitudinality);

    return EventFeatureInfo(eventShoweryness, eventEnergy, eventVolume, longitudinality, nHits, nClusters, vertexVector.size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::IncrementShoweryParameters(const ClusterList &clusterList, unsigned int &nShoweryHits, unsigned int &nHits,
    float &eventEnergy) const
{
    for (const Cluster *const pCluster : clusterList)
    {
        if (this->IsClusterShowerLike(pCluster))
            nShoweryHits += pCluster->GetNCaloHits();

        eventEnergy += pCluster->GetElectromagneticEnergy();
        nHits += pCluster->GetNCaloHits();

        OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        CaloHitList caloHitList;
        orderedCaloHitList.FillCaloHitList(caloHitList);    
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool SvmVertexSelectionAlgorithm::IsClusterShowerLike(const Cluster *const pCluster) const
{
    return (pCluster->GetParticleId() == E_MINUS && LArClusterHelper::GetLength(pCluster) < m_minShowerSpineLength);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::GetEventShapeFeatures(const ClusterList &clusterList, float &eventVolume, float &longitudinality) const
{
    InputFloat xMin, yMin, zMin, xMax, yMax, zMax;

    for (const Cluster *const pCluster : clusterList)
    {
        CartesianVector minPosition(0.f, 0.f, 0.f), maxPosition(0.f, 0.f, 0.f);
        LArClusterHelper::GetClusterBoundingBox(pCluster, minPosition, maxPosition);

        this->UpdateSpanCoordinate(minPosition.GetX(), maxPosition.GetX(), xMin, xMax);
        this->UpdateSpanCoordinate(minPosition.GetY(), maxPosition.GetY(), yMin, yMax);
        this->UpdateSpanCoordinate(minPosition.GetZ(), maxPosition.GetZ(), zMin, zMax);
    }

    const float xSpan(this->GetCoordinateSpan(xMax, xMin));
    const float ySpan(this->GetCoordinateSpan(yMax, zMin));
    const float zSpan(this->GetCoordinateSpan(yMax, zMin));

    // Calculate the volume and longitudinality of the event (ySpan often 0 - to be investigated).
    if ((xSpan > std::numeric_limits<float>::epsilon()) && (ySpan > std::numeric_limits<float>::epsilon()))
    {
        eventVolume     = xSpan * ySpan * zSpan;
        longitudinality = zSpan / (xSpan + ySpan);
    }

    else if (ySpan > std::numeric_limits<float>::epsilon())
    {
        eventVolume     = ySpan * ySpan * zSpan;
        longitudinality = zSpan / (ySpan + ySpan);
    }

    else if (xSpan > std::numeric_limits<float>::epsilon())
    {
        eventVolume     = xSpan * xSpan * zSpan;
        longitudinality = zSpan / (xSpan + xSpan);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void SvmVertexSelectionAlgorithm::UpdateSpanCoordinate(const float minPositionCoord, const float maxPositionCoord, InputFloat &minCoord,
    InputFloat &maxCoord) const
{
    if (!minCoord.IsInitialized() || minPositionCoord < minCoord.Get())
        minCoord = minPositionCoord;

    if (!maxCoord.IsInitialized() || maxPositionCoord > maxCoord.Get())
        maxCoord = maxPositionCoord;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SvmVertexSelectionAlgorithm::GetCoordinateSpan(const InputFloat &minCoord, const InputFloat &maxCoord) const
{
   if (minCoord.IsInitialized() && maxCoord.IsInitialized())
        return std::fabs(maxCoord.Get() - minCoord.Get());

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::AddEventFeaturesToVector(const EventFeatureInfo &eventFeatureInfo,
    LArMvaHelper::MvaFeatureVector &featureVector) const
{
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_eventShoweryness));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_eventEnergy));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_eventVolume));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_longitudinality));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_nHits));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_nClusters));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_nCandidates));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::PopulateVertexFeatureInfoMap(const BeamConstants &beamConstants, float &longestClusterDirectionProbability, float &minChiSquaredPerHit, float &longestClusterLength, pandora::CartesianVector &beginpointPosition, const ClusterListMap &clusterListMap,
    const SlidingFitDataListMap &slidingFitDataListMap, const ShowerClusterListMap &showerClusterListMap, const KDTreeMap &kdTreeMap,
    const Vertex *const pVertex, VertexFeatureInfoMap &vertexFeatureInfoMap) const
{
    float bestFastScore(-std::numeric_limits<float>::max()); // not actually used - artefact of toolizing RPhi score and still using performance trick

    const double beamDeweighting(this->GetBeamDeweightingScore(beamConstants, pVertex));

    const double energyKick(LArMvaHelper::CalculateFeaturesOfType<EnergyKickFeatureTool>(m_featureToolVector, this, pVertex,
        slidingFitDataListMap, clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0).Get());

    const double localAsymmetry(LArMvaHelper::CalculateFeaturesOfType<LocalAsymmetryFeatureTool>(m_featureToolVector, this, pVertex,
        slidingFitDataListMap, clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0).Get());

    const double globalAsymmetry(LArMvaHelper::CalculateFeaturesOfType<GlobalAsymmetryFeatureTool>(m_featureToolVector, this, pVertex,
        slidingFitDataListMap, clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0).Get());

    const double showerAsymmetry(LArMvaHelper::CalculateFeaturesOfType<ShowerAsymmetryFeatureTool>(m_featureToolVector, this, pVertex,
        slidingFitDataListMap, clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0).Get());

    //const double rPhiFeature(LArMvaHelper::CalculateFeaturesOfType<RPhiFeatureTool>(m_featureToolVector, this, pVertex,
    //    slidingFitDataListMap, clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0).Get());

    double directionFeature(0), directionBeginpointDistance(0), directionFitQuality(0);

    if (m_enableDirectionFeatures)
    {
        const pandora::CartesianVector vertexProjection(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W));

        directionFeature = longestClusterDirectionProbability;
        directionBeginpointDistance = (vertexProjection - beginpointPosition).GetMagnitude()/longestClusterLength;
        directionFitQuality = minChiSquaredPerHit;
    }

    VertexFeatureInfo vertexFeatureInfo(beamDeweighting, 0.f, directionFeature, directionBeginpointDistance, directionFitQuality, energyKick, localAsymmetry, globalAsymmetry, showerAsymmetry);
    vertexFeatureInfoMap.emplace(pVertex, vertexFeatureInfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::PopulateInitialScoreList(VertexFeatureInfoMap &vertexFeatureInfoMap, const Vertex *const pVertex,
                                                           VertexScoreList &initialScoreList) const
{
    VertexFeatureInfo vertexFeatureInfo = vertexFeatureInfoMap.at(pVertex);

    const float beamDeweightingScore(vertexFeatureInfo.m_beamDeweighting / m_beamDeweightingConstant);
    const float energyKickScore(-vertexFeatureInfo.m_energyKick / m_energyKickConstant);
    const float localAsymmetryScore(vertexFeatureInfo.m_localAsymmetry / m_localAsymmetryConstant);
    const float globalAsymmetryScore(vertexFeatureInfo.m_globalAsymmetry / m_globalAsymmetryConstant);
    const float showerAsymmetryScore(vertexFeatureInfo.m_showerAsymmetry / m_showerAsymmetryConstant);

    initialScoreList.emplace_back(pVertex, beamDeweightingScore + energyKickScore + localAsymmetryScore + globalAsymmetryScore + showerAsymmetryScore);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::GetBestRegionVertices(VertexScoreList &initialScoreList, VertexVector &bestRegionVertices) const
{
    std::sort(initialScoreList.begin(), initialScoreList.end());

    for (const VertexScore &vertexScore : initialScoreList)
    {
        const Vertex *const pVertex(vertexScore.GetVertex());
        bool farEnoughAway(true);

        for (const Vertex *const pRegionVertex: bestRegionVertices)
        {
            if (pRegionVertex == pVertex)
            {
                farEnoughAway = false;
                break;
            }

            const float distance = (pRegionVertex->GetPosition() - pVertex->GetPosition()).GetMagnitude();

            if (distance <= m_regionRadius)
            {
                farEnoughAway = false;
                break;
            }
        }

        if (farEnoughAway)
            bestRegionVertices.push_back(pVertex);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::ProduceTrainingSets(const VertexVector &vertexVector, const VertexVector &bestRegionVertices,
    VertexFeatureInfoMap &vertexFeatureInfoMap, const LArMvaHelper::MvaFeatureVector &eventFeatureList, const KDTreeMap &kdTreeMap) const
{
    if (vertexVector.empty())
        return;

    // Create a distribution for random coin flips.
    std::random_device device;
    std::mt19937 generator(device());
    std::bernoulli_distribution coinFlip(0.5);

    const std::string interactionType(this->GetInteractionType());

    // Produce training examples for the vertices representing regions.
    const Vertex *const pBestRegionVertex(this->ProduceTrainingExamples(bestRegionVertices, vertexFeatureInfoMap, coinFlip, generator,
        interactionType, m_trainingOutputFileRegion, eventFeatureList, m_regionRadius, m_useRPhiFeatureForRegion, m_enableDirectionFeatures));

    // Get all the vertices in the best region.
    VertexVector regionalVertices{pBestRegionVertex};
    for (const Vertex *const pVertex : vertexVector)
    {
        if (pVertex == pBestRegionVertex)
            continue;

        if ((pBestRegionVertex->GetPosition() - pVertex->GetPosition()).GetMagnitude() < m_regionRadius)
            regionalVertices.push_back(pVertex);
    }

    this->CalculateRPhiScores(regionalVertices, vertexFeatureInfoMap, kdTreeMap);

    //if (m_enableIntersectionFeatures)
    //    this->CalculateIntersectionScores(regionalVertices, vertexFeatureInfoMap, clustersW);

    // Produce training examples for the final vertices within the best region.
    if (!regionalVertices.empty())
    {
        this->ProduceTrainingExamples(regionalVertices, vertexFeatureInfoMap, coinFlip, generator, interactionType, m_trainingOutputFileVertex,
            eventFeatureList, m_maxTrueVertexRadius, true, false);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::CalculateRPhiScores(VertexVector &vertexVector, VertexFeatureInfoMap &vertexFeatureInfoMap,
    const KDTreeMap &kdTreeMap) const
{
    float bestFastScore(-std::numeric_limits<float>::max());

    for (auto iter = vertexVector.begin(); iter != vertexVector.end(); /* no increment */)
    {
        VertexFeatureInfo &vertexFeatureInfo = vertexFeatureInfoMap.at(*iter);
        vertexFeatureInfo.m_rPhiFeature = static_cast<float>(LArMvaHelper::CalculateFeaturesOfType<RPhiFeatureTool>(m_featureToolVector, this, *iter,
            SlidingFitDataListMap(), ClusterListMap(), kdTreeMap, ShowerClusterListMap(), vertexFeatureInfo.m_beamDeweighting,
            bestFastScore).at(0).Get());

        if (m_dropFailedRPhiFastScoreCandidates && (vertexFeatureInfo.m_rPhiFeature <= std::numeric_limits<float>::epsilon()))
            iter = vertexVector.erase(iter);

        else
            ++iter;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
/*
void SvmVertexSelectionAlgorithm::CalculateIntersectionScores(VertexVector &vertexVector, VertexFeatureInfoMap &vertexFeatureInfoMap,
    ClusterList &clusterList) const
{
    //PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, "All Clusters", BLACK));
    //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    int nViableClusters(CountViableClusters(clusterList));
    //std::cout << "nViableClusters: " << nViableClusters << std::endl;

    CartesianVector intersectionPoint(GetIntersectionPoint(clusterList));

    //PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, "All Clusters", BLACK));

    for (const auto pVertex : vertexVector)
    {
        float intersectionFeature(0.f);

        const pandora::CartesianVector vertexProjection(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W));
        //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjection, "Vertex Projection", RED, 1));

        if  (nViableClusters > 1 && intersectionPoint.GetMagnitude() > 0)
            intersectionFeature = (vertexProjection - intersectionPoint).GetMagnitude();
        else
            intersectionFeature = GetMinimalExtremalZDistance(clusterList, vertexProjection);

        std::cout << "intersectionFeature: " << intersectionFeature << std::endl;

        VertexFeatureInfo &vertexFeatureInfo = vertexFeatureInfoMap.at(pVertex);
        vertexFeatureInfo.m_intersectionFeature = intersectionFeature; 
    }

    //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

int SvmVertexSelectionAlgorithm::CountViableClusters(ClusterList &clusterList) const
{
    int nViableClusters(0);

    for (const auto pCluster : clusterList)
    {
        //std::cout << "L: " << LArClusterHelper::GetLength(pCluster) << std::endl;
        if (LArClusterHelper::GetLength(pCluster) > 5.0)
            ++nViableClusters;
    }

    return nViableClusters;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float SvmVertexSelectionAlgorithm::GetMinimalExtremalZDistance(ClusterList &clusterList, const CartesianVector &vertexProjection) const
{
    pandora::ClusterList longestClusterList;
    float longestClusterLength(0.f);

    for (const auto pCluster : clusterList)
    {
        if (LArClusterHelper::GetLength(pCluster) > longestClusterLength)
        {
            longestClusterLength = LArClusterHelper::GetLength(pCluster);
            longestClusterList.clear();
            longestClusterList.insert(longestClusterList.begin(), pCluster);
        }
    }

    pandora::CartesianVector innerCoordinate(0.f, 0.f, 0.f), outerCoordinate(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(longestClusterList, innerCoordinate, outerCoordinate);

    float innerDistance((vertexProjection - innerCoordinate).GetMagnitude()), outerDistance((vertexProjection - outerCoordinate).GetMagnitude());

    return (innerDistance < outerDistance ? innerDistance : outerDistance); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector SvmVertexSelectionAlgorithm::GetIntersectionPoint(ClusterList &clusterList) const
{
    int nViableClusters(CountViableClusters(clusterList));
    CartesianVector origin(0.f, 0.f, 0.f);
    
    if (nViableClusters < 2)
        return origin;

    ClusterList twoLargestClusters;
    GetTwoLargestClusters(clusterList, twoLargestClusters);

    //PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &twoLargestClusters, "Two Largest Clusters", BLACK));
    //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    std::vector<CartesianPointVector> spacePointsVector;

    for (const auto pCluster : twoLargestClusters)
    {
        CartesianPointVector spacepoints;
        GetSpacepoints(pCluster, spacepoints);
        spacePointsVector.push_back(spacepoints);
    }

    float closestApproach(1e6);
    CartesianVector intersectionPoint(0.f, 0.f, 0.f);

    for (const auto &position1 : spacePointsVector.front())
    {
        //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position1, "Extrapolation", BLACK, 1));

        for (const auto &position2 : spacePointsVector.back())
        {
            //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position2, "Extrapolation", BLACK, 1));

            if ((position1 - position2).GetMagnitude() < closestApproach)     
            {
                closestApproach = (position1 - position2).GetMagnitude();
                intersectionPoint.SetValues((position1.GetX() + position2.GetX())/2, (position1.GetY() + position2.GetY())/2, (position1.GetZ() + position2.GetZ())/2);
            }
        }
    }

    //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    std::cout << "closestApproach: " << closestApproach << std::endl;

    if (closestApproach < 1.0)
        return intersectionPoint;
    else
        return origin;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::GetTwoLargestClusters(ClusterList &clusterList, ClusterList &twoLargestClusters) const
{
    ClusterVector clusterVector(clusterList.begin(), clusterList.end());

    std::sort(clusterVector.begin(), clusterVector.end(), ClusterLengthComparison);
    twoLargestClusters.insert(twoLargestClusters.begin(), clusterVector.at(0));
    twoLargestClusters.insert(twoLargestClusters.begin(), clusterVector.at(1));

    //std::cout << "L1: " << LArClusterHelper::GetLength(clusterVector.at(0)) << std::endl;
    //std::cout << "L2: " << LArClusterHelper::GetLength(clusterVector.at(1)) << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SvmVertexSelectionAlgorithm::ClusterLengthComparison(const Cluster *const pCluster1, const Cluster *const pCluster2)
{
    return LArClusterHelper::GetLength(pCluster1) > LArClusterHelper::GetLength(pCluster2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::GetSpacepoints(const Cluster *const pCluster, CartesianPointVector &spacepoints) const
{
    LArClusterHelper::GetCoordinateVector(pCluster, spacepoints);

    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult fitResult(pCluster, m_slidingFitWindow, slidingFitPitch);

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
*/

std::string SvmVertexSelectionAlgorithm::GetInteractionType() const
{
    // Extract input collections
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    // ATTN Assumes single neutrino is parent of all neutrino-induced mc particles
    LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, LArMCParticleHelper::PrimaryParameters(),
        LArMCParticleHelper::IsBeamNeutrinoFinalState, nuMCParticlesToGoodHitsMap);

    MCParticleList mcPrimaryList;
    for (const auto &mapEntry : nuMCParticlesToGoodHitsMap) mcPrimaryList.push_back(mapEntry.first);
    mcPrimaryList.sort(LArMCParticleHelper::SortByMomentum);

    const LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(mcPrimaryList));
    return LArInteractionTypeHelper::ToString(interactionType);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Vertex * SvmVertexSelectionAlgorithm::ProduceTrainingExamples(const VertexVector &vertexVector,
    const VertexFeatureInfoMap &vertexFeatureInfoMap, std::bernoulli_distribution &coinFlip, std::mt19937 &generator,
    const std::string &interactionType, const std::string &trainingOutputFile, const LArMvaHelper::MvaFeatureVector &eventFeatureList,
    const float maxRadius, const bool useRPhi, const bool useDirection) const
{
    const Vertex *pBestVertex(nullptr);
    float bestVertexDr(std::numeric_limits<float>::max());

    LArMvaHelper::MvaFeatureVector bestVertexFeatureList;
    this->GetBestVertex(vertexVector, pBestVertex, bestVertexDr);

    VertexFeatureInfo bestVertexFeatureInfo(vertexFeatureInfoMap.at(pBestVertex));
    this->AddVertexFeaturesToVector(bestVertexFeatureInfo, bestVertexFeatureList, useRPhi, useDirection);

    std::string fileString(std::to_string(m_fileIdentifier)); 

    /*
    std::string directionString("UNKNOWN");

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));
    
    for (const auto pMCParticle : *pMCParticleList)
    {
        if (pMCParticle->GetParticleId() == 13 && (pMCParticle->GetEndpoint() - pMCParticle->GetVertex()).GetMagnitude() > 30.0)
        {
            if (pMCParticle->GetVertex().GetZ() + 20.0 < pMCParticle->GetEndpoint().GetZ())
                directionString = "FORWARDS";
            else if (pMCParticle->GetVertex().GetZ() > pMCParticle->GetEndpoint().GetZ() + 20.0)
                directionString = "BACKWARDS";
            else
                directionString = "UNKNOWN";
        }
    }
    */

    for (const Vertex *const pVertex : vertexVector)
    {
        if (pVertex == pBestVertex)
            continue;

        LArMvaHelper::MvaFeatureVector featureList;
        VertexFeatureInfo vertexFeatureInfo(vertexFeatureInfoMap.at(pVertex));
        this->AddVertexFeaturesToVector(vertexFeatureInfo, featureList, useRPhi, useDirection);

        if (pBestVertex && (bestVertexDr < maxRadius))
        {
            if (coinFlip(generator))
            {
                LArMvaHelper::ProduceTrainingExample("/usera/jjd49/pandora_direction/CondorUtilities/txts/" + fileString + "_" + trainingOutputFile + "_" + interactionType + ".txt", true, eventFeatureList,
                    bestVertexFeatureList, featureList);
            }

            else
            {
                LArMvaHelper::ProduceTrainingExample("/usera/jjd49/pandora_direction/CondorUtilities/txts/" + fileString + "_" + trainingOutputFile + "_" + interactionType + ".txt", false, eventFeatureList, featureList,
                    bestVertexFeatureList);
            }
        }
    }

    return pBestVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::GetBestVertex(const VertexVector &vertexVector, const Vertex *&pBestVertex, float &bestVertexDr) const
{
    // Extract input collections
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    // Obtain vector: true neutrinos
    MCParticleVector mcNeutrinoVector;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, mcNeutrinoVector);

    for (const Vertex *const pVertex : vertexVector)
    {
        CartesianVector vertexPosition(pVertex->GetPosition());
        const CartesianVector correctedVertexPosition(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(vertexPosition));
 
        float mcVertexDr(std::numeric_limits<float>::max());
        for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
        {
            const CartesianVector mcNeutrinoPosition(pMCNeutrino->GetEndpoint().GetX() + m_mcVertexXCorrection, pMCNeutrino->GetEndpoint().GetY(),
                pMCNeutrino->GetEndpoint().GetZ());

            const float dr = (mcNeutrinoPosition - correctedVertexPosition).GetMagnitude();

            if (dr < mcVertexDr)
                mcVertexDr = dr;
        }

        if (std::abs(mcVertexDr - 0.5) < bestVertexDr) //Need to account for 0.5cm MC/reco offset
        {
            bestVertexDr = mcVertexDr;
            pBestVertex = pVertex;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::AddVertexFeaturesToVector(const VertexFeatureInfo &vertexFeatureInfo,
    LArMvaHelper::MvaFeatureVector &featureVector, const bool useRPhi, const bool useDirection) const
{
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_beamDeweighting));
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_energyKick));
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_globalAsymmetry));
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_localAsymmetry));
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_showerAsymmetry));

    if (useRPhi)
        featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_rPhiFeature));

    //if (m_enableIntersectionFeatures)
    //    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_intersectionFeature));

    if (useDirection)
    {
        featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_directionFeature));
        featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_directionBeginpointDistance));
        featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_directionFitQuality));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Vertex * SvmVertexSelectionAlgorithm::CompareVertices(const VertexVector &vertexVector, const VertexFeatureInfoMap &vertexFeatureInfoMap,
    const LArMvaHelper::MvaFeatureVector &eventFeatureList, const SupportVectorMachine &supportVectorMachine, const bool useRPhi, const bool useDirection) const
{
    const Vertex *pBestVertex(vertexVector.front());
    LArMvaHelper::MvaFeatureVector chosenFeatureList;

    VertexFeatureInfo chosenVertexFeatureInfo(vertexFeatureInfoMap.at(pBestVertex));
    this->AddVertexFeaturesToVector(chosenVertexFeatureInfo, chosenFeatureList, useRPhi, useDirection);

    for (const Vertex *const pVertex : vertexVector)
    {
        if (pVertex == pBestVertex)
            continue;

        LArMvaHelper::MvaFeatureVector featureList;
        VertexFeatureInfo vertexFeatureInfo(vertexFeatureInfoMap.at(pVertex));
        this->AddVertexFeaturesToVector(vertexFeatureInfo, featureList, useRPhi, useDirection);

        if (LArMvaHelper::Classify(supportVectorMachine, eventFeatureList, featureList, chosenFeatureList))
        {
            pBestVertex = pVertex;
            chosenFeatureList = featureList;
        }
    }

    return pBestVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Vertex * SvmVertexSelectionAlgorithm::ScoreVertices(const VertexVector &vertexVector, ClusterList &clusterList, const VertexFeatureInfoMap &vertexFeatureInfoMap,
    const LArMvaHelper::MvaFeatureVector &eventFeatureList, const SupportVectorMachine &supportVectorMachine, const bool useRPhi, const bool useDirection) const 
{
    const Vertex *pBestVertex(vertexVector.front());

    if (m_onlyCorrectDirection)
    {
        for (const Vertex *const pVertex : vertexVector)
        {
            pandora::CartesianVector vertexPosition(pVertex->GetPosition());
            std::function< TrackDirectionTool::DirectionFitObject (const pandora::Cluster*) > lambda = [this](const pandora::Cluster* const pCluster){ return this->m_pTrackDirectionTool->GetClusterDirection(pCluster); };
            float directionFlowProbability(m_pDirectionFlowProbabilityTool->GetDirectionFlowProbability(lambda, vertexPosition, clusterList)); 

            if (directionFlowProbability > 0.9)
            {
                pBestVertex = pVertex;
                break;
            }
        }
    }

    LArMvaHelper::MvaFeatureVector chosenFeatureList;

    VertexFeatureInfo chosenVertexFeatureInfo(vertexFeatureInfoMap.at(pBestVertex));
    this->AddVertexFeaturesToVector(chosenVertexFeatureInfo, chosenFeatureList, useRPhi, useDirection);

    if (m_onlyCorrectDirection)
        std::cout << "Only correct direction vertex candidates are being considered" << std::endl;

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    pandora::MCParticleVector trueNeutrinos;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);

    CartesianVector trueVertexPosition(trueNeutrinos.front()->GetVertex());

    //std::cout << "(" << trueVertexPosition.GetX() << ", " << trueVertexPosition.GetY() << ", " << trueVertexPosition.GetZ() << ")" << std::endl;

    //for (const Vertex *const pVertex : vertexVector)
    //    std::cout << "(" << pVertex->GetPosition().GetX() << ", " << pVertex->GetPosition().GetY() << ", " << pVertex->GetPosition().GetZ() << ")" << std::endl;

    for (const Vertex *const pVertex : vertexVector)
    {
        if (pVertex == pBestVertex)
            continue;

        LArMvaHelper::MvaFeatureVector featureList;
        VertexFeatureInfo vertexFeatureInfo(vertexFeatureInfoMap.at(pVertex));
        this->AddVertexFeaturesToVector(vertexFeatureInfo, featureList, useRPhi, useDirection);

        float newVertexProbability(LArMvaHelper::CalculateProbability(supportVectorMachine, eventFeatureList, featureList, chosenFeatureList));
        float bestVertexProbability(LArMvaHelper::CalculateProbability(supportVectorMachine, eventFeatureList, chosenFeatureList, featureList));

        pandora::CartesianVector newVertexPosition(pVertex->GetPosition()), bestVertexPosition(pBestVertex->GetPosition());

        std::function< TrackDirectionTool::DirectionFitObject (const pandora::Cluster*) > lambda = [this](const pandora::Cluster* const pCluster){ return this->m_pTrackDirectionTool->GetClusterDirection(pCluster); };

        float newDirectionFlowProbability(m_pDirectionFlowProbabilityTool->GetDirectionFlowProbability(lambda, newVertexPosition, clusterList)); 
        float bestDirectionFlowProbability(m_pDirectionFlowProbabilityTool->GetDirectionFlowProbability(lambda, bestVertexPosition, clusterList)); 

        std::cout << "Vertex DR: " << GetVertexDR(pBestVertex, true) << " with probability (P_svm * P_df): " << bestVertexProbability << " * " <<  bestDirectionFlowProbability << " = " << bestVertexProbability * bestDirectionFlowProbability << std::endl;
        std::cout << "Vertex DR: " << GetVertexDR(pVertex, true) << " with probability (P_svm * P_df): " << newVertexProbability << " * " <<  newDirectionFlowProbability << " = " << newVertexProbability * newDirectionFlowProbability << std::endl;
        std::cout << "----------------------------------------------------" << std::endl;

        if (m_onlyCorrectDirection && newDirectionFlowProbability < 0.1)
            continue;

        if (newVertexProbability * newDirectionFlowProbability > bestVertexProbability * bestDirectionFlowProbability)
        {
            pBestVertex = pVertex;
            chosenFeatureList = featureList;
        }
    }

    return pBestVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float SvmVertexSelectionAlgorithm::GetVertexDR(const pandora::Vertex* const pVertex, bool enableSpaceChargeCorrection) const
{
    // Extract input collections
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    pandora::MCParticleVector trueNeutrinos;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);

    if (trueNeutrinos.size() != 1)
        return -1.f;

    CartesianVector trueVertexPosition(trueNeutrinos.front()->GetVertex());
    CartesianVector vertexCandidatePosition(pVertex->GetPosition());
        
    if (enableSpaceChargeCorrection)
        vertexCandidatePosition = LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(vertexCandidatePosition);

    return (vertexCandidatePosition - trueVertexPosition).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::PopulateFinalVertexScoreList(const VertexFeatureInfoMap &vertexFeatureInfoMap, const Vertex *const pFavouriteVertex,
    const VertexVector &vertexVector, VertexScoreList &finalVertexScoreList) const
{
    if (pFavouriteVertex)
    {
        const CartesianVector vertexPosition(pFavouriteVertex->GetPosition());

        for (const Vertex *const pVertex : vertexVector)
        {
            if ((pVertex->GetPosition() - vertexPosition).GetMagnitude() < m_rPhiFineTuningRadius)
            {
                const float rPhiScore(vertexFeatureInfoMap.at(pVertex).m_rPhiFeature);
                finalVertexScoreList.emplace_back(pVertex, rPhiScore);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SvmVertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    LArSpaceChargeHelper::Configure("/usera/jjd49/pandora_direction/PandoraPFA/LArContent-origin/vertex_direction/larpandoracontent/LArDirection/SCEoffsets_MicroBooNE_E273.root");

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureTools", algorithmToolVector));

    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVector));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SvmFileName", m_svmFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RegionSvmName", m_regionSvmName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexSvmName", m_vertexSvmName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingSetMode", m_trainingSetMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AllowClassifyDuringTraining", m_allowClassifyDuringTraining));

    if ((!m_trainingSetMode || m_allowClassifyDuringTraining))
    {
        if (m_svmFileName.empty() || m_regionSvmName.empty() || m_vertexSvmName.empty())
        {
            std::cout << "SvmVertexSelectionAlgorithm: SvmFileName, RegionSvmName and VertexSvmName must be set if training set mode is" <<
                         "off or we allow classification during training" << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }

        const std::string fullSvmFileName(LArFileHelper::FindFileInPath(m_svmFileName, m_filePathEnvironmentVariable));
        m_svMachineRegion.Initialize(fullSvmFileName, m_regionSvmName);
        m_svMachineVertex.Initialize(fullSvmFileName, m_vertexSvmName);
    }


    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCVertexXCorrection", m_mcVertexXCorrection));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingOutputFileRegion", m_trainingOutputFileRegion));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingOutputFileVertex", m_trainingOutputFileVertex));

    if (m_trainingSetMode && (m_trainingOutputFileRegion.empty() || m_trainingOutputFileVertex.empty()))
    {
        std::cout << "SvmVertexSelectionAlgorithm: TrainingOutputFileRegion and TrainingOutputFileVertex are required for training set " <<
                     "mode" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CaloHitListName", m_caloHitListName));

    if (m_trainingSetMode && (m_mcParticleListName.empty() || m_caloHitListName.empty()))
    {
        std::cout << "SvmVertexSelectionAlgorithm: MCParticleListName and CaloHitListName are required for training set mode" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerSpineLength", m_minShowerSpineLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BeamDeweightingConstant", m_beamDeweightingConstant));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LocalAsymmetryConstant", m_localAsymmetryConstant));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GlobalAsymmetryConstant", m_globalAsymmetryConstant));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerAsymmetryConstant", m_showerAsymmetryConstant));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnergyKickConstant", m_energyKickConstant));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerClusteringDistance", m_showerClusteringDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerClusterHits", m_minShowerClusterHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseShowerClusteringApproximation", m_useShowerClusteringApproximation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RegionRadius", m_regionRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RPhiFineTuningRadius", m_rPhiFineTuningRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTrueVertexRadius", m_maxTrueVertexRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseRPhiFeatureForRegion", m_useRPhiFeatureForRegion));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DropFailedRPhiFastScoreCandidates", m_dropFailedRPhiFastScoreCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableIntersectionFeatures", m_enableIntersectionFeatures));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableDirectionFeatures", m_enableDirectionFeatures));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DirectionScoreReweighting", m_directionScoreReweighting));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OnlyCorrectDirection", m_onlyCorrectDirection));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ExtrapolationNumberSteps", m_extrapolationNSteps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ExtrapolationStepSize", m_extrapolationStepSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualiseEvent", m_visualiseEvent));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FileIdentifier", m_fileIdentifier));

    AlgorithmTool *pAlgorithmTool(nullptr);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "TrackDirection", pAlgorithmTool));

    if (!(this->m_pTrackDirectionTool = dynamic_cast<TrackDirectionTool *>(pAlgorithmTool)))
        throw STATUS_CODE_FAILURE;

    AlgorithmTool *pAnotherAlgorithmTool(nullptr);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "DirectionFlowProbability", pAnotherAlgorithmTool));

    if (!(this->m_pDirectionFlowProbabilityTool = dynamic_cast<DirectionFlowProbabilityTool*>(pAnotherAlgorithmTool)))
        throw STATUS_CODE_FAILURE;

    return VertexSelectionBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
