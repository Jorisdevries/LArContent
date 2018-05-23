/**
 *  @file   larpandoracontent/LArVertex/DirectionFlowProbabilityAlgorithm.h
 *
 *  @brief  Header file for the svm vertex selection algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DIRECTION_FLOW_PROBABILITY_ALGORITHM_H
#define LAR_DIRECTION_FLOW_PROBABILITY_ALGORITHM_H 1

#include "Api/PandoraContentApi.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArDirection/TrackDirectionTool.h"
#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

namespace lar_content
{
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DirectionFlowProbabilityAlgorithm class
 */
class DirectionFlowProbabilityAlgorithm : public VertexSelectionBaseAlgorithm
{
public:
    /**
     *  @brief Emerging cluster info class
     */
    class EmergingCluster
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  beamDeweighting the beam deweighting feature
         *  @param  rPhiFeature the r/phi feature
         */
        EmergingCluster(const pandora::CartesianVector originPosition, const EmergingCluster* pEmergingCluster, TrackDirectionTool::DirectionFitObject* pDirectionFit);

        /**
         *  @brief  Constructor
         */
        const EmergingCluster* ParentEmergingCluster() const;

        /**
         *  @brief  Constructor
         */
        const pandora::CartesianVector* Origin() const;

        /**
         *  @brief  Constructor
         */
        TrackDirectionTool::DirectionFitObject* DirectionFit() const;

        const pandora::CartesianVector*             m_origin;           ///< The point of origin.
        const EmergingCluster*                      m_parent;           ///< The parent EmergingCluster 
        //const pandora::Cluster*                     m_cluster;          ///< The underlying cluster 
        TrackDirectionTool::DirectionFitObject*     m_directionfit;     ///< Direction fit.
    };

    typedef std::vector<EmergingCluster*> EmergingClusterVector;


    //--------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Default constructor
     */
    DirectionFlowProbabilityAlgorithm();

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    typedef std::unordered_map<const pandora::Cluster*, pandora::CartesianPointVector> ClusterToSpacepointsMap; 
    typedef std::map<const pandora::Cluster*, TrackDirectionTool::DirectionFitObject*> DirectionFitMap;

    void GetVertexScoreList();

    void Test();

    float GetDirectionFlowProbability(const pandora::CartesianVector &positionVector, const pandora::ClusterList *pClusterList);

    /** 
     *  @brief  Select a subset of input clusters (contained in the input list names) for processing in this algorithm
     *
     *  @param  pClusterList input cluster list 
     *  @param  selectedClusterVector to receive the selected clusters 
     */
    void SelectClusters(const pandora::ClusterList *pClusterList, pandora::ClusterVector &selectedClusterVector);


    bool BeginpointIsParentOrigin(const EmergingCluster* const pEmergingCluster, const EmergingCluster* const pParentEmergingCluster) const;

    EmergingClusterVector CreateEmergingClusters(const pandora::CartesianVector &vertexPosition, const pandora::ClusterVector &inputClusterVector);

    pandora::ClusterVector GetPrimaryClusters(const pandora::CartesianVector &positionVector, const pandora::ClusterVector &inputClusterVector) const;

    /**
     *  @brief  Get the vertex score list
     *
     *  @param  vertexVector the vector of vertices
     */
    pandora::ClusterVector GetOrderedDaughters(const pandora::CartesianVector &positionVector, const pandora::Cluster* const pParentCluster, const pandora::ClusterVector &inputClusterVector) const;

    bool ClusterPointsToPosition(const pandora::Cluster *const pCluster, const pandora::CartesianVector &positionVector, ClusterToSpacepointsMap &clusterToSpacepointsMap) const;

    /** 
     *  @brief  Identify where (extrapolated) clusters plausibly cross in 2D
     *
     *  @param  clusterVector the input clusters
     *  @param  crossingPoints to receive the 2D crossing points
     */
    ClusterToSpacepointsMap FillClusterToSpacepointsMap(const pandora::ClusterVector &clusterVector) const;

    /** 
     *  @brief  Get a list of spacepoints representing cluster 2D hit positions and extrapolated positions
     *
     *  @param  pCluster address of the cluster
     *  @param  spacePoints to receive the list of spacepoints
     */
    void GetSpacepoints(const pandora::Cluster *const pCluster, pandora::CartesianPointVector &spacePoints) const;

    void AddToSlidingFitCache(const pandora::Cluster *const pCluster);

    const TwoDSlidingFitResult &GetCachedSlidingFit(const pandora::Cluster *const pCluster) const;

    TrackDirectionTool::DirectionFitObject* GetCachedDirectionFit(const pandora::Cluster *const pCluster);

    unsigned int            m_slidingFitWindow;                 ///< The layer window for the sliding linear fits 
    TwoDSlidingFitResultMap m_slidingFitResultMap;              ///< The sliding fit result map
    DirectionFitMap         m_directionFitMap;                  ///< The direction fit map
    float                   m_impactRadius;                     ///< The impact radius determining whether a sliding fit extrapolation points to a position 
    unsigned int            m_extrapolationNSteps;              ///< The number of steps used in the sliding fit extrapolation method
    float                   m_extrapolationStepSize;            ///< The extrapolation step size.
    float                   m_minimumClusterLength;             ///< The minimum length a cluster must be in order to be considered 
    TrackDirectionTool      *m_pTrackDirectionTool;             ///< The track direction tool 
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const DirectionFlowProbabilityAlgorithm::EmergingCluster* DirectionFlowProbabilityAlgorithm::EmergingCluster::ParentEmergingCluster() const
{
    return m_parent;
} 

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector* DirectionFlowProbabilityAlgorithm::EmergingCluster::Origin() const
{
    return m_origin;
} 

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::DirectionFitObject* DirectionFlowProbabilityAlgorithm::EmergingCluster::DirectionFit() const
{
    return m_directionfit;
} 

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // #ifndef LAR_DIRECTION_FLOW_PROBABILITY_ALGORITHM_H
