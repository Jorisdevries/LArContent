/**
 *  @file   larpandoracontent/LArTwoDReco/ClusterSplitting/ClusterSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the cluster splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ClusterSplittingAlgorithm::Run()
{
    if (m_inputClusterListNames.empty())
        return this->RunUsingCurrentList();

    std::string originalListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, originalListName));

    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const StatusCode listChangeStatusCode(PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

        if (STATUS_CODE_NOT_FOUND == listChangeStatusCode)
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ClusterSplittingAlgorithm: cluster list not found " << clusterListName << std::endl;

            continue;
        }

        if (STATUS_CODE_SUCCESS != listChangeStatusCode)
            return listChangeStatusCode;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunUsingCurrentList());
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, originalListName));
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterSplittingAlgorithm::RunUsingCurrentList() const
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    ClusterList internalClusterList(pClusterList->begin(), pClusterList->end());
    internalClusterList.sort(LArClusterHelper::SortByNHits);

    for (ClusterList::iterator iter = internalClusterList.begin(); iter != internalClusterList.end(); ++iter)
    {
        const Cluster *const pCluster = *iter;
        ClusterList clusterSplittingList;

        if (STATUS_CODE_SUCCESS != this->SplitCluster(pCluster, clusterSplittingList))
            continue;

        internalClusterList.splice(internalClusterList.end(), clusterSplittingList);
        *iter = NULL;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterSplittingAlgorithm::SplitCluster(const Cluster *const pCluster, ClusterList &clusterSplittingList) const
{
    // Split cluster into two CaloHit lists
    PandoraContentApi::Cluster::Parameters firstParameters, secondParameters;

    if (STATUS_CODE_SUCCESS != this->DivideCaloHits(pCluster, firstParameters.m_caloHitList, secondParameters.m_caloHitList))
        return STATUS_CODE_NOT_FOUND;

    if (firstParameters.m_caloHitList.empty() || secondParameters.m_caloHitList.empty())
        return STATUS_CODE_NOT_ALLOWED;

    // Begin cluster fragmentation operations
    const ClusterList clusterList(1, pCluster);
    std::string clusterListToSaveName, clusterListToDeleteName;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName,
        clusterListToSaveName));

    // Create new clusters
    const Cluster *pFirstCluster(NULL), *pSecondCluster(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, firstParameters, pFirstCluster));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, secondParameters, pSecondCluster));

    clusterSplittingList.push_back(pFirstCluster);
    clusterSplittingList.push_back(pSecondCluster);

    // End cluster fragmentation operations
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNames", m_inputClusterListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
