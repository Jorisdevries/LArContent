/**
 *  @file   larpandoracontent/LArControlFlow/MasterAlgorithm.cc
 *
 *  @brief  Implementation of the master algorithm class.
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArContent.h"

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArDirectionHelper.h"
#include "larpandoracontent/LArHelpers/LArSpaceChargeHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#include "larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MasterAlgorithm::MasterAlgorithm() :
    m_shouldRunAllHitsCosmicReco(true),
    m_shouldRunStitching(true),
    m_shouldRunCosmicHitRemoval(true),
    m_cheatTopFaceCosmicRemoval(false),
    m_recoTopFaceCosmicRemoval(false),
    m_shouldRunSlicing(true),
    m_shouldRunNeutrinoRecoOption(true),
    m_shouldRunCosmicRecoOption(true),
    m_shouldPerformSliceId(true),
    m_printOverallRecoStatus(false),
    m_visualizeOverallRecoStatus(false),
    m_pSlicingWorkerInstance(nullptr),
    m_pSliceNuWorkerInstance(nullptr),
    m_pSliceCRWorkerInstance(nullptr),
    m_fullWidthCRWorkerWireGaps(true),
    m_passMCParticlesToWorkerInstances(true),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH"),
    m_inTimeMaxX0(1.f),
    m_writeToTree(false)
{
}

    bool eventContainsTargetPfo(false);
    bool eventContainsTrueTargetCR(false);

//------------------------------------------------------------------------------------------------------------------------------------------

MasterAlgorithm::~MasterAlgorithm()
{
    if (m_writeToTree)
    {    
        try  
        {    
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "TopFaceCosmicRemoval", m_fileName.c_str(), "UPDATE"));
        }    
        catch (const StatusCodeException &)
        {    
            std::cout << "EventValidationAlgorithm: Unable to write tree " << m_treeName << " to file " << m_fileName << std::endl;
        }    
    }    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MasterAlgorithm::ShiftPfoHierarchy(const ParticleFlowObject *const pParentPfo, const PfoToLArTPCMap &pfoToLArTPCMap, const float x0) const
{
    if (!pParentPfo->GetParentPfoList().empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    PfoToLArTPCMap::const_iterator larTPCIter(pfoToLArTPCMap.find(pParentPfo));

    if (pfoToLArTPCMap.end() == larTPCIter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    PfoList pfoList;
    LArPfoHelper::GetAllDownstreamPfos(pParentPfo, pfoList);
    const float signedX0(larTPCIter->second->IsDriftInPositiveX() ? -x0 : x0);

    if (m_visualizeOverallRecoStatus)
    {
        std::cout << "ShiftPfoHierarchy: signedX0 " << signedX0 << std::endl;
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &pfoList, "BeforeShiftCRPfos", GREEN));
    }

    for (const ParticleFlowObject *const pDaughterPfo : pfoList)
    {
        CaloHitList caloHitList;
        for (const Cluster *const pCluster : pDaughterPfo->GetClusterList())
        {
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
            caloHitList.insert(caloHitList.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
        }

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            PandoraContentApi::CaloHit::Metadata metadata;
            metadata.m_x0 = signedX0;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::AlterMetadata(*this, pCaloHit, metadata));
        }

        for (const Vertex *const pVertex : pDaughterPfo->GetVertexList())
        {
            PandoraContentApi::Vertex::Metadata metadata;
            metadata.m_x0 = signedX0;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::AlterMetadata(*this, pVertex, metadata));
        }
    }

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &pfoList, "AfterShiftCRPfos", RED));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MasterAlgorithm::StitchPfos(const ParticleFlowObject *const pPfoToEnlarge, const ParticleFlowObject *const pPfoToDelete,
    PfoToLArTPCMap &pfoToLArTPCMap) const
{
    if (pPfoToEnlarge == pPfoToDelete)
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    // ATTN Remove pfos from pfo to lar tpc map here, avoiding problems if stitching across multiple tpcs.
    pfoToLArTPCMap.erase(pPfoToEnlarge);
    pfoToLArTPCMap.erase(pPfoToDelete);

    const PfoList daughterPfos(pPfoToDelete->GetDaughterPfoList());
    const ClusterVector daughterClusters(pPfoToDelete->GetClusterList().begin(), pPfoToDelete->GetClusterList().end());
    const VertexVector daughterVertices(pPfoToDelete->GetVertexList().begin(), pPfoToDelete->GetVertexList().end());

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete, m_recreatedPfoListName));

    for (const ParticleFlowObject *const pDaughterPfo : daughterPfos)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pPfoToEnlarge, pDaughterPfo));

    for (const  Vertex *const pDaughterVertex : daughterVertices)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pDaughterVertex, m_recreatedVertexListName));

    for (const Cluster *const pDaughterCluster : daughterClusters)
    {
        const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
        const Cluster *pParentCluster(PfoMopUpBaseAlgorithm::GetParentCluster(pPfoToEnlarge->GetClusterList(), daughterHitType));

        if (pParentCluster)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pDaughterCluster,
                m_recreatedClusterListName, m_recreatedClusterListName));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfoToEnlarge, pDaughterCluster));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Initialize()
{
    try
    {
        const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
        const DetectorGapList &gapList(this->GetPandora().GetGeometry()->GetDetectorGapList());

        for (const LArTPCMap::value_type &mapEntry : larTPCMap)
            m_crWorkerInstances.push_back(this->CreateWorkerInstance(*(mapEntry.second), gapList, m_crSettingsFile));

        if (m_shouldRunSlicing)
            m_pSlicingWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_slicingSettingsFile);

        if (m_shouldRunNeutrinoRecoOption)
            m_pSliceNuWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_nuSettingsFile);

        if (m_shouldRunCosmicRecoOption)
            m_pSliceCRWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_crSettingsFile);
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cout << "MasterAlgorithm: Exception during initialization of worker instances " << statusCodeException.ToString() << std::endl;
        return statusCodeException.GetStatusCode();
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Run()
{
    ++m_eventNumber;

    if (m_passMCParticlesToWorkerInstances)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CopyMCParticles());

    PfoToFloatMap stitchedPfosToX0Map;
    VolumeIdToHitListMap volumeIdToHitListMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetVolumeIdToHitListMap(volumeIdToHitListMap));

    if (m_shouldRunAllHitsCosmicReco)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunCosmicRayReconstruction(volumeIdToHitListMap));

        PfoToLArTPCMap pfoToLArTPCMap;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RecreateCosmicRayPfos(pfoToLArTPCMap));

        if (m_shouldRunStitching)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->StitchCosmicRayPfos(pfoToLArTPCMap, stitchedPfosToX0Map));
    }

    if (m_shouldRunCosmicHitRemoval)
    {
        PfoList clearCosmicRayPfos, ambiguousPfos;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->TagCosmicRayPfos(stitchedPfosToX0Map, clearCosmicRayPfos, ambiguousPfos));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunCosmicRayHitRemoval(ambiguousPfos));
    }

    SliceVector sliceVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunSlicing(volumeIdToHitListMap, sliceVector));

    if (m_shouldRunNeutrinoRecoOption || m_shouldRunCosmicRecoOption)
    {
        SliceHypotheses nuSliceHypotheses, crSliceHypotheses;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunSliceReconstruction(sliceVector, nuSliceHypotheses, crSliceHypotheses));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SelectBestSliceHypotheses(nuSliceHypotheses, crSliceHypotheses));
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Reset());
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::CopyMCParticles() const
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputMCParticleListName, pMCParticleList));

    PandoraInstanceList pandoraWorkerInstances(m_crWorkerInstances);
    if (m_pSlicingWorkerInstance) pandoraWorkerInstances.push_back(m_pSlicingWorkerInstance);
    if (m_pSliceNuWorkerInstance) pandoraWorkerInstances.push_back(m_pSliceNuWorkerInstance);
    if (m_pSliceCRWorkerInstance) pandoraWorkerInstances.push_back(m_pSliceCRWorkerInstance);

    LArMCParticleFactory mcParticleFactory;

    for (const Pandora *const pPandoraWorker : pandoraWorkerInstances)
    {
        for (const MCParticle *const pMCParticle : *pMCParticleList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(pPandoraWorker, pMCParticle, &mcParticleFactory));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::GetVolumeIdToHitListMap(VolumeIdToHitListMap &volumeIdToHitListMap) const
{
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const unsigned int nLArTPCs(larTPCMap.size());

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputHitListName, pCaloHitList));

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit*>(pCaloHit));

        if (!pLArCaloHit && (1 != nLArTPCs))
            return STATUS_CODE_INVALID_PARAMETER;

        const unsigned int volumeId(pLArCaloHit ? pLArCaloHit->GetLArTPCVolumeId() : 0);
        const LArTPC *const pLArTPC(larTPCMap.at(volumeId));

        LArTPCHitList &larTPCHitList(volumeIdToHitListMap[volumeId]);
        larTPCHitList.m_allHitList.push_back(pCaloHit);

        if (((pCaloHit->GetPositionVector().GetX() >= (pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX())) &&
            (pCaloHit->GetPositionVector().GetX() <= (pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX()))))
        {
            larTPCHitList.m_truncatedHitList.push_back(pCaloHit);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunCosmicRayReconstruction(const VolumeIdToHitListMap &volumeIdToHitListMap) const
{
    unsigned int workerCounter(0);

    for (const Pandora *const pCRWorker : m_crWorkerInstances)
    {
        const LArTPC &larTPC(pCRWorker->GetGeometry()->GetLArTPC());
        VolumeIdToHitListMap::const_iterator iter(volumeIdToHitListMap.find(larTPC.GetLArTPCVolumeId()));

        if (volumeIdToHitListMap.end() == iter)
            continue;

        for (const CaloHit *const pCaloHit : iter->second.m_allHitList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(pCRWorker, pCaloHit));

        if (m_printOverallRecoStatus)
            std::cout << "Running cosmic-ray reconstruction worker instance " << ++workerCounter << " of " << m_crWorkerInstances.size() << std::endl;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pCRWorker));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RecreateCosmicRayPfos(PfoToLArTPCMap &pfoToLArTPCMap) const
{
    for (const Pandora *const pCRWorker : m_crWorkerInstances)
    {
        const PfoList *pCRPfos(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pCRWorker, pCRPfos));

        PfoList newPfoList;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Recreate(*pCRPfos, newPfoList));

        const LArTPC &larTPC(pCRWorker->GetGeometry()->GetLArTPC());

        for (const Pfo *const pNewPfo : newPfoList)
            pfoToLArTPCMap[pNewPfo] = &larTPC;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::StitchCosmicRayPfos(PfoToLArTPCMap &pfoToLArTPCMap, PfoToFloatMap &stitchedPfosToX0Map) const
{
    const PfoList *pRecreatedCRPfos(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(this->GetPandora(), pRecreatedCRPfos));

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pRecreatedCRPfos, "RecreatedCRPfos", GREEN));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    for (StitchingBaseTool *const pStitchingTool : m_stitchingToolVector)
        pStitchingTool->Run(this, pRecreatedCRPfos, pfoToLArTPCMap, stitchedPfosToX0Map);

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pRecreatedCRPfos, "AfterStitchingCRPfos", RED));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::TagCosmicRayPfos(const PfoToFloatMap &stitchedPfosToX0Map, PfoList &clearCosmicRayPfos, PfoList &ambiguousPfos) const
{
    const PfoList *pRecreatedCRPfos(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(this->GetPandora(), pRecreatedCRPfos));

    PfoList nonStitchedParentCosmicRayPfos;
    for (const Pfo *const pPfo : *pRecreatedCRPfos)
    {
        if (!pPfo->GetParentPfoList().empty())
            continue;

        float pfoCosmicProbability(0.f), deltaChiSquaredUpDownPerHit(0.f), deltaChiSquaredUpDown(0.f), pfoPolarAngle(0.f);
        bool isTargetTopFacePfo(false), pfoIntersectsTopFace(false), pfoFiducialLowY(false), pfoHighTopY(false), isTargetCosmic(false);
        
        if (m_recoTopFaceCosmicRemoval || m_writeToTree)
        {
            const auto pMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));
            const bool isTrueCosmicRay(LArMCParticleHelper::IsCosmicRay(pMCParticle));

            const pandora::CartesianVector lowYVector(pMCParticle->GetVertex().GetY() < pMCParticle->GetEndpoint().GetY() ? pMCParticle->GetVertex() : pMCParticle->GetEndpoint());
            const pandora::CartesianVector highYVector(pMCParticle->GetVertex().GetY() > pMCParticle->GetEndpoint().GetY() ? pMCParticle->GetVertex() : pMCParticle->GetEndpoint());

            const bool mcFiducialLowY(LArDirectionHelper::IsInFiducialVolume(lowYVector));
            const bool mcIntersectsTopFace(highYVector.GetZ() > 116.5 && LArDirectionHelper::IntersectsYFace(lowYVector, highYVector));
            isTargetCosmic = (isTrueCosmicRay && mcFiducialLowY && mcIntersectsTopFace);

            try
            {
                CartesianVector yAxis(0.f, 1.f, 0.f);
                pfoPolarAngle = LArDirectionHelper::GetAngleWithVector(this->GetPandora(), pPfo, yAxis);

                TrackDirectionTool::DirectionFitObject directionFit = m_pTrackDirectionTool->GetPfoDirection(pPfo);
                float fitEndpointEnergy(directionFit.GetFitParameters().GetParameterZero());
                float endpointFitChargeRatio(directionFit.GetEndpointFitChargeRatio(5));

                pfoCosmicProbability = LArDirectionHelper::CalculateCosmicProbability(directionFit);
                deltaChiSquaredUpDownPerHit = directionFit.GetUpDownDeltaChiSquaredPerHit();
                deltaChiSquaredUpDown = directionFit.GetUpDownDeltaChiSquared();

                std::vector<pandora::CartesianVector> lowHighYPositions(LArDirectionHelper::GetLowHighYPoints(this->GetPandora(), pPfo));
                pandora::CartesianVector correctedPfoLowYVector(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(lowHighYPositions.at(0)));
                pandora::CartesianVector correctedPfoHighYVector(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(lowHighYPositions.at(1)));

                pfoIntersectsTopFace = LArDirectionHelper::IntersectsYFace(correctedPfoLowYVector, correctedPfoHighYVector);
                pfoFiducialLowY = LArDirectionHelper::IsInFiducialVolume(correctedPfoLowYVector);
                pfoHighTopY = (correctedPfoHighYVector.GetY() >= 100.0 ? true : false); 
                isTargetTopFacePfo = (pfoHighTopY && pfoIntersectsTopFace && pfoFiducialLowY);

                pandora::HitType view3D(TPC_3D), viewW(TPC_VIEW_W);
                pandora::CaloHitList hits3D, hitsW;
                LArPfoHelper::GetCaloHits(pPfo, view3D, hits3D);
                LArPfoHelper::GetCaloHits(pPfo, viewW, hits3D);

                if (pfoIntersectsTopFace && pfoFiducialLowY && pfoHighTopY)
                    eventContainsTargetPfo = true;
                if (isTargetCosmic)
                    eventContainsTrueTargetCR = true;

                //whether this particular cosmic would have been removed based on out-of-time information
                PfoToFloatMap::const_iterator pfoToX0Iter = stitchedPfosToX0Map.find(pPfo);
                const float x0Shift((pfoToX0Iter != stitchedPfosToX0Map.end()) ? pfoToX0Iter->second : 0.f);
                bool removedByRegularTagging(std::fabs(x0Shift) > m_inTimeMaxX0);

                /*
                //event 1185:2 is an example where endpoint filtering helps
                std::cout << ">>>>>>>>>>>>> deltaChiSquaredUpDownPerHit: " << deltaChiSquaredUpDownPerHit << std::endl;
                std::cout << "MC vertex x: " << pMCParticle->GetVertex().GetX() << std::endl;

                if (pMCParticle->GetVertex().GetX() > 19.2 && pMCParticle->GetVertex().GetX() < 19.3)
                {
                    std::cout << "MC vertex x: " << pMCParticle->GetVertex().GetX() << std::endl;
                    directionFit.DrawFit();
                }
                */

                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "SufficientThreeDHits", hits3D.size() >= 50 ? 1 : 0));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "NumberThreeDHits", static_cast<int>(hits3D.size())));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "NumberWHits", static_cast<int>(hitsW.size())));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "RemovedByRegularTagging", removedByRegularTagging ? 1 : 0));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "MCPDG", pMCParticle->GetParticleId()));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "FitEndpointEnergy", fitEndpointEnergy));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "EndpointFitChargeRatio", endpointFitChargeRatio));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "MCTargetCosmic", isTargetCosmic ? 1 : 0));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "MCCosmicRay", isTrueCosmicRay ? 1 : 0));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "MCFiducialLowY", mcFiducialLowY? 1 : 0));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "MCStoppingTopFace", LArDirectionHelper::IsStoppingTopFaceMCParticle(pMCParticle) ? 1 : 0));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "MCIntersectsTopFace", mcIntersectsTopFace ? 1 : 0));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "MCLowY", lowYVector.GetY()));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "MCHighY", highYVector.GetY()));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "TargetPFO", isTargetTopFacePfo ? 1 : 0));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "PfoIntersectsTopFace", pfoIntersectsTopFace ? 1 : 0));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "PfoFiducialLowY", pfoFiducialLowY ? 1 : 0));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "PfoHighTopY", pfoHighTopY ? 1 : 0));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "PfoSuperHighTopY", correctedPfoHighYVector.GetY() >= 115.5 ? 1 : 0));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "PfoPolarAngle", pfoPolarAngle));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "PfoPolarAngleCriteriaMet", (pfoPolarAngle < 3.1415/2.0 && correctedPfoHighYVector.GetZ() > correctedPfoLowYVector.GetZ()) ? 1 : 0));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "HighYZLargerThanLowYZ", correctedPfoHighYVector.GetZ() > correctedPfoLowYVector.GetZ() ? 1 : 0));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "PfoLowY", correctedPfoLowYVector.GetY()));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "PfoHighY", correctedPfoHighYVector.GetY()));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "PfoLowYEndpointZ", correctedPfoLowYVector.GetZ()));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "PfoHighYEndpointZ", correctedPfoHighYVector.GetZ()));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "TruePfoPolarAngle", (pMCParticle->GetEndpoint() - pMCParticle->GetVertex()).GetOpeningAngle(yAxis)));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "TrueExtent", (pMCParticle->GetEndpoint() - pMCParticle->GetVertex()).GetMagnitude()));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "TrueDownwards", pMCParticle->GetEndpoint().GetY() < pMCParticle->GetVertex().GetY() ? 1 : 0 ));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "CosmicProbability", pfoCosmicProbability));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "DeltaChiSquaredUpDownPerHit", deltaChiSquaredUpDownPerHit));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "DeltaChiSquaredUpDown", deltaChiSquaredUpDown));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "FileIdentifier", m_fileIdentifier));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TopFaceCosmicRemoval", "EventNumber", m_eventNumber));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), "TopFaceCosmicRemoval"));

                //std::cout << "No exception 1" << std::endl;
            }
            catch (...)
            {
                //std::cout << "Exception 1" << std::endl;
                //continue;
            }
        }

        PfoToFloatMap::const_iterator pfoToX0Iter = stitchedPfosToX0Map.find(pPfo);
        const float x0Shift((pfoToX0Iter != stitchedPfosToX0Map.end()) ? pfoToX0Iter->second : 0.f);
        PfoList &targetList(((std::fabs(x0Shift) > m_inTimeMaxX0) || (m_cheatTopFaceCosmicRemoval && isTargetCosmic) || (m_recoTopFaceCosmicRemoval && isTargetTopFacePfo && deltaChiSquaredUpDownPerHit <= -0.25)) ? clearCosmicRayPfos : nonStitchedParentCosmicRayPfos);
        targetList.push_back(pPfo);

    }

    for (CosmicRayTaggingBaseTool *const pCosmicRayTaggingTool : m_cosmicRayTaggingToolVector)
        pCosmicRayTaggingTool->FindAmbiguousPfos(nonStitchedParentCosmicRayPfos, ambiguousPfos, this);

    for (const Pfo *const pPfo : nonStitchedParentCosmicRayPfos)
    {
        if (ambiguousPfos.end() == std::find(ambiguousPfos.begin(), ambiguousPfos.end(), pPfo))
            clearCosmicRayPfos.push_back(pPfo);
    }

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &clearCosmicRayPfos, "ClearCRPfos", RED));
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &ambiguousPfos, "AmbiguousCRPfos", BLUE));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunCosmicRayHitRemoval(const PfoList &ambiguousPfos) const
{
    PfoList allPfosToDelete;
    LArPfoHelper::GetAllConnectedPfos(ambiguousPfos, allPfosToDelete);

    for (const Pfo *const pPfoToDelete : allPfosToDelete)
    {
        const ClusterList clusterList(pPfoToDelete->GetClusterList());
        const VertexList vertexList(pPfoToDelete->GetVertexList());
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &clusterList));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &vertexList));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunSlicing(const VolumeIdToHitListMap &volumeIdToHitListMap, SliceVector &sliceVector) const
{
    for (const VolumeIdToHitListMap::value_type &mapEntry : volumeIdToHitListMap)
    {
        for (const CaloHit *const pCaloHit : mapEntry.second.m_truncatedHitList)
        {
            if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
                continue;

            const HitType hitType(pCaloHit->GetHitType());
            if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
                continue;

            if (m_shouldRunSlicing)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSlicingWorkerInstance, pCaloHit));
            }
            else
            {
                if (sliceVector.empty()) sliceVector.push_back(CaloHitList());
                sliceVector.back().push_back(pCaloHit);
            }
        }
    }

    if (m_shouldRunSlicing)
    {
        if (m_printOverallRecoStatus)
            std::cout << "Running slicing worker instance" << std::endl;

        const PfoList *pSlicePfos(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSlicingWorkerInstance));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSlicingWorkerInstance, pSlicePfos));

        if (m_visualizeOverallRecoStatus)
        {
            PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pSlicePfos, "OnePfoPerSlice", BLUE));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }

        for (const Pfo *const pSlicePfo : *pSlicePfos)
        {
            sliceVector.push_back(CaloHitList());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_U, sliceVector.back());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_V, sliceVector.back());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_W, sliceVector.back());
        }
    }

    if (m_printOverallRecoStatus)
        std::cout << "Identified " << sliceVector.size() << " slice(s)" << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunSliceReconstruction(SliceVector &sliceVector, SliceHypotheses &nuSliceHypotheses, SliceHypotheses &crSliceHypotheses) const
{
    unsigned int sliceCounter(0);

    for (const CaloHitList &sliceHits : sliceVector)
    {
        for (const CaloHit *const pSliceCaloHit : sliceHits)
        {
            // ATTN Must ensure we copy the hit actually owned by master instance; access differs with/without slicing enabled
            const CaloHit *const pCaloHitInMaster(m_shouldRunSlicing ? static_cast<const CaloHit*>(pSliceCaloHit->GetParentAddress()) : pSliceCaloHit);

            if (m_shouldRunNeutrinoRecoOption)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSliceNuWorkerInstance, pCaloHitInMaster));

            if (m_shouldRunCosmicRecoOption)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSliceCRWorkerInstance, pCaloHitInMaster));
        }

        ++sliceCounter;

        if (m_shouldRunNeutrinoRecoOption)
        {
            if (m_printOverallRecoStatus)
                std::cout << "Running nu worker instance for slice " << sliceCounter << " of " << sliceVector.size() << std::endl;

            const PfoList *pSliceNuPfos(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSliceNuWorkerInstance));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSliceNuWorkerInstance, pSliceNuPfos));
            nuSliceHypotheses.push_back(*pSliceNuPfos);
        }

        if (m_shouldRunCosmicRecoOption)
        {
            if (m_printOverallRecoStatus)
                std::cout << "Running cr worker instance for slice " << sliceCounter << " of " << sliceVector.size() << std::endl;

            const PfoList *pSliceCRPfos(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSliceCRWorkerInstance));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSliceCRWorkerInstance, pSliceCRPfos));
            crSliceHypotheses.push_back(*pSliceCRPfos);
        }
    }

    if (m_shouldRunNeutrinoRecoOption && m_shouldRunCosmicRecoOption && (nuSliceHypotheses.size() != crSliceHypotheses.size()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::SelectBestSliceHypotheses(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses) const
{
    if (m_printOverallRecoStatus)
        std::cout << "Select best slice hypotheses" << std::endl;

    PfoList selectedSlicePfos;

    if (m_shouldPerformSliceId)
    {
        for (SliceIdBaseTool *const pSliceIdTool : m_sliceIdToolVector)
            pSliceIdTool->SelectOutputPfos(this, nuSliceHypotheses, crSliceHypotheses, selectedSlicePfos);
    }
    else if (m_shouldRunNeutrinoRecoOption != m_shouldRunCosmicRecoOption)
    {
        const SliceHypotheses &sliceHypotheses(m_shouldRunNeutrinoRecoOption ? nuSliceHypotheses : crSliceHypotheses);

        for (const PfoList &slice : sliceHypotheses)
            selectedSlicePfos.insert(selectedSlicePfos.end(), slice.begin(), slice.end());
    }

    /*
    //NEW
    PfoVector nuSlicePfos(nuSliceHypotheses.front().begin(), nuSliceHypotheses.front().end());
    std::sort(nuSlicePfos.begin(), nuSlicePfos.end(), [](const pandora::ParticleFlowObject* pPfo1, const pandora::ParticleFlowObject* pPfo2){return LArPfoHelper::GetThreeDLengthSquared(pPfo1) > LArPfoHelper::GetThreeDLengthSquared(pPfo2);}); //sort in descending order
    const auto pNuPfo(nuSlicePfos.front());

    std::cout << ">>>>>>>>>>>>>>>>>> " << LArPfoHelper::GetVertex(selectedSlicePfos.front())->GetPosition().GetX() << std::endl;

    bool nuRecoIsTargetPfo(false), nuRecoIsTrueTargetCR(false), nuRecoTrueCosmicRay(false), nuRecoIntersectsTopFace(false), nuRecoFiducialLowY(false), nuRecoHighTopY(false);
    float nuRecoDeltaChiSquared(0.f), nuRecoCosmicProbability(0.f), nuRecoPolarAngle(0.f);

    const PfoList *pRecreatedCRPfos(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(this->GetPandora(), pRecreatedCRPfos));

    PfoList nonStitchedParentCosmicRayPfos;
    for (const Pfo *const pPfo : *pRecreatedCRPfos)
    {
        if (!pPfo->GetParentPfoList().empty())
            continue;

        try
        {
            if ( LArPfoHelper::GetVertex(selectedSlicePfos.front())->GetPosition().GetX() ==  LArPfoHelper::GetVertex(pPfo)->GetPosition().GetX())
                std::cout << "FOUND NU RECO" << std::endl;

            const auto pMCParticle(LArMCParticleHelper::GetMainMCParticle(pNuPfo));
            const bool isTrueCosmicRay(LArMCParticleHelper::IsCosmicRay(pMCParticle));
            const bool isTargetCosmic(isTrueCosmicRay && this->IsStoppingTopFaceMCParticle(pMCParticle));

            CartesianVector yAxis(0.f, 1.f, 0.f);
            float pfoPolarAngle(GetAngleWithVector(pNuPfo, yAxis));

            TrackDirectionTool::DirectionFitObject directionFit = m_pTrackDirectionTool->GetPfoDirection(pNuPfo);
            float pfoCosmicProbability = this->CalculateCosmicProbability(directionFit);
            float deltaChiSquaredUpDownPerHit(directionFit.GetUpDownDeltaChiSquaredPerHit());

            bool pfoIntersectsTopFace(this->IntersectsYFace(directionFit)), pfoFiducialLowY(this->HasFiducialLowY(directionFit)), pfoHighTopY(this->HasHighTopY(directionFit, 100.0));
            bool isTargetTopFacePfo(pfoHighTopY && pfoIntersectsTopFace && pfoFiducialLowY);

            nuRecoTrueCosmicRay = isTrueCosmicRay;
            nuRecoIntersectsTopFace = pfoIntersectsTopFace;
            nuRecoFiducialLowY = pfoFiducialLowY;
            nuRecoHighTopY = pfoHighTopY;
            nuRecoPolarAngle = pfoPolarAngle;

            nuRecoDeltaChiSquared = deltaChiSquaredUpDownPerHit;
            nuRecoCosmicProbability = pfoCosmicProbability;

            nuRecoIsTargetPfo = isTargetTopFacePfo;
            nuRecoIsTrueTargetCR = isTargetCosmic;

            std::cout << "No exception 2" << std::endl;
        }
        catch (...)
        {
            std::cout << "Exception 2" << std::endl;
            //continue;
        }
    }

    ////////

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputMCParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputHitListName, pCaloHitList));

    LArMCParticleHelper::MCRelationMap mcPrimaryMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);

    LArMCParticleHelper::CaloHitToMCMap hitToMCMap;
    LArMCParticleHelper::MCContributionMap mcToTrueHitListMap; 
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcPrimaryMap, hitToMCMap, mcToTrueHitListMap);

    //Front PFO in selectedSlicePfos is NuReco (assuming max number neutrinos is 1)
    pandora::CaloHitList nuRecoPfoHits;
    LArPfoHelper::GetCaloHits(selectedSlicePfos.front(), TPC_3D, nuRecoPfoHits);

    int nNeutrinoInducedHits(1), nTotalHits(1);

    for (const auto pCaloHit : nuRecoPfoHits)
    {
        const pandora::CaloHit* pParentCaloHit(static_cast<const CaloHit*>(pCaloHit->GetParentAddress()));

        ++nTotalHits;

        if (hitToMCMap.find(pParentCaloHit) != hitToMCMap.end())
        {
            if (LArMCParticleHelper::IsBeamNeutrinoFinalState(hitToMCMap.at(pParentCaloHit)))
                ++nNeutrinoInducedHits;
        }
    }

    float neutrinoHitPurityFraction(static_cast<float>(nNeutrinoInducedHits)/nTotalHits);

    int interactionTypeInt(-1);

    try
    {
        LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, LArMCParticleHelper::PrimaryParameters(),
            LArMCParticleHelper::IsBeamNeutrinoFinalState, nuMCParticlesToGoodHitsMap);

        MCParticleList mcPrimaryList;
        for (const auto &mapEntry : nuMCParticlesToGoodHitsMap) mcPrimaryList.push_back(mapEntry.first);
        mcPrimaryList.sort(LArMCParticleHelper::SortByMomentum);

        const LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(mcPrimaryList));
        interactionTypeInt = interactionType;
    }
    catch (...)
    {
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "EventContainsTargetPfo", eventContainsTargetPfo ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "EventContainsTrueTargetCR", eventContainsTrueTargetCR ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "NuRecoIsTargetPfo", nuRecoIsTargetPfo ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "NuRecoIsTrueTargetCR", nuRecoIsTrueTargetCR? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "NuRecoTrueCosmicRay", nuRecoTrueCosmicRay ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "NuRecoIntersectsTopFace", nuRecoIntersectsTopFace ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "NuRecoFiducialLowY", nuRecoFiducialLowY ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "NuRecoHighTopY", nuRecoHighTopY ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "NuRecoPolarAngle", nuRecoPolarAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "NuRecoHitPurityFraction", neutrinoHitPurityFraction));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "NuRecoDeltaChiSquared", nuRecoDeltaChiSquared));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "NuRecoCosmicProbability", nuRecoCosmicProbability));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "NeutrinoIdCorrect", neutrinoHitPurityFraction > 0.9 ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "TrueInteractionType", interactionTypeInt));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "FileIdentifier", m_fileIdentifier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "NeutrinoIdCheck", "EventNumber", m_eventNumber));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "NeutrinoIdCheck"));
    */

    PfoList newSlicePfoList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Recreate(selectedSlicePfos, newSlicePfoList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Reset()
{
    for (const Pandora *const pCRWorker : m_crWorkerInstances)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pCRWorker));

    if (m_pSlicingWorkerInstance)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSlicingWorkerInstance));

    if (m_pSliceNuWorkerInstance)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSliceNuWorkerInstance));

    if (m_pSliceCRWorkerInstance)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSliceCRWorkerInstance));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Copy(const Pandora *const pPandora, const CaloHit *const pCaloHit) const
{
    // TODO Protection to ensure that only calo hits owned by the master instance are copied?
    // TODO Might ultimately want to create LArCaloHits in worker instances, rather than CaloHits, just as we create LArMCParticles, below
    PandoraApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = pCaloHit->GetPositionVector();
    parameters.m_expectedDirection = pCaloHit->GetExpectedDirection();
    parameters.m_cellNormalVector = pCaloHit->GetCellNormalVector();
    parameters.m_cellGeometry = pCaloHit->GetCellGeometry();
    parameters.m_cellSize0 = pCaloHit->GetCellSize0();
    parameters.m_cellSize1 = pCaloHit->GetCellSize1();
    parameters.m_cellThickness = pCaloHit->GetCellThickness();
    parameters.m_nCellRadiationLengths = pCaloHit->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pCaloHit->GetNCellInteractionLengths();
    parameters.m_time = pCaloHit->GetTime();
    parameters.m_inputEnergy = pCaloHit->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pCaloHit->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pCaloHit->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pCaloHit->GetHadronicEnergy();
    parameters.m_isDigital = pCaloHit->IsDigital();
    parameters.m_hitType = pCaloHit->GetHitType();
    parameters.m_hitRegion = pCaloHit->GetHitRegion();
    parameters.m_layer = pCaloHit->GetLayer();
    parameters.m_isInOuterSamplingLayer = pCaloHit->IsInOuterSamplingLayer();
    // ATTN Parent of calo hit in worker is corresponding calo hit in master
    parameters.m_pParentAddress = static_cast<const void*>(pCaloHit);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPandora, parameters));

    if (m_passMCParticlesToWorkerInstances)
    {
        MCParticleVector mcParticleVector;
        for (const auto &weightMapEntry : pCaloHit->GetMCParticleWeightMap()) mcParticleVector.push_back(weightMapEntry.first);
        std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(*pPandora, pCaloHit,
                pMCParticle, pCaloHit->GetMCParticleWeightMap().at(pMCParticle)));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Copy(const Pandora *const pPandora, const MCParticle *const pMCParticle, const LArMCParticleFactory *const pMCParticleFactory) const
{
    LArMCParticleParameters parameters;
    const LArMCParticle *const pLArMCParticle = dynamic_cast<const LArMCParticle*>(pMCParticle);

    if (!pLArMCParticle)
    {
        std::cout << "MasterAlgorithm::Copy - Expect to pass only LArMCParticles to Pandora worker instances." << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    parameters.m_nuanceCode = pLArMCParticle->GetNuanceCode();
    parameters.m_energy = pMCParticle->GetEnergy();
    parameters.m_momentum = pMCParticle->GetMomentum();
    parameters.m_vertex = pMCParticle->GetVertex();
    parameters.m_endpoint = pMCParticle->GetEndpoint();
    parameters.m_particleId = pMCParticle->GetParticleId();
    parameters.m_mcParticleType = pMCParticle->GetMCParticleType();
    // ATTN Parent of mc particle in worker is corresponding mc particle in master
    parameters.m_pParentAddress = static_cast<const void*>(pMCParticle);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPandora, parameters, *pMCParticleFactory));

    for (const MCParticle *const pDaughterMCParticle : pMCParticle->GetDaughterList())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*pPandora, pMCParticle, pDaughterMCParticle));

    for (const MCParticle *const pParentMCParticle : pMCParticle->GetParentList())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*pPandora, pParentMCParticle, pMCParticle));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Recreate(const PfoList &inputPfoList, PfoList &newPfoList) const
{
    if (inputPfoList.empty())
        return STATUS_CODE_SUCCESS;

    // TODO if no pfo in input list is primary - raise exception

    std::string clusterListName;
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, clusterListName));

    std::string vertexListName;
    const VertexList *pVertexList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    std::string pfoListName;
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (const Pfo *const pPfo : inputPfoList)
    {
        if (pPfo->GetParentPfoList().empty())
            this->Recreate(pPfo, nullptr, newPfoList);
    }

    if (!pClusterList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_recreatedClusterListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_recreatedClusterListName));
    }

    if (!pVertexList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_recreatedVertexListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_recreatedVertexListName));
    }

    if (!pPfoList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<ParticleFlowObject>(*this, m_recreatedPfoListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<ParticleFlowObject>(*this, m_recreatedPfoListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Recreate(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject *const pNewParentPfo, PfoList &newPfoList) const
{
    ClusterList inputClusterList2D, inputClusterList3D, newClusterList;
    LArPfoHelper::GetTwoDClusterList(pInputPfo, inputClusterList2D);
    LArPfoHelper::GetThreeDClusterList(pInputPfo, inputClusterList3D);

    for (const Cluster *const pInputCluster : inputClusterList2D)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
            newCaloHitList.push_back(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
            newIsolatedCaloHitList.push_back(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));

        if (!newCaloHitList.empty())
            newClusterList.push_back(this->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }

    for (const Cluster *const pInputCluster : inputClusterList3D)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
        {
            const CaloHit *const pWorkerParentCaloHit(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));
            const CaloHit *const pMasterParentCaloHit(static_cast<const CaloHit*>(pWorkerParentCaloHit->GetParentAddress()));
            newCaloHitList.push_back(this->CreateCaloHit(pInputCaloHit, pMasterParentCaloHit));
        }

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
        {
            const CaloHit *const pWorkerParentCaloHit(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));
            const CaloHit *const pMasterParentCaloHit(static_cast<const CaloHit*>(pWorkerParentCaloHit->GetParentAddress()));
            newIsolatedCaloHitList.push_back(this->CreateCaloHit(pInputCaloHit, pMasterParentCaloHit));
        }

        if (!newCaloHitList.empty())
            newClusterList.push_back(this->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }

    VertexList newVertexList;

    for (const Vertex *const pInputVertex : pInputPfo->GetVertexList())
        newVertexList.push_back(this->CreateVertex(pInputVertex));

    const ParticleFlowObject *const pNewPfo(this->CreatePfo(pInputPfo, newClusterList, newVertexList));
    newPfoList.push_back(pNewPfo);

    if (pNewParentPfo)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNewParentPfo, pNewPfo))

    for (const ParticleFlowObject *const pInputDaughterPfo : pInputPfo->GetDaughterPfoList())
        this->Recreate(pInputDaughterPfo, pNewPfo, newPfoList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *MasterAlgorithm::CreateCaloHit(const CaloHit *const pInputCaloHit, const CaloHit *const pParentCaloHit) const
{
    PandoraContentApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = pInputCaloHit->GetPositionVector();
    parameters.m_expectedDirection = pInputCaloHit->GetExpectedDirection();
    parameters.m_cellNormalVector = pInputCaloHit->GetCellNormalVector();
    parameters.m_cellGeometry = pInputCaloHit->GetCellGeometry();
    parameters.m_cellSize0 = pInputCaloHit->GetCellSize0();
    parameters.m_cellSize1 = pInputCaloHit->GetCellSize1();
    parameters.m_cellThickness = pInputCaloHit->GetCellThickness();
    parameters.m_nCellRadiationLengths = pInputCaloHit->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pInputCaloHit->GetNCellInteractionLengths();
    parameters.m_time = pInputCaloHit->GetTime();
    parameters.m_inputEnergy = pInputCaloHit->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pInputCaloHit->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pInputCaloHit->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pInputCaloHit->GetHadronicEnergy();
    parameters.m_isDigital = pInputCaloHit->IsDigital();
    parameters.m_hitType = pInputCaloHit->GetHitType();
    parameters.m_hitRegion = pInputCaloHit->GetHitRegion();
    parameters.m_layer = pInputCaloHit->GetLayer();
    parameters.m_isInOuterSamplingLayer = pInputCaloHit->IsInOuterSamplingLayer();
    parameters.m_pParentAddress = static_cast<const void*>(pParentCaloHit);

    const CaloHit *pNewCaloHit(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pNewCaloHit));

    PandoraContentApi::CaloHit::Metadata metadata;
    metadata.m_isIsolated = pInputCaloHit->IsIsolated();
    metadata.m_isPossibleMip = pInputCaloHit->IsPossibleMip();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::AlterMetadata(*this, pNewCaloHit, metadata));

    return pNewCaloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *MasterAlgorithm::CreateCluster(const Cluster *const pInputCluster, const CaloHitList &newCaloHitList,
    const CaloHitList &newIsolatedCaloHitList) const
{
    PandoraContentApi::Cluster::Parameters parameters;
    parameters.m_caloHitList = newCaloHitList;
    parameters.m_isolatedCaloHitList = newIsolatedCaloHitList;

    const Cluster *pNewCluster(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));

    PandoraContentApi::Cluster::Metadata metadata;
    metadata.m_particleId = pInputCluster->GetParticleId();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pNewCluster, metadata));

    return pNewCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Vertex *MasterAlgorithm::CreateVertex(const Vertex *const pInputVertex) const
{
    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = pInputVertex->GetPosition();
    parameters.m_vertexLabel = pInputVertex->GetVertexLabel();
    parameters.m_vertexType = pInputVertex->GetVertexType();

    const Vertex *pNewVertex(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pNewVertex));

    return pNewVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *MasterAlgorithm::CreatePfo(const ParticleFlowObject *const pInputPfo, const ClusterList &newClusterList,
    const VertexList &newVertexList) const
{
    PandoraContentApi::ParticleFlowObject::Parameters parameters;
    parameters.m_particleId = pInputPfo->GetParticleId();
    parameters.m_charge = pInputPfo->GetCharge();
    parameters.m_mass = pInputPfo->GetMass();
    parameters.m_energy = pInputPfo->GetEnergy();
    parameters.m_momentum = pInputPfo->GetMomentum();
    parameters.m_clusterList = newClusterList;
    parameters.m_trackList.clear();
    parameters.m_vertexList = newVertexList;
    parameters.m_propertiesToAdd = pInputPfo->GetPropertiesMap();

    const ParticleFlowObject *pNewPfo(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, parameters, pNewPfo));

    return pNewPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Pandora *MasterAlgorithm::CreateWorkerInstance(const LArTPC &larTPC, const DetectorGapList &gapList, const std::string &settingsFile) const
{
    // The Pandora instance
    const Pandora *const pPandora(new Pandora());
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPandora, new lar_content::LArRotationalTransformationPlugin));
    MultiPandoraApi::AddDaughterPandoraInstance(&(this->GetPandora()), pPandora);

    // The LArTPC
    PandoraApi::Geometry::LArTPC::Parameters larTPCParameters;
    larTPCParameters.m_larTPCVolumeId = larTPC.GetLArTPCVolumeId();
    larTPCParameters.m_centerX = larTPC.GetCenterX();
    larTPCParameters.m_centerY = larTPC.GetCenterY();
    larTPCParameters.m_centerZ = larTPC.GetCenterZ();
    larTPCParameters.m_widthX = larTPC.GetWidthX();
    larTPCParameters.m_widthY = larTPC.GetWidthY();
    larTPCParameters.m_widthZ = larTPC.GetWidthZ();
    larTPCParameters.m_wirePitchU = larTPC.GetWirePitchU();
    larTPCParameters.m_wirePitchV = larTPC.GetWirePitchV();
    larTPCParameters.m_wirePitchW = larTPC.GetWirePitchW();
    larTPCParameters.m_wireAngleU = larTPC.GetWireAngleU();
    larTPCParameters.m_wireAngleV = larTPC.GetWireAngleV();
    larTPCParameters.m_wireAngleW = larTPC.GetWireAngleW();
    larTPCParameters.m_sigmaUVW = larTPC.GetSigmaUVW();
    larTPCParameters.m_isDriftInPositiveX = larTPC.IsDriftInPositiveX();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPandora, larTPCParameters));

    const float tpcMinX(larTPC.GetCenterX() - 0.5f * larTPC.GetWidthX()), tpcMaxX(larTPC.GetCenterX() + 0.5f * larTPC.GetWidthX());

    // The Gaps
    for (const DetectorGap *const pGap : gapList)
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap*>(pGap));

        if (pLineGap && (((pLineGap->GetLineEndX() >= tpcMinX) && (pLineGap->GetLineEndX() <= tpcMaxX)) ||
            ((pLineGap->GetLineStartX() >= tpcMinX) && (pLineGap->GetLineStartX() <= tpcMaxX))) )
        {
            PandoraApi::Geometry::LineGap::Parameters lineGapParameters;
            const LineGapType lineGapType(pLineGap->GetLineGapType());
            lineGapParameters.m_lineGapType = lineGapType;
            lineGapParameters.m_lineStartX = pLineGap->GetLineStartX();
            lineGapParameters.m_lineEndX = pLineGap->GetLineEndX();

            if (m_fullWidthCRWorkerWireGaps && ((lineGapType == TPC_WIRE_GAP_VIEW_U) || (lineGapType == TPC_WIRE_GAP_VIEW_V) || (lineGapType == TPC_WIRE_GAP_VIEW_W)))
            {
                lineGapParameters.m_lineStartX = -std::numeric_limits<float>::max();
                lineGapParameters.m_lineEndX = std::numeric_limits<float>::max();
            }

            lineGapParameters.m_lineStartZ = pLineGap->GetLineStartZ();
            lineGapParameters.m_lineEndZ = pLineGap->GetLineEndZ();
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LineGap::Create(*pPandora, lineGapParameters));
        }
    }

    // Configuration
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, settingsFile));
    return pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Pandora *MasterAlgorithm::CreateWorkerInstance(const LArTPCMap &larTPCMap, const DetectorGapList &gapList, const std::string &settingsFile) const
{
    if (larTPCMap.empty())
    {
        std::cout << "MasterAlgorithm::CreateWorkerInstance - no LArTPC details provided" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    // The Pandora instance
    const Pandora *const pPandora(new Pandora());
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPandora, new lar_content::LArRotationalTransformationPlugin));
    MultiPandoraApi::AddDaughterPandoraInstance(&(this->GetPandora()), pPandora);

    // The Parent LArTPC
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);
    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    PandoraApi::Geometry::LArTPC::Parameters larTPCParameters;
    larTPCParameters.m_larTPCVolumeId = 0;
    larTPCParameters.m_centerX = 0.5f * (parentMaxX + parentMinX);
    larTPCParameters.m_centerY = 0.5f * (parentMaxY + parentMinY);
    larTPCParameters.m_centerZ = 0.5f * (parentMaxZ + parentMinZ);
    larTPCParameters.m_widthX = parentMaxX - parentMinX;
    larTPCParameters.m_widthY = parentMaxY - parentMinY;
    larTPCParameters.m_widthZ = parentMaxZ - parentMinZ;
    larTPCParameters.m_wirePitchU = std::max(pFirstLArTPC->GetWirePitchU(), pFirstLArTPC->GetWirePitchV());
    larTPCParameters.m_wirePitchV = std::max(pFirstLArTPC->GetWirePitchU(), pFirstLArTPC->GetWirePitchV());
    larTPCParameters.m_wirePitchW = pFirstLArTPC->GetWirePitchW();
    larTPCParameters.m_wireAngleU = pFirstLArTPC->GetWireAngleU();
    larTPCParameters.m_wireAngleV = pFirstLArTPC->GetWireAngleV();
    larTPCParameters.m_wireAngleW = pFirstLArTPC->GetWireAngleW();
    larTPCParameters.m_sigmaUVW = pFirstLArTPC->GetSigmaUVW();
    larTPCParameters.m_isDriftInPositiveX = pFirstLArTPC->IsDriftInPositiveX();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPandora, larTPCParameters));

    // The Gaps
    for (const DetectorGap *const pGap : gapList)
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap*>(pGap));

        if (pLineGap)
        {
            PandoraApi::Geometry::LineGap::Parameters lineGapParameters;
            lineGapParameters.m_lineGapType = pLineGap->GetLineGapType();
            lineGapParameters.m_lineStartX = pLineGap->GetLineStartX();
            lineGapParameters.m_lineEndX = pLineGap->GetLineEndX();
            lineGapParameters.m_lineStartZ = pLineGap->GetLineStartZ();
            lineGapParameters.m_lineEndZ = pLineGap->GetLineEndZ();
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LineGap::Create(*pPandora, lineGapParameters));
        }
    }

    // Configuration
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, settingsFile));
    return pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    LArSpaceChargeHelper::Configure("/usera/jjd49/pandora_direction/PandoraPFA/LArContent-origin/vertex_direction/larpandoracontent/LArDirection/SCEoffsets_MicroBooNE_E273.root");

    ExternalSteeringParameters *pExternalParameters(nullptr);

    if (this->ExternalParametersPresent())
    {
        pExternalParameters = dynamic_cast<ExternalSteeringParameters*>(this->GetExternalParameters());
        if (!pExternalParameters) return STATUS_CODE_FAILURE;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunAllHitsCosmicReco, xmlHandle, "ShouldRunAllHitsCosmicReco", m_shouldRunAllHitsCosmicReco));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunStitching, xmlHandle, "ShouldRunStitching", m_shouldRunStitching));

    if (m_shouldRunStitching && !m_shouldRunAllHitsCosmicReco)
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldRunStitching requires ShouldRunAllHitsCosmicReco to be true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldRunStitching)
    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "StitchingTools", algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            StitchingBaseTool *const pStitchingTool(dynamic_cast<StitchingBaseTool*>(pAlgorithmTool));
            if (!pStitchingTool) return STATUS_CODE_INVALID_PARAMETER;
            m_stitchingToolVector.push_back(pStitchingTool);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunCosmicHitRemoval, xmlHandle, "ShouldRunCosmicHitRemoval", m_shouldRunCosmicHitRemoval));

    if (m_shouldRunCosmicHitRemoval && !m_shouldRunAllHitsCosmicReco)
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldRunCosmicHitRemoval requires ShouldRunAllHitsCosmicReco to be true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldRunCosmicHitRemoval)
    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "CosmicRayTaggingTools", algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            CosmicRayTaggingBaseTool *const pCosmicRayTaggingTool(dynamic_cast<CosmicRayTaggingBaseTool*>(pAlgorithmTool));
            if (!pCosmicRayTaggingTool) return STATUS_CODE_INVALID_PARAMETER;
            m_cosmicRayTaggingToolVector.push_back(pCosmicRayTaggingTool);
        }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------------

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
    "CheatTopFaceCosmicRemoval", m_cheatTopFaceCosmicRemoval));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
    "RecoTopFaceCosmicRemoval", m_recoTopFaceCosmicRemoval));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
    "WriteToTree", m_writeToTree));

    if (m_writeToTree)
    {    
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
            "FileIdentifier", m_fileIdentifier));
    } 

    //Track direction tool
    AlgorithmTool *pTrackDirectionTool(nullptr);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "TrackDirection", pTrackDirectionTool));

    if (!(this->m_pTrackDirectionTool = dynamic_cast<TrackDirectionTool *>(pTrackDirectionTool)))
        throw STATUS_CODE_FAILURE;

    //-----------------------------------------------------------------------------------------------------------------------------------------------

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunSlicing, xmlHandle, "ShouldRunSlicing", m_shouldRunSlicing));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunNeutrinoRecoOption, xmlHandle, "ShouldRunNeutrinoRecoOption", m_shouldRunNeutrinoRecoOption));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunCosmicRecoOption, xmlHandle, "ShouldRunCosmicRecoOption", m_shouldRunCosmicRecoOption));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldPerformSliceId, xmlHandle, "ShouldPerformSliceId", m_shouldPerformSliceId));

    if (m_shouldPerformSliceId && (!m_shouldRunSlicing || !m_shouldRunNeutrinoRecoOption || !m_shouldRunCosmicRecoOption))
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldPerformSliceId requires ShouldRunSlicing and both neutrino and cosmic reconstruction options" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldPerformSliceId)
    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "SliceIdTools", algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            SliceIdBaseTool *const pSliceIdIdTool(dynamic_cast<SliceIdBaseTool*>(pAlgorithmTool));
            if (!pSliceIdIdTool) return STATUS_CODE_INVALID_PARAMETER;
            m_sliceIdToolVector.push_back(pSliceIdIdTool);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_printOverallRecoStatus, xmlHandle, "PrintOverallRecoStatus", m_printOverallRecoStatus));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualizeOverallRecoStatus", m_visualizeOverallRecoStatus));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FullWidthCRWorkerWireGaps", m_fullWidthCRWorkerWireGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PassMCParticlesToWorkerInstances", m_passMCParticlesToWorkerInstances));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CRSettingsFile", m_crSettingsFile));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NuSettingsFile", m_nuSettingsFile));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SlicingSettingsFile", m_slicingSettingsFile));
    m_crSettingsFile = LArFileHelper::FindFileInPath(m_crSettingsFile, m_filePathEnvironmentVariable);
    m_nuSettingsFile = LArFileHelper::FindFileInPath(m_nuSettingsFile, m_filePathEnvironmentVariable);
    m_slicingSettingsFile = LArFileHelper::FindFileInPath(m_slicingSettingsFile, m_filePathEnvironmentVariable);

    if (m_passMCParticlesToWorkerInstances)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputMCParticleListName", m_inputMCParticleListName));
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputHitListName", m_inputHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedPfoListName", m_recreatedPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedClusterListName", m_recreatedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedVertexListName", m_recreatedVertexListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InTimeMaxX0", m_inTimeMaxX0));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::ReadExternalSettings(const ExternalSteeringParameters *const pExternalParameters, const InputBool inputBool,
    const TiXmlHandle xmlHandle, const std::string &xmlTag, bool &outputBool)
{
    if (pExternalParameters && inputBool.IsInitialized())
    {
        outputBool = inputBool.Get();
    }
    else
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, xmlTag, outputBool));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
