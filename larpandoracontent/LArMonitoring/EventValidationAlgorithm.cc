/**
 *  @file   larpandoracontent/LArMonitoring/EventValidationAlgorithm.cc
 *
 *  @brief  Implementation of the event validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"

#include "larpandoracontent/LArMonitoring/EventValidationAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/SvmPfoCharacterisationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArSpaceChargeHelper.h"

#include <sstream>
#include <math.h> 

using namespace pandora;

namespace lar_content
{

EventValidationAlgorithm::EventValidationAlgorithm() :
    m_useTrueNeutrinosOnly(false),
    m_testBeamMode(false),
    m_selectInputHits(true),
    m_minHitSharingFraction(0.9f),
    m_maxPhotonPropagation(2.5f),
    m_printAllToScreen(false),
    m_printMatchingToScreen(true),
    m_writeToTree(false),
    m_useSmallPrimaries(true),
    m_matchingMinSharedHits(5),
    m_matchingMinCompleteness(0.1f),
    m_matchingMinPurity(0.5f),
    m_slidingFitWindow(20),
    m_viewEvent(false),
    m_fileIdentifier(0),
    m_eventNumber(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventValidationAlgorithm::~EventValidationAlgorithm()
{
    if (m_writeToTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "EventSelection", m_fileName.c_str(), "UPDATE"));
            //PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "HitInformation", m_fileName.c_str(), "UPDATE"));
            //PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "BraggHitInformation", m_fileName.c_str(), "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "EventValidationAlgorithm: Unable to write tree " << m_treeName << " to file " << m_fileName << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationAlgorithm::Run()
{
    ++m_eventNumber;

    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    const PfoList *pPfoList = nullptr;
    (void) PandoraContentApi::GetList(*this, m_pfoListName, pPfoList);

    ValidationInfo validationInfo;
    this->FillValidationInfo(pMCParticleList, pCaloHitList, pPfoList, validationInfo);

    if (m_printAllToScreen)
        this->PrintAllMatches(validationInfo);

    if (m_printMatchingToScreen)
        this->PrintInterpretedMatches(validationInfo);

    if (m_writeToTree)
        this->WriteInterpretedMatches(validationInfo);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::FillValidationInfo(const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList,
    const PfoList *const pPfoList, ValidationInfo &validationInfo) const
{
    if (pMCParticleList && pCaloHitList)
    {
        LArMCParticleHelper::PrimaryParameters parameters;

        parameters.m_selectInputHits = m_selectInputHits;
        parameters.m_minHitSharingFraction = m_minHitSharingFraction;
        parameters.m_maxPhotonPropagation = m_maxPhotonPropagation;
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamParticle, targetMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, targetMCParticleToHitsMap);

        parameters.m_minPrimaryGoodHits = 0;
        parameters.m_minHitsForGoodView = 0;
        parameters.m_minHitSharingFraction = 0.f;
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, allMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamParticle, allMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, allMCParticleToHitsMap);

        validationInfo.SetTargetMCParticleToHitsMap(targetMCParticleToHitsMap);
        validationInfo.SetAllMCParticleToHitsMap(allMCParticleToHitsMap);
    }
    
    if (pPfoList)
    {
        PfoList allConnectedPfos;
        LArPfoHelper::GetAllConnectedPfos(*pPfoList, allConnectedPfos);

        PfoList finalStatePfos;
        for (const ParticleFlowObject *const pPfo : allConnectedPfos)
        {
            if ((!m_testBeamMode && LArPfoHelper::IsFinalState(pPfo)) || (m_testBeamMode && pPfo->GetParentPfoList().empty()))
                finalStatePfos.push_back(pPfo);
        }

        LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, validationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap);
        validationInfo.SetPfoToHitsMap(pfoToHitsMap);
    }

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(validationInfo.GetPfoToHitsMap(), {validationInfo.GetAllMCParticleToHitsMap()}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    validationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);
    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::ProcessOutput(const ValidationInfo &validationInfo, const bool useInterpretedMatching, const bool printToScreen, const bool fillTree) const
{
    if (printToScreen && useInterpretedMatching) std::cout << "---INTERPRETED-MATCHING-OUTPUT------------------------------------------------------------------" << std::endl;
    else if (printToScreen) std::cout << "---RAW-MATCHING-OUTPUT--------------------------------------------------------------------------" << std::endl;

    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap(useInterpretedMatching ?
        validationInfo.GetInterpretedMCToPfoHitSharingMap() : validationInfo.GetMCToPfoHitSharingMap());

    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetTargetMCParticleToHitsMap()}, mcPrimaryVector);

    int nNeutrinoPrimaries(0);
    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) && validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary)) ++nNeutrinoPrimaries;

    PfoVector primaryPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(validationInfo.GetPfoToHitsMap(), primaryPfoVector);

    int pfoIndex(0), neutrinoPfoIndex(0);
    PfoToIdMap pfoToIdMap, neutrinoPfoToIdMap;
    for (const Pfo *const pPrimaryPfo : primaryPfoVector)
    {
        pfoToIdMap.insert(PfoToIdMap::value_type(pPrimaryPfo, ++pfoIndex));
        const Pfo *const pRecoNeutrino(LArPfoHelper::IsNeutrinoFinalState(pPrimaryPfo) ? LArPfoHelper::GetParentNeutrino(pPrimaryPfo) : nullptr);

        if (pRecoNeutrino && !neutrinoPfoToIdMap.count(pRecoNeutrino))
            neutrinoPfoToIdMap.insert(PfoToIdMap::value_type(pRecoNeutrino, ++neutrinoPfoIndex));
    }

    PfoSet recoNeutrinos;
    MCParticleList associatedMCPrimaries;

    int nCorrectNu(0), nTotalNu(0), nCorrectTB(0), nTotalTB(0), nCorrectCR(0), nTotalCR(0), nFakeNu(0), nFakeCR(0), nSplitNu(0), nSplitCR(0), nLost(0);
    int mcPrimaryIndex(0), nTargetMatches(0), nTargetNuMatches(0), nTargetCRMatches(0), nTargetGoodNuMatches(0), nTargetNuSplits(0), nTargetNuLosses(0);
    IntVector mcPrimaryId, mcPrimaryPdg, nMCHitsTotal, nMCHitsU, nMCHitsV, nMCHitsW;
    FloatVector mcPrimaryE, mcPrimaryPX, mcPrimaryPY, mcPrimaryPZ;
    FloatVector mcPrimaryVtxX, mcPrimaryVtxY, mcPrimaryVtxZ, mcPrimaryEndX, mcPrimaryEndY, mcPrimaryEndZ;
    IntVector nPrimaryMatchedPfos, nPrimaryMatchedNuPfos, nPrimaryMatchedCRPfos;
    IntVector bestMatchPfoId, bestMatchPfoPdg, bestMatchPfoIsRecoNu, bestMatchPfoRecoNuId;
    IntVector bestMatchPfoNHitsTotal, bestMatchPfoNHitsU, bestMatchPfoNHitsV, bestMatchPfoNHitsW;
    IntVector bestMatchPfoNSharedHitsTotal, bestMatchPfoNSharedHitsU, bestMatchPfoNSharedHitsV, bestMatchPfoNSharedHitsW;

    //NEW CODE
    IntVector nNuRecoCosmicRayMatches;

    std::stringstream targetSS;

    //NEW CODE
    int nNuParticleMatches(0), nRecoNuGoodParticleMatches(0), nRecoNuParticleMatches(0), nMuonParticleMatches(0), nProtonParticleMatches(0);
    int nNuTrackMatches(0), nRecoNuTrackMatches(0), nMuonTrackMatches(0), nProtonTrackMatches(0);
    int interactionTypeCopy(167), modifiedInteractionTypeCopy(167), nuanceCodeCopy(-1);
    int nRecoNuSplits(-1);
    int nMuons(0), nProtons(0), nNeutrons(0), nPhotons(0), nPiPlus(0), nPiMinus(0), nPiZero(0), nElectrons(0), nOther(0);
    int nPrimaries(0);
    int nRecoNuCosmicRays(0);

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        const bool hasMatch(mcToPfoHitSharingMap.count(pMCPrimary) && !mcToPfoHitSharingMap.at(pMCPrimary).empty());
        const bool isTargetPrimary(validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary)); //is reconstructable

        if (!hasMatch && !isTargetPrimary)
            continue;

        //NEW CODE
        int mcParticlePdg(pMCPrimary->GetParticleId());

        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) && isTargetPrimary)
        {
            ++nPrimaries;
            if (std::abs(mcParticlePdg) == 13) ++nMuons;
            else if (std::abs(mcParticlePdg) == 2212) ++nProtons;
            else if (std::abs(mcParticlePdg) == 11) ++nElectrons;
            else if (mcParticlePdg == 2112) ++nNeutrons;
            else if (mcParticlePdg == 22) ++nPhotons;
            else if (mcParticlePdg == 211) ++nPiPlus;
            else if (mcParticlePdg == -211) ++nPiMinus;
            else if (mcParticlePdg == 111) ++nPiZero;
            else ++nOther;
        }

        associatedMCPrimaries.push_back(pMCPrimary);
        const int nTargetPrimaries(associatedMCPrimaries.size());
        const bool isLastNeutrinoPrimary(++mcPrimaryIndex == nNeutrinoPrimaries);
        const CaloHitList &mcPrimaryHitList(validationInfo.GetAllMCParticleToHitsMap().at(pMCPrimary));

        const int mcNuanceCode(LArMCParticleHelper::GetNuanceCode(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)));

        if (LArMCParticleHelper::IsNeutrino(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)))
            nuanceCodeCopy = mcNuanceCode;

        const int isBeamNeutrinoFinalState(LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary));
        const int isBeamParticle(LArMCParticleHelper::IsBeamParticle(pMCPrimary));
        const int isCosmicRay(LArMCParticleHelper::IsCosmicRay(pMCPrimary));
#ifdef MONITORING
        const CartesianVector &targetVertex(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)->GetVertex());
        const float targetVertexX(targetVertex.GetX()), targetVertexY(targetVertex.GetY()), targetVertexZ(targetVertex.GetZ());
#endif
        targetSS << (!isTargetPrimary ? "(Non target) " : "")
                 << "PrimaryId " << mcPrimaryIndex
                 << ", Nu " << isBeamNeutrinoFinalState
                 << ", TB " << isBeamParticle
                 << ", CR " << isCosmicRay
                 << ", MCPDG " << pMCPrimary->GetParticleId()
                 << ", Energy " << pMCPrimary->GetEnergy()
                 << ", Dist. " << (pMCPrimary->GetEndpoint() - pMCPrimary->GetVertex()).GetMagnitude()
                 << ", nMCHits " << mcPrimaryHitList.size()
                 << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, mcPrimaryHitList)
                 << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, mcPrimaryHitList)
                 << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, mcPrimaryHitList) << ")" << std::endl;

        mcPrimaryId.push_back(mcPrimaryIndex);
        mcPrimaryPdg.push_back(pMCPrimary->GetParticleId());
        mcPrimaryE.push_back(pMCPrimary->GetEnergy());
        mcPrimaryPX.push_back(pMCPrimary->GetMomentum().GetX());
        mcPrimaryPY.push_back(pMCPrimary->GetMomentum().GetY());
        mcPrimaryPZ.push_back(pMCPrimary->GetMomentum().GetZ());
        mcPrimaryVtxX.push_back(pMCPrimary->GetVertex().GetX());
        mcPrimaryVtxY.push_back(pMCPrimary->GetVertex().GetY());
        mcPrimaryVtxZ.push_back(pMCPrimary->GetVertex().GetZ());
        mcPrimaryEndX.push_back(pMCPrimary->GetEndpoint().GetX());
        mcPrimaryEndY.push_back(pMCPrimary->GetEndpoint().GetY());
        mcPrimaryEndZ.push_back(pMCPrimary->GetEndpoint().GetZ());
        nMCHitsTotal.push_back(mcPrimaryHitList.size());
        nMCHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, mcPrimaryHitList));
        nMCHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, mcPrimaryHitList));
        nMCHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, mcPrimaryHitList));

        int matchIndex(0), nPrimaryMatches(0), nPrimaryNuMatches(0), nPrimaryCRMatches(0), nPrimaryGoodNuMatches(0), nPrimaryNuSplits(0);
        //NEW CODE
        int nPrimaryTrackMatches(0), nPrimaryNuTrackMatches(0);
        int nPrimaryCustomCRMatches(0);
#ifdef MONITORING
        float recoVertexX(std::numeric_limits<float>::max()), recoVertexY(std::numeric_limits<float>::max()), recoVertexZ(std::numeric_limits<float>::max());
#endif
        //Loop over PFOs matched to the MC primary
        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pMCPrimary))
        {
            const CaloHitList &sharedHitList(pfoToSharedHits.second);
            const CaloHitList &pfoHitList(validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first));

            const bool isRecoNeutrinoFinalState(LArPfoHelper::IsNeutrinoFinalState(pfoToSharedHits.first));
            const bool isGoodMatch(this->IsGoodMatch(mcPrimaryHitList, pfoHitList, sharedHitList));

            const int pfoId(pfoToIdMap.at(pfoToSharedHits.first));
            const int recoNuId(isRecoNeutrinoFinalState ? neutrinoPfoToIdMap.at(LArPfoHelper::GetParentNeutrino(pfoToSharedHits.first)) : -1);

            if (0 == matchIndex++)
            {
                bestMatchPfoId.push_back(pfoId);
                bestMatchPfoPdg.push_back(pfoToSharedHits.first->GetParticleId());
                bestMatchPfoIsRecoNu.push_back(isRecoNeutrinoFinalState ? 1 : 0);
                bestMatchPfoRecoNuId.push_back(recoNuId);
                bestMatchPfoNHitsTotal.push_back(pfoHitList.size());
                bestMatchPfoNHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList));
                bestMatchPfoNHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList));
                bestMatchPfoNHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList));
                bestMatchPfoNSharedHitsTotal.push_back(sharedHitList.size());
                bestMatchPfoNSharedHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList));
                bestMatchPfoNSharedHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList));
                bestMatchPfoNSharedHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList));
#ifdef MONITORING
                try
                {
                    const Vertex *const pRecoVertex(LArPfoHelper::GetVertex(isRecoNeutrinoFinalState ? LArPfoHelper::GetParentNeutrino(pfoToSharedHits.first) : pfoToSharedHits.first));
                    pandora::CartesianVector recoVertexPosition(pRecoVertex->GetPosition());
                    const CartesianVector correctedVertexPosition(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(recoVertexPosition));

                    recoVertexX = correctedVertexPosition.GetX();
                    recoVertexY = correctedVertexPosition.GetY();
                    recoVertexZ = correctedVertexPosition.GetZ();
                }
                catch (const StatusCodeException &) {}
#endif
            }

            if (isRecoNeutrinoFinalState && isGoodMatch) ++nPrimaryMatches;
            if (isRecoNeutrinoFinalState && isGoodMatch && LArPfoHelper::IsTrack(pfoToSharedHits.first)) ++nPrimaryTrackMatches;

            if (isRecoNeutrinoFinalState)
            {
                const Pfo *const pRecoNeutrino(LArPfoHelper::GetParentNeutrino(pfoToSharedHits.first));
                const bool isSplitRecoNeutrino(!recoNeutrinos.empty() && !recoNeutrinos.count(pRecoNeutrino));
                if (!isSplitRecoNeutrino && isGoodMatch) ++nPrimaryGoodNuMatches;
                if (isSplitRecoNeutrino && isBeamNeutrinoFinalState && isGoodMatch) ++nPrimaryNuSplits;
                recoNeutrinos.insert(pRecoNeutrino);
            }

            if (!m_testBeamMode)
            {
                if (isRecoNeutrinoFinalState && isGoodMatch) ++nPrimaryNuMatches;
                if (isRecoNeutrinoFinalState && LArPfoHelper::IsTrack(pfoToSharedHits.first) && isGoodMatch) ++nPrimaryNuTrackMatches;
                if (isRecoNeutrinoFinalState && isGoodMatch && LArMCParticleHelper::IsCosmicRay(LArMCParticleHelper::GetParentMCParticle(this->GetMainMCParticle(pfoToSharedHits.first))))
                    ++nPrimaryCustomCRMatches;
                if (!isRecoNeutrinoFinalState && isGoodMatch) ++nPrimaryCRMatches;
            }
            else
            {
                bool isTestBeam(LArPfoHelper::IsTestBeam(pfoToSharedHits.first));
                if (isTestBeam && isGoodMatch) ++nPrimaryNuMatches;
                if (!isTestBeam && isGoodMatch) ++nPrimaryCRMatches;
            }

            targetSS << "-" << (!isGoodMatch ? "(Below threshold) " : "")
                     << "MatchedPfoId " << pfoId
                     << ", Nu " << isRecoNeutrinoFinalState;
            if (isRecoNeutrinoFinalState) targetSS << " [NuId: " << recoNuId << "]";
            targetSS << ", CR " << !isRecoNeutrinoFinalState
                     << ", PDG " << pfoToSharedHits.first->GetParticleId()
                     << ", nMatchedHits " << sharedHitList.size()
                     << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList)
                     << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList)
                     << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList) << ")"
                     << ", nPfoHits " << pfoHitList.size()
                     << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList)
                     << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList)
                     << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList) << ")" << std::endl;
        }

        if (std::abs(pMCPrimary->GetParticleId()) == 13 && LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) && isTargetPrimary)
        {
            nMuonParticleMatches = nPrimaryMatches;
            nMuonTrackMatches = nPrimaryTrackMatches;
        }

        if (std::abs(pMCPrimary->GetParticleId()) == 2212 && LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) && isTargetPrimary)
        {
            nProtonParticleMatches = nPrimaryMatches;
            nProtonTrackMatches = nPrimaryTrackMatches;
        }

        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) && isTargetPrimary)
        {
            nNuParticleMatches += nPrimaryMatches;
            nNuTrackMatches += nPrimaryTrackMatches;
        }

        if (mcToPfoHitSharingMap.at(pMCPrimary).empty())
        {
            targetSS << "-No matched Pfo" << std::endl;
            bestMatchPfoId.push_back(-1); bestMatchPfoPdg.push_back(0); bestMatchPfoIsRecoNu.push_back(0); bestMatchPfoRecoNuId.push_back(-1);
            bestMatchPfoNHitsTotal.push_back(0); bestMatchPfoNHitsU.push_back(0); bestMatchPfoNHitsV.push_back(0); bestMatchPfoNHitsW.push_back(0);
            bestMatchPfoNSharedHitsTotal.push_back(0); bestMatchPfoNSharedHitsU.push_back(0); bestMatchPfoNSharedHitsV.push_back(0); bestMatchPfoNSharedHitsW.push_back(0);
        }

        nPrimaryMatchedPfos.push_back(nPrimaryMatches);
        nPrimaryMatchedNuPfos.push_back(nPrimaryNuMatches);
        nPrimaryMatchedCRPfos.push_back(nPrimaryCRMatches);
        nTargetMatches += nPrimaryMatches;
        nTargetNuMatches += nPrimaryNuMatches;
        nTargetCRMatches += nPrimaryCRMatches;
        nTargetGoodNuMatches += nPrimaryGoodNuMatches;
        nTargetNuSplits += nPrimaryNuSplits;
        if (0 == nPrimaryMatches) ++nTargetNuLosses;

        //NEW CODE
        nRecoNuParticleMatches += nPrimaryNuMatches;
        nRecoNuCosmicRays += nPrimaryCustomCRMatches;

        nRecoNuGoodParticleMatches += nPrimaryGoodNuMatches; //nu reco not split
        nRecoNuTrackMatches += nPrimaryTrackMatches;
        nRecoNuSplits += nPrimaryNuSplits;

        if (fillTree)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "fileIdentifier", m_fileIdentifier));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber - 1));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNuanceCode", mcNuanceCode));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isNeutrino", isBeamNeutrinoFinalState));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isBeamParticle", isBeamParticle));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCosmicRay", isCosmicRay));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetPrimaries", nTargetPrimaries));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "targetVertexX", targetVertexX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "targetVertexY", targetVertexY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "targetVertexZ", targetVertexZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexX", recoVertexX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexY", recoVertexY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexZ", recoVertexZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryId", &mcPrimaryId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPdg", &mcPrimaryPdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryE", &mcPrimaryE));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPX", &mcPrimaryPX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPY", &mcPrimaryPY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPZ", &mcPrimaryPZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxX", &mcPrimaryVtxX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxY", &mcPrimaryVtxY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxZ", &mcPrimaryVtxZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndX", &mcPrimaryEndX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndY", &mcPrimaryEndY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndZ", &mcPrimaryEndZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsTotal", &nMCHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsU", &nMCHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsV", &nMCHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsW", &nMCHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nNuRecoCosmicRayMatches", &nNuRecoCosmicRayMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedPfos", &nPrimaryMatchedPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedNuPfos", &nPrimaryMatchedNuPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedCRPfos", &nPrimaryMatchedCRPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoId", &bestMatchPfoId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoPdg", &bestMatchPfoPdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoIsRecoNu", &bestMatchPfoIsRecoNu));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoRecoNuId", &bestMatchPfoRecoNuId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsTotal", &bestMatchPfoNHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsU", &bestMatchPfoNHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsV", &bestMatchPfoNHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsW", &bestMatchPfoNHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsTotal", &bestMatchPfoNSharedHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsU", &bestMatchPfoNSharedHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsV", &bestMatchPfoNSharedHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsW", &bestMatchPfoNSharedHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetMatches", nTargetMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetNuMatches", nTargetNuMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetCRMatches", nTargetCRMatches));

            if (!m_testBeamMode)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetGoodNuMatches", nTargetGoodNuMatches));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetNuSplits", nTargetNuSplits));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetNuLosses", nTargetNuLosses));
            }
        }

        if (isLastNeutrinoPrimary || isBeamParticle || isCosmicRay)
        {
            const LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(associatedMCPrimaries));
#ifdef MONITORING
            const int interactionTypeInt(static_cast<int>(interactionType));
            int modifiedInteractionTypeInt(static_cast<int>(interactionType));

            //NEW CODE
            //If no neutrino primaries present, default value for modifiedInteractionTypeCopy is 167: NO_NEUTRINO
            //isLastNeutrinoPrimary is false when there are no neutrino primaries
            //isLastNeutrinoPrimary is when all MC primaries are in associatedMCPrimaries, which determines interactionType, but before any CRs (MC primary list is sorted)
            if (isLastNeutrinoPrimary) 
            {
                //These custom interaction types serve to describe what is inside the *reconstructed* neutrino for the purpose of making interaction type tables for a certain particle multiplicity
                if (nRecoNuCosmicRays == nRecoNuParticleMatches)
                    modifiedInteractionTypeInt = 165; //custom CR interaction type

                const MCParticleList *pMCParticleList(nullptr);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

                const CaloHitList *pCaloHitList(nullptr);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

                LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap;
                LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, LArMCParticleHelper::PrimaryParameters(),
                LArMCParticleHelper::IsBeamNeutrinoFinalState, nuMCParticlesToGoodHitsMap);

                MCParticleList mcPrimaryList;
                for (const auto &mapEntry : nuMCParticlesToGoodHitsMap) mcPrimaryList.push_back(mapEntry.first);

                if (mcPrimaryList.size() == 0 && nRecoNuCosmicRays < nRecoNuParticleMatches)
                    modifiedInteractionTypeInt = 166; //often there will be no reconstructable particles left in nu reco due to mistagging, but if a non-CR particle is selected and there are no reconstructable MC particles in the event this constitutes a different class of event: this is a custom interaction type for this case

                interactionTypeCopy = interactionTypeInt;
                modifiedInteractionTypeCopy = modifiedInteractionTypeInt;
            }
#endif
            // ATTN Some redundancy introduced to contributing variables
            const int isCorrectNu(isBeamNeutrinoFinalState && (nTargetGoodNuMatches == nTargetNuMatches) && (nTargetGoodNuMatches == nTargetPrimaries) && (nTargetCRMatches == 0) && (nTargetNuSplits == 0) && (nTargetNuLosses == 0));
            const int isCorrectTB(isBeamParticle && (nTargetNuMatches == 1) && (nTargetCRMatches == 0));
            const int isCorrectCR(isCosmicRay && (nTargetNuMatches == 0) && (nTargetCRMatches == 1));
            const int isFakeNu(isCosmicRay && (nTargetNuMatches > 0));
            const int isFakeCR(!isCosmicRay && (nTargetCRMatches > 0));
            const int isSplitNu(!isCosmicRay && ((nTargetNuMatches > nTargetPrimaries) || (nTargetNuSplits > 0)));
            const int isSplitCR(isCosmicRay && (nTargetCRMatches > 1));
            const int isLost(nTargetMatches == 0);

            std::stringstream outcomeSS;
            outcomeSS << LArInteractionTypeHelper::ToString(interactionType) << " (Nuance " << mcNuanceCode << ", Nu " << isBeamNeutrinoFinalState << ", TB " << isBeamParticle << ", CR " << isCosmicRay << ")" << std::endl;

            if (isLastNeutrinoPrimary) ++nTotalNu;
            if (isBeamParticle) ++nTotalTB;
            if (isCosmicRay) ++nTotalCR;
            if (isCorrectNu) ++nCorrectNu;
            if (isCorrectTB) ++nCorrectTB;
            if (isCorrectCR) ++nCorrectCR;
            if (isFakeNu) ++nFakeNu;
            if (isFakeCR) ++nFakeCR;
            if (isSplitNu) ++nSplitNu;
            if (isSplitCR) ++nSplitCR;
            if (isLost) ++nLost;

            if (isCorrectNu) outcomeSS << "IsCorrectNu ";
            if (isCorrectTB) outcomeSS << "IsCorrectTB ";
            if (isCorrectCR) outcomeSS << "IsCorrectCR ";
            if (isFakeNu) outcomeSS << "IsFakeNu ";
            if (isFakeCR) outcomeSS << "IsFakeCR ";
            if (isSplitNu) outcomeSS << "isSplitNu ";
            if (isSplitCR) outcomeSS << "IsSplitCR ";
            if (isLost) outcomeSS << "IsLost ";
            if (nTargetNuMatches > 0) outcomeSS << "(NNuMatches: " << nTargetNuMatches << ") ";
            if (nTargetNuLosses > 0) outcomeSS << "(NNuLosses: " << nTargetNuLosses << ") ";
            if (nTargetNuSplits > 0) outcomeSS << "(NNuSplits: " << nTargetNuSplits << ") ";
            if (nTargetCRMatches > 0) outcomeSS << "(NCRMatches: " << nTargetCRMatches << ") ";
            if (printToScreen) std::cout << outcomeSS.str() << std::endl << targetSS.str() << std::endl;

            if (fillTree)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "interactionType", interactionTypeInt));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectNu", isCorrectNu));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectTB", isCorrectTB));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectCR", isCorrectCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isFakeNu", isFakeNu));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isFakeCR", isFakeCR));
                if (!m_testBeamMode)
                {
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isSplitNu", isSplitNu));
                }
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isSplitCR", isSplitCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isLost", isLost));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
            }

            targetSS.str(std::string()); targetSS.clear();
            recoNeutrinos.clear(); associatedMCPrimaries.clear();
            nTargetMatches = 0; nTargetNuMatches = 0; nTargetCRMatches = 0; nTargetGoodNuMatches = 0; nTargetNuSplits = 0; nTargetNuLosses = 0;
            mcPrimaryId.clear(); mcPrimaryPdg.clear(); nMCHitsTotal.clear(); nMCHitsU.clear(); nMCHitsV.clear(); nMCHitsW.clear();
            mcPrimaryE.clear(); mcPrimaryPX.clear(); mcPrimaryPY.clear(); mcPrimaryPZ.clear();
            mcPrimaryVtxX.clear(); mcPrimaryVtxY.clear(); mcPrimaryVtxZ.clear(); mcPrimaryEndX.clear(); mcPrimaryEndY.clear(); mcPrimaryEndZ.clear();
            nPrimaryMatchedPfos.clear(); nPrimaryMatchedNuPfos.clear(); nPrimaryMatchedCRPfos.clear();
            bestMatchPfoId.clear(); bestMatchPfoPdg.clear(); bestMatchPfoIsRecoNu.clear(); bestMatchPfoRecoNuId.clear();
            bestMatchPfoNHitsTotal.clear(); bestMatchPfoNHitsU.clear(); bestMatchPfoNHitsV.clear(); bestMatchPfoNHitsW.clear();
            bestMatchPfoNSharedHitsTotal.clear(); bestMatchPfoNSharedHitsU.clear(); bestMatchPfoNSharedHitsV.clear(); bestMatchPfoNSharedHitsW.clear();
        }
    }

    if (fillTree)
    {
        //NEW CODE
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

        const PfoList *pPfoList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

        this->WriteEventDescription(pMCParticleList, pPfoList, pCaloHitList);
        this->WriteVariables(pPfoList);

        std::cout << "nRecoNuCosmicRays: " << nRecoNuCosmicRays << std::endl;
        std::cout << "nRecoNuParticleMatches: " << nRecoNuParticleMatches << std::endl;
        std::cout << "modifiedInteractionTypeCopy: " << modifiedInteractionTypeCopy << std::endl;

        int signal(((nPrimaries == 1 && nMuons == 1 && nRecoNuParticleMatches == 1 && nRecoNuCosmicRays == 0) || (nPrimaries == 2 && nMuons == 1 && nProtons == 1 && nRecoNuParticleMatches == 2 && nRecoNuCosmicRays == 0 && nMuonParticleMatches == 1 && nProtonParticleMatches == 1)) ? 1 : 0);

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "fileIdentifier", m_fileIdentifier));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "eventNumber", m_eventNumber - 1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NeutrinoNuanceCode", nuanceCodeCopy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "TrueInteractionType", interactionTypeCopy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "ModifiedInteractionType", modifiedInteractionTypeCopy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "TrueNeutrinoNumberAssociatedParticles", nNuParticleMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "TrueNeutrinoNumberAssociatedTracks", nNuTrackMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "TrueMuonNumberAssociatedParticles", nMuonParticleMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "TrueMuonNumberAssociatedTracks", nMuonTrackMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "TrueProtonNumberAssociatedParticles", nProtonParticleMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "TrueProtonNumberAssociatedTracks", nProtonTrackMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "RecoNeutrinoNumberAssociatedParticles", nRecoNuParticleMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "RecoNeutrinoNumberGoodAssociatedParticles", nRecoNuGoodParticleMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "RecoNeutrinoNumberAssociatedTracks", nRecoNuTrackMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "RecoNeutrinoNumberSplits", nRecoNuSplits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "RecoNeutrinoNumberCosmicRays", nRecoNuCosmicRays));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NumberTruePrimaries", nPrimaries));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NumberTrueMuons", nMuons));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NumberTrueProtons", nProtons));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NumberTrueNeutrons", nNeutrons));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NumberTrueElectrons", nElectrons));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NumberTruePhotons", nPhotons));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NumberTruePiPlus", nPiPlus));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NumberTruePiMinus", nPiMinus));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NumberTruePiZero", nPiZero));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NumberTrueOther", nOther));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "Signal", signal));

        PANDORA_MONITORING_API(FillTree(this->GetPandora(), "EventSelection"));
    }

    /*
    std::cout << "ModifiedInteractionType: " << modifiedInteractionTypeCopy << std::endl;
    std::cout << "nMuonParticleMatches: " << nMuonParticleMatches << std::endl;
    std::cout << "nRecoNuParticleMatches: " << nRecoNuParticleMatches << std::endl;

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        if (LArMCParticleHelper::IsNeutrino(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)))
            std::cout << "Neutrino PDG: " << LArMCParticleHelper::GetParentMCParticle(pMCPrimary)->GetParticleId() << std::endl;
    }
    */

    //NEW CODE
    /*
    std::cout << "NeutrinoNuanceCode: " << nuanceCodeCopy<< std::endl;
    std::cout << "TrueInteractionType: " << interactionTypeCopy << std::endl;
    std::cout << "nNuParticleMatches: " << nNuParticleMatches << std::endl;
    std::cout << "nNuTrackMatches: " << nNuTrackMatches << std::endl;
    std::cout << "nMuonParticleMatches: " << nMuonParticleMatches << std::endl;
    std::cout << "nMuonTrackMatches: " << nMuonTrackMatches << std::endl;
    std::cout << "nProtonParticleMatches: " << nProtonParticleMatches << std::endl;
    std::cout << "nProtonTrackMatches: " << nProtonTrackMatches << std::endl;
    */

    if (useInterpretedMatching)
    {
        std::stringstream summarySS;
        summarySS << "---SUMMARY--------------------------------------------------------------------------------------" << std::endl;
        if (nTotalNu > 0) summarySS << "#CorrectNu: " << nCorrectNu << "/" << nTotalNu << ", Fraction: " << (nTotalNu > 0 ? static_cast<float>(nCorrectNu) / static_cast<float>(nTotalNu) : 0.f) << std::endl;
        if (nTotalTB > 0) summarySS << "#CorrectTB: " << nCorrectTB << "/" << nTotalTB << ", Fraction: " << (nTotalTB > 0 ? static_cast<float>(nCorrectTB) / static_cast<float>(nTotalTB) : 0.f) << std::endl;
        if (nTotalCR > 0) summarySS << "#CorrectCR: " << nCorrectCR << "/" << nTotalCR << ", Fraction: " << (nTotalCR > 0 ? static_cast<float>(nCorrectCR) / static_cast<float>(nTotalCR) : 0.f) << std::endl;
        if (nFakeNu > 0) summarySS << "#FakeNu: " << nFakeNu << " ";
        if (nFakeCR > 0) summarySS << "#FakeCR: " << nFakeCR << " ";
        if (nSplitNu > 0) summarySS << "#SplitNu: " << nSplitNu << " ";
        if (nSplitCR > 0) summarySS << "#SplitCR: " << nSplitCR << " ";
        if (nLost > 0) summarySS << "#Lost: " << nLost << " ";
        if (nFakeNu || nFakeCR || nSplitNu || nSplitCR || nLost) summarySS << std::endl;
        if (printToScreen) std::cout << summarySS.str();
    }

    if (printToScreen) std::cout << "------------------------------------------------------------------------------------------------" << std::endl << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::InterpretMatching(const ValidationInfo &validationInfo, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetAllMCParticleToHitsMap()}, mcPrimaryVector);

    PfoSet usedPfos;
    while (this->GetStrongestPfoMatch(validationInfo, mcPrimaryVector, usedPfos, interpretedMCToPfoHitSharingMap)) {}
    this->GetRemainingPfoMatches(validationInfo, mcPrimaryVector, usedPfos, interpretedMCToPfoHitSharingMap);

    // Ensure all primaries have an entry, and sorting is as desired
    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        LArMCParticleHelper::PfoToSharedHitsVector &pfoHitPairs(interpretedMCToPfoHitSharingMap[pMCPrimary]);
        std::sort(pfoHitPairs.begin(), pfoHitPairs.end(), [] (const LArMCParticleHelper::PfoCaloHitListPair &a, const LArMCParticleHelper::PfoCaloHitListPair &b) -> bool {
            return ((a.second.size() != b.second.size()) ? a.second.size() > b.second.size() : LArPfoHelper::SortByNHits(a.first, b.first)); });
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::GetStrongestPfoMatch(const ValidationInfo &validationInfo, const MCParticleVector &mcPrimaryVector,
    PfoSet &usedPfos, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    const MCParticle *pBestMCParticle(nullptr);
    LArMCParticleHelper::PfoCaloHitListPair bestPfoHitPair(nullptr, CaloHitList());

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        if (interpretedMCToPfoHitSharingMap.count(pMCPrimary))
            continue;

        if (!m_useSmallPrimaries && !validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary))
            continue;

        if (!validationInfo.GetMCToPfoHitSharingMap().count(pMCPrimary))
            continue;

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : validationInfo.GetMCToPfoHitSharingMap().at(pMCPrimary))
        {
            if (usedPfos.count(pfoToSharedHits.first))
                continue;

            if (!this->IsGoodMatch(validationInfo.GetAllMCParticleToHitsMap().at(pMCPrimary), validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first), pfoToSharedHits.second))
                continue;

            if (pfoToSharedHits.second.size() > bestPfoHitPair.second.size())
            {
                pBestMCParticle = pMCPrimary;
                bestPfoHitPair = pfoToSharedHits;
            }
        }
    }

    if (!pBestMCParticle || !bestPfoHitPair.first)
        return false;

    interpretedMCToPfoHitSharingMap[pBestMCParticle].push_back(bestPfoHitPair);
    usedPfos.insert(bestPfoHitPair.first);
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetRemainingPfoMatches(const ValidationInfo &validationInfo, const MCParticleVector &mcPrimaryVector,
    const PfoSet &usedPfos, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        if (!m_useSmallPrimaries && !validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary))
            continue;

        if (!validationInfo.GetMCToPfoHitSharingMap().count(pMCPrimary))
            continue;

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : validationInfo.GetMCToPfoHitSharingMap().at(pMCPrimary))
        {
            if (usedPfos.count(pfoToSharedHits.first))
                continue;

            const LArMCParticleHelper::MCParticleCaloHitListPair mcParticleToHits(pMCPrimary, pfoToSharedHits.second);
            LArMCParticleHelper::PfoToMCParticleHitSharingMap::iterator iter(pfoToMCParticleHitSharingMap.find(pfoToSharedHits.first));

            if (pfoToMCParticleHitSharingMap.end() == iter)
            {
                pfoToMCParticleHitSharingMap[pfoToSharedHits.first].push_back(mcParticleToHits);
            }
            else
            {
                if (1 != iter->second.size())
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                LArMCParticleHelper::MCParticleCaloHitListPair &originalMCParticleToHits(iter->second.at(0));

                if (mcParticleToHits.second.size() > originalMCParticleToHits.second.size())
                    originalMCParticleToHits = mcParticleToHits;
            }
        }
    }

    for (const auto &mapEntry : pfoToMCParticleHitSharingMap)
    {
        const LArMCParticleHelper::MCParticleCaloHitListPair &mcParticleToHits(mapEntry.second.at(0));
        interpretedMCToPfoHitSharingMap[mcParticleToHits.first].push_back(LArMCParticleHelper::PfoCaloHitListPair(mapEntry.first, mcParticleToHits.second));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::IsGoodMatch(const CaloHitList &trueHits, const CaloHitList &recoHits, const CaloHitList &sharedHits) const
{
    const float purity((recoHits.size() > 0) ? static_cast<float>(sharedHits.size()) / static_cast<float>(recoHits.size()) : 0.f);
    const float completeness((trueHits.size() > 0) ? static_cast<float>(sharedHits.size()) / static_cast<float>(trueHits.size()) : 0.f);

    return ((sharedHits.size() >= m_matchingMinSharedHits) && (purity >= m_matchingMinPurity) && (completeness >= m_matchingMinCompleteness));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::WriteEventDescription(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList, const pandora::CaloHitList* pCaloHitList) const
{
    LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, LArMCParticleHelper::PrimaryParameters(),
        LArMCParticleHelper::IsBeamNeutrinoFinalState, nuMCParticlesToGoodHitsMap);

    int trueInteractionType(this->GetInteractionType(nuMCParticlesToGoodHitsMap));

    pandora::MCParticleVector trueNeutrinos;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);

    if (trueNeutrinos.size() != 1) return;

    //int neutrinoNuanceCode(LArMCParticleHelper::GetNuanceCode(trueNeutrinos.front()));

    pandora::PfoList trueNeutrinoPfos(this->GetTrueNeutrinoAssociatedPfos(pPfoList));  
    int nTrueNeutrinoAssociatedTracks(this->CountTrackPfos(trueNeutrinoPfos));

    int firstNeutrinoPfoMCPDG(-1), secondNeutrinoPfoMCPDG(-1);
    int firstNeutrinoPfoNumberHits(-1), secondNeutrinoPfoNumberHits(-1);
    int firstNeutrinoPfoParentMCPDG(-1), secondNeutrinoPfoParentMCPDG(-1);

    try
    {
        if (trueNeutrinoPfos.size() > 0)
        {
            firstNeutrinoPfoMCPDG = LArMCParticleHelper::GetMainMCParticle(trueNeutrinoPfos.front())->GetParticleId();
            firstNeutrinoPfoNumberHits = this->CountPfoHits(trueNeutrinoPfos.front());
            firstNeutrinoPfoParentMCPDG = -1;

            if ((LArMCParticleHelper::GetMainMCParticle(trueNeutrinoPfos.front()))->GetParentList().size() > 0)
                firstNeutrinoPfoParentMCPDG = (LArMCParticleHelper::GetMainMCParticle(trueNeutrinoPfos.front()))->GetParentList().front()->GetParticleId();
        }

        if (trueNeutrinoPfos.size() > 1)
        {
            secondNeutrinoPfoMCPDG = LArMCParticleHelper::GetMainMCParticle(*(std::next(trueNeutrinoPfos.begin())))->GetParticleId();
            secondNeutrinoPfoNumberHits = this->CountPfoHits(*(std::next(trueNeutrinoPfos.begin())));
            secondNeutrinoPfoParentMCPDG = -1;

            if ((LArMCParticleHelper::GetMainMCParticle(*(std::next(trueNeutrinoPfos.begin()))))->GetParentList().size() > 0)
                secondNeutrinoPfoParentMCPDG = (LArMCParticleHelper::GetMainMCParticle(*(std::next(trueNeutrinoPfos.begin()))))->GetParentList().front()->GetParticleId();
        }
    }
    catch (...)
    {
        std::cout << "Could not determine neutrino PFO PDGs" << std::endl;
    }

    PfoList neutrinoPfos;
    LArPfoHelper::GetRecoNeutrinos(pPfoList, neutrinoPfos);

    if (neutrinoPfos.size() != 1) return;

    pandora::PfoList recoNeutrinoPrimaryDaughters(this->GetPrimaryDaughters(neutrinoPfos));
    //int nRecoNeutrinoAssociatedTracks(this->CountTrackPfos(recoNeutrinoPrimaryDaughters));

    //int interactionTypeResemblance(this->GetInteractionTypeResemblance(recoNeutrinoPrimaryDaughters, nuMCParticlesToGoodHitsMap));

    //int chosenSliceContainsNeutrino(trueInteractionType == interactionTypeResemblance ? 1 : 0);
    int nSliceCosmicRays(this->GetNumberCosmicRaysChosenSlice(recoNeutrinoPrimaryDaughters));

    pandora::PfoList allConnectedPfos;
    LArPfoHelper::GetAllConnectedPfos(neutrinoPfos, allConnectedPfos);

    pandora::CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(allConnectedPfos, TPC_3D, caloHitList);
    int chosenSliceNumberHits(caloHitList.size());

    std::vector<int> trueNeutrinoHitsVector(this->GetNeutrinoInducedHits(pMCParticleList, pPfoList, recoNeutrinoPrimaryDaughters, pCaloHitList));
 
    int nTrueNeutrinoInducedHits(trueNeutrinoHitsVector.at(0)), nPfoClustersTrueNeutrinoHits(trueNeutrinoHitsVector.at(1)), nPfoIsolatedTrueNeutrinoHits(trueNeutrinoHitsVector.at(2)), nPfoAllTrueNeutrinoHits(trueNeutrinoHitsVector.at(3)), nRecoNuClustersTrueNeutrinoHits(trueNeutrinoHitsVector.at(4)),  nRecoNuIsolatedTrueNeutrinoHits(trueNeutrinoHitsVector.at(5)), nRecoNuAllTrueNeutrinoHits(trueNeutrinoHitsVector.at(6));

    bool taggingFailure(this->IsTaggingFailure(pMCParticleList, pPfoList, pCaloHitList, nRecoNuAllTrueNeutrinoHits));

    /*
    std::cout << "Reconstructed neutrino particle multiplicity: " << recoNeutrinoPrimaryDaughters.size() << std::endl;
    std::cout << "True interaction type: " << LArInteractionTypeHelper::ToString(static_cast<LArInteractionTypeHelper::InteractionType>(trueInteractionType)) << std::endl;
    std::cout << "Number of cosmic rays: " << nSliceCosmicRays << std::endl;
    std::cout << "Tagging failure: " << taggingFailure << std::endl;
    std::cout << "firstNeutrinoPfoMCPDG: " << firstNeutrinoPfoMCPDG << std::endl;
    std::cout << "firstNeutrinoPfoParentMCPDG: " << firstNeutrinoPfoParentMCPDG << std::endl;
    std::cout << "Reconstructable: " << (trueInteractionType == -1 ? 0 : 1) << std::endl;
    */

    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "TrueInteractionType", trueInteractionType));
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NeutrinoNuanceCode", neutrinoNuanceCode));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "Reconstructable", (trueInteractionType == -1 ? 0 : 1)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "FirstNeutrinoPfoMCPDG", firstNeutrinoPfoMCPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "SecondNeutrinoPfoMCPDG", secondNeutrinoPfoMCPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "FirstNeutrinoPfoNumberHits", firstNeutrinoPfoNumberHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "SecondNeutrinoPfoNumberHits", secondNeutrinoPfoNumberHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "FirstNeutrinoPfoParentMCPDG", firstNeutrinoPfoParentMCPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "SecondNeutrinoPfoParentMCPDG", secondNeutrinoPfoParentMCPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "TrueNeutrinoNumberAssociatedParticles", (int)trueNeutrinoPfos.size()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "TrueNeutrinoNumberAssociatedTracks", nTrueNeutrinoAssociatedTracks));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "TrueNeutrinoNumberInducedHits", nTrueNeutrinoInducedHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "PfoClustersTrueNeutrinoHits", nPfoClustersTrueNeutrinoHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "PfoIsolatedTrueNeutrinoHits", nPfoIsolatedTrueNeutrinoHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "PfoAllTrueNeutrinoHits", nPfoAllTrueNeutrinoHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "RecoNuClustersTrueNeutrinoHits", nRecoNuClustersTrueNeutrinoHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "RecoNuIsolatedTrueNeutrinoHits", nRecoNuIsolatedTrueNeutrinoHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "RecoNuAllTrueNeutrinoHits", nRecoNuAllTrueNeutrinoHits));
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "RecoNeutrinoNumberAssociatedParticles", (int)recoNeutrinoPrimaryDaughters.size()));
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "RecoNeutrinoNumberAssociatedTracks", nRecoNeutrinoAssociatedTracks));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "TaggingFailure", taggingFailure ? 1 : 0));
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "ChosenSliceInteractionType", interactionTypeResemblance));
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "ChosenSliceContainsTrueNeutrino", chosenSliceContainsNeutrino));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "ChosenSliceNumberCosmicRays", nSliceCosmicRays));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "ChosenSliceNumberHits", chosenSliceNumberHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "FileIdentifier", m_fileIdentifier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "EventNumber", m_eventNumber));
    //PANDORA_MONITORING_API(FillTree(this->GetPandora(), "EventSelection"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

int EventValidationAlgorithm::GetInteractionType(LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap) const
{
    MCParticleList mcPrimaryList;
    for (const auto &mapEntry : nuMCParticlesToGoodHitsMap) mcPrimaryList.push_back(mapEntry.first);
    mcPrimaryList.sort(LArMCParticleHelper::SortByMomentum);

    if (mcPrimaryList.size() == 0)
        return -1;

    const LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(mcPrimaryList));
    return static_cast<int>(static_cast<unsigned int>(interactionType));
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<int> EventValidationAlgorithm::GetNeutrinoInducedHits(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList, pandora::PfoList &recoNeutrinoPrimaryDaughters, const pandora::CaloHitList* pCaloHitList) const
{
    LArMCParticleHelper::MCRelationMap mcPrimaryMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);

    LArMCParticleHelper::CaloHitToMCMap hitToMCMap;
    LArMCParticleHelper::MCContributionMap mcToTrueHitListMap; 
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcPrimaryMap, hitToMCMap, mcToTrueHitListMap);

    pandora::PfoList allPfos;

    for (const auto pPfo : *pPfoList)
        allPfos.insert(allPfos.begin(), pPfo);

    pandora::CaloHitList allCaloHits, allTrueNeutrinoHits;

    for (const auto pCaloHit : *pCaloHitList)
        allCaloHits.insert(allCaloHits.begin(), pCaloHit);

    allTrueNeutrinoHits = this->GetTrueNeutrinoHits(allCaloHits, hitToMCMap);

    pandora::CaloHitList pfosAllHits, pfosClusterHits, pfosIsolatedHits;
    pandora::CaloHitList recoNuAllHits, recoNuClusterHits, recoNuIsolatedHits;

    this->GetClusterHits(allPfos, pfosClusterHits);
    this->GetIsolatedHits(allPfos, pfosIsolatedHits);
    this->MergeHitLists(pfosClusterHits, pfosIsolatedHits, pfosAllHits);

    this->GetClusterHits(recoNeutrinoPrimaryDaughters, recoNuClusterHits);
    this->GetIsolatedHits(recoNeutrinoPrimaryDaughters, recoNuIsolatedHits);
    this->MergeHitLists(recoNuClusterHits, recoNuIsolatedHits, recoNuAllHits);

    int nTrueNeutrinoInducedHits(this->GetTrueNeutrinoHits(allCaloHits, hitToMCMap).size()); 

    int nPfoClustersTrueNeutrinoHits(this->GetTrueNeutrinoHits(pfosClusterHits, hitToMCMap).size()), nPfoIsolatedTrueNeutrinoHits(this->GetTrueNeutrinoHits(pfosIsolatedHits, hitToMCMap).size()), nPfoAllTrueNeutrinoHits(this->GetTrueNeutrinoHits(pfosAllHits, hitToMCMap).size());

    int nRecoNuClustersTrueNeutrinoHits(this->GetTrueNeutrinoHits(recoNuClusterHits, hitToMCMap).size()), nRecoNuIsolatedTrueNeutrinoHits(this->GetTrueNeutrinoHits(recoNuIsolatedHits, hitToMCMap).size()), nRecoNuAllTrueNeutrinoHits(this->GetTrueNeutrinoHits(recoNuAllHits, hitToMCMap).size());

    /*
    std::cout << "------------------------------------" << std::endl;
    std::cout << "nTrueNeutrinoInducedHits: " << nTrueNeutrinoInducedHits << std::endl;
    //std::cout << "nPfoClustersTrueNeutrinoHits: " << nPfoClustersTrueNeutrinoHits << std::endl;
    //std::cout << "nPfoIsolatedTrueNeutrinoHits: " << nPfoIsolatedTrueNeutrinoHits << std::endl;
    std::cout << "nPfoAllTrueNeutrinoHits: " << nPfoAllTrueNeutrinoHits << std::endl;
    //std::cout << "nRecoNuClustersTrueNeutrinoHits: " << nRecoNuClustersTrueNeutrinoHits << std::endl;
    //std::cout << "nRecoNuIsolatedTrueNeutrinoHits: " << nRecoNuIsolatedTrueNeutrinoHits << std::endl;
    std::cout << "nRecoNuAllTrueNeutrinoHits: " << nRecoNuAllTrueNeutrinoHits << std::endl;
    std::cout << "------------------------------------" << std::endl;
    */

   if (m_viewEvent)
   {
       pandora::CaloHitList allCaloHitsToDraw, trueNeutrinoHitsToDraw, recoNuHitsToDraw;

        for (const auto pCaloHit : *pCaloHitList)
        {
            if (pCaloHit->GetHitType() != TPC_VIEW_W)
                continue;

            allCaloHitsToDraw.insert(allCaloHitsToDraw.begin(), pCaloHit);

            try
            {
                const pandora::MCParticle* pMCParticle(hitToMCMap.at(pCaloHit));
                const pandora::MCParticle* pParentMCParticle = LArMCParticleHelper::GetParentMCParticle(pMCParticle);

                if (LArMCParticleHelper::IsNeutrino(pParentMCParticle))
                    trueNeutrinoHitsToDraw.insert(trueNeutrinoHitsToDraw.begin(), pCaloHit);
            }
            catch (...)
            {
                continue;
            }
        }

       LArPfoHelper::GetCaloHits(recoNeutrinoPrimaryDaughters, TPC_VIEW_W, recoNuHitsToDraw);
       LArPfoHelper::GetIsolatedCaloHits(recoNeutrinoPrimaryDaughters, TPC_VIEW_W, recoNuHitsToDraw);

       PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &allCaloHitsToDraw, "CaloHits", BLACK)); 
       PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trueNeutrinoHitsToDraw, "TrueNeutrinoHits", RED)); 
       PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &recoNuHitsToDraw, "RecoNeutrinoHits", BLUE)); 
       PANDORA_MONITORING_API(ViewEvent(this->GetPandora())); 
   }

    std::vector<int> trueNeutrinoHitsVector = {nTrueNeutrinoInducedHits, nPfoClustersTrueNeutrinoHits, nPfoIsolatedTrueNeutrinoHits, nPfoAllTrueNeutrinoHits, nRecoNuClustersTrueNeutrinoHits, nRecoNuIsolatedTrueNeutrinoHits, nRecoNuAllTrueNeutrinoHits};

    return trueNeutrinoHitsVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetClusterHits(pandora::PfoList &pfoList, pandora::CaloHitList &caloHitList) const
{
    LArPfoHelper::GetCaloHits(pfoList, TPC_VIEW_U, caloHitList);
    LArPfoHelper::GetCaloHits(pfoList, TPC_VIEW_V, caloHitList);
    LArPfoHelper::GetCaloHits(pfoList, TPC_VIEW_W, caloHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetIsolatedHits(pandora::PfoList &pfoList, pandora::CaloHitList &caloHitList) const
{
    LArPfoHelper::GetIsolatedCaloHits(pfoList, TPC_VIEW_U, caloHitList);
    LArPfoHelper::GetIsolatedCaloHits(pfoList, TPC_VIEW_V, caloHitList);
    LArPfoHelper::GetIsolatedCaloHits(pfoList, TPC_VIEW_W, caloHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::MergeHitLists(pandora::CaloHitList &clusterHitList, pandora::CaloHitList &isolatedHitList, pandora::CaloHitList &combinedHitList) const
{
    combinedHitList = clusterHitList;
    pandora::CaloHitList isolatedCaloHitsCopy = isolatedHitList;
    combinedHitList.splice(combinedHitList.begin(), isolatedCaloHitsCopy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CaloHitList EventValidationAlgorithm::GetTrueNeutrinoHits(pandora::CaloHitList &hitList, LArMCParticleHelper::CaloHitToMCMap &hitToMCMap) const
{
    pandora::CaloHitList trueNeutrinoHits;

    for (const auto pCaloHit : hitList)
    {
        try
        {
            const pandora::MCParticle* pMCParticle(hitToMCMap.at(pCaloHit));
            const pandora::MCParticle* pParentMCParticle = LArMCParticleHelper::GetParentMCParticle(pMCParticle);

            if (LArMCParticleHelper::IsNeutrino(pParentMCParticle))
                trueNeutrinoHits.insert(trueNeutrinoHits.begin(), pCaloHit);
        }
        catch (...)
        {
            continue;
        }
    }

    return trueNeutrinoHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::IsTaggingFailure(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList, const pandora::CaloHitList* pCaloHitList, int nRecoNuAllTrueNeutrinoHits) const 
{
    LArMCParticleHelper::MCRelationMap mcPrimaryMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);

    LArMCParticleHelper::CaloHitToMCMap hitToMCMap;
    LArMCParticleHelper::MCContributionMap mcToTrueHitListMap; 
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcPrimaryMap, hitToMCMap, mcToTrueHitListMap);

    bool taggingFailure(false), vocal(false);

    std::cout << "------------------------------------" << std::endl;

    for (const auto pPfo : *pPfoList)
    {
        pandora::PfoList thisPfo;
        thisPfo.insert(thisPfo.begin(), pPfo);

        if (vocal)
        {
            if (LArPfoHelper::IsNeutrino(pPfo))
                std::cout << "PFO hypothesis: NEUTRINO" << std::endl;
            else
                std::cout << "PFO hypothesis: COSMIC RAY" << std::endl;
        }

        pandora::CaloHitList thisPfoAllHits, thisPfoClusterHits, thisPfoIsolatedHits;

        this->GetClusterHits(thisPfo, thisPfoClusterHits);
        this->GetIsolatedHits(thisPfo, thisPfoIsolatedHits);
        this->MergeHitLists(thisPfoClusterHits, thisPfoIsolatedHits, thisPfoAllHits);
        
        int nThisPfoClustersTrueNeutrinoHits(this->GetTrueNeutrinoHits(thisPfoClusterHits, hitToMCMap).size()), nThisPfoIsolatedTrueNeutrinoHits(this->GetTrueNeutrinoHits(thisPfoIsolatedHits, hitToMCMap).size()), nThisPfoAllTrueNeutrinoHits(this->GetTrueNeutrinoHits(thisPfoAllHits, hitToMCMap).size());

        if (!LArPfoHelper::IsNeutrino(pPfo) && nThisPfoAllTrueNeutrinoHits > nRecoNuAllTrueNeutrinoHits)
            taggingFailure = true;

        if (vocal)
        {
            std::cout << "This PFO total number true nu hits: " << nThisPfoClustersTrueNeutrinoHits << std::endl;
            std::cout << "This PFO number isolated true nu hits: " << nThisPfoIsolatedTrueNeutrinoHits << std::endl;
            std::cout << "This PFO number cluster true nu hits: " << nThisPfoAllTrueNeutrinoHits << std::endl;
            std::cout << "********************************************" << std::endl;
        }
    } 

    return taggingFailure;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::PfoList EventValidationAlgorithm::GetTrueNeutrinoAssociatedPfos(const pandora::PfoList* pPfoList) const
{
    PfoList neutrinoPfos, trueNeutrinoPfos;
    LArPfoHelper::GetRecoNeutrinos(pPfoList, neutrinoPfos);

    for (const auto pPfo : *pPfoList)
    {
        if (pPfo == neutrinoPfos.front() || !LArPfoHelper::IsFinalState(pPfo))
            continue;

        const pandora::MCParticle* pMCParticle(this->GetMainMCParticle(pPfo));
        const pandora::MCParticle* pParentMCParticle = LArMCParticleHelper::GetParentMCParticle(pMCParticle);

        if (LArMCParticleHelper::IsNeutrino(pParentMCParticle))
            trueNeutrinoPfos.insert(trueNeutrinoPfos.begin(), pPfo);
    }

    return trueNeutrinoPfos;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int EventValidationAlgorithm::CountTrackPfos(pandora::PfoList &pfoList) const
{
    int nTrackPfos(0);    

    for (const auto pPfo : pfoList)
    {
        if (LArPfoHelper::IsTrack(pPfo))
            ++nTrackPfos;
    }

    return nTrackPfos;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::PfoList EventValidationAlgorithm::GetPrimaryDaughters(pandora::PfoList &neutrinoPfos) const
{
    pandora::PfoList allConnectedPfos, primaryDaughters;
    LArPfoHelper::GetAllConnectedPfos(neutrinoPfos, allConnectedPfos);

    for (const auto pPfo : allConnectedPfos)
    {
        if (pPfo == neutrinoPfos.front() || !LArPfoHelper::IsFinalState(pPfo))
            continue;

        primaryDaughters.insert(primaryDaughters.begin(), pPfo);
    }

    return primaryDaughters;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int EventValidationAlgorithm::GetNumberCosmicRaysChosenSlice(pandora::PfoList &primaryDaughters) const
{
    int nCosmicRays(0);

    //pandora::PfoList allConnectedPfos;
    //LArPfoHelper::GetAllConnectedPfos(neutrinoPfos, allConnectedPfos);

    for (const auto pPfo : primaryDaughters)
    {
        //if (pPfo == neutrinoPfos.front() || !LArPfoHelper::IsFinalState(pPfo))
        //    continue;

        const pandora::MCParticle* pMCParticle(this->GetMainMCParticle(pPfo));
        const pandora::MCParticle* pParentMCParticle = LArMCParticleHelper::GetParentMCParticle(pMCParticle);

        if (LArMCParticleHelper::IsCosmicRay(pParentMCParticle))
            ++nCosmicRays; 
    }

    return nCosmicRays;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* EventValidationAlgorithm::GetMainMCParticle(const ParticleFlowObject *const pPfo)
{
    ClusterList clusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);

    CaloHitList caloHitList;

    for (const Cluster *const pCluster : clusterList)
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    MCParticleWeightMap mcParticleWeightMap;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const MCParticleWeightMap &hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());
        MCParticleVector mcParticleVector;

        for (const MCParticleWeightMap::value_type &mapEntry : hitMCParticleWeightMap)
            mcParticleVector.push_back(mapEntry.first);

        std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            try
            {
                const MCParticle *const pMCPrimary = LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);
                mcParticleWeightMap[pMCPrimary] += hitMCParticleWeightMap.at(pMCParticle);
            }

            catch (...)
            {
                continue;
            }
        }
    }

    float bestWeight(0.f);
    const MCParticle *pBestMCParticle(nullptr);

    MCParticleVector mcParticleVector;

    for (const MCParticleWeightMap::value_type &mapEntry : mcParticleWeightMap)
        mcParticleVector.push_back(mapEntry.first);

    std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

    for (const MCParticle *const pCurrentMCParticle : mcParticleVector)
    {
        const float currentWeight(mcParticleWeightMap.at(pCurrentMCParticle));

        if (currentWeight > bestWeight)
        {
            pBestMCParticle = pCurrentMCParticle;
            bestWeight = currentWeight;
        }
    }

    if (!pBestMCParticle)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return pBestMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::WriteVariables(const pandora::PfoList* pPfoList) const
{
    PfoList neutrinoPfos, connectedPfos; 
    LArPfoHelper::GetRecoNeutrinos(pPfoList, neutrinoPfos);
    LArPfoHelper::GetAllDownstreamPfos(neutrinoPfos, connectedPfos);

    if (neutrinoPfos.size() != 1) return;

    //Write containment definitions
    this->WriteContainmentDefinitions(neutrinoPfos, connectedPfos);

    //Get primary nu reco daughters
    pandora::PfoList recoNeutrinoPrimaryDaughters(this->GetPrimaryDaughters(neutrinoPfos));

    //Event variables
    this->WriteEventVariables(recoNeutrinoPrimaryDaughters);

    //Get shortest and longest PFO
    const pandora::ParticleFlowObject* pLongestPfo(GetLongestPfo(recoNeutrinoPrimaryDaughters));
    const pandora::ParticleFlowObject* pShortestPfo(GetShortestPfo(recoNeutrinoPrimaryDaughters));

    //Topological variables
    this->WriteTopologicalVariables(pLongestPfo, "LongestPfo");
    this->WriteTopologicalVariables(pShortestPfo, "ShortestPfo");

    //Direction fit variables 
    this->WriteDirectionFitVariables(pLongestPfo, "LongestPfo");
    this->WriteDirectionFitVariables(pShortestPfo, "ShortestPfo");

    //PID variables
    this->WritePIDVariables(pLongestPfo, "LongestPfo");
    this->WritePIDVariables(pShortestPfo, "ShortestPfo");

    //Cosmic variables
    this->WriteCosmicVariables(pLongestPfo, "LongestPfo");
    this->WriteCosmicVariables(pShortestPfo, "ShortestPfo");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::WriteContainmentDefinitions(PfoList &neutrinoPfos, PfoList &connectedPfos) const
{
    //Check if nu reco contained
    float containedHits(0.f), uncontainedHits(0.f);

    pandora::CaloHitList pfoHits;
    LArPfoHelper::GetCaloHits(connectedPfos, TPC_3D, pfoHits);

    for (const auto pCaloHit : pfoHits)
    {
        CartesianVector correctedPosition(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(pCaloHit->GetPositionVector()));

        if (!this->IsInFiducialVolume(correctedPosition))
            uncontainedHits += 1.0;
        else
            containedHits += 1.0;
    }

    float containmentFraction(containedHits/(containedHits + uncontainedHits));

    bool everyPfoEndpointContained(true);

    for (const auto pPfo : connectedPfos)
    {
        if (pPfo == neutrinoPfos.front())
            continue;

        if (!this->IsContained(pPfo))
            everyPfoEndpointContained = false;
    }

    bool everyPfoSlidingFitEndpointContained(true);

    for (const auto pPfo : connectedPfos)
    {
        if (pPfo == neutrinoPfos.front())
            continue;

        if (!this->EndpointsContained(pPfo))
            everyPfoSlidingFitEndpointContained = false;
    }

    const pandora::CartesianVector recoNeutrinoVertexPosition(LArPfoHelper::GetVertex(neutrinoPfos.front())->GetPosition());
    const pandora::CartesianVector correctedRecoNeutrinoVertexPosition(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(recoNeutrinoVertexPosition));
    bool recoNeutrinoVertexContained(this->IsInFiducialVolume(correctedRecoNeutrinoVertexPosition));

    //std::cout << "correctedRecoNeutrinoVertexPosition : (" << correctedRecoNeutrinoVertexPosition.GetX() << ", " << correctedRecoNeutrinoVertexPosition.GetY() << ", " << correctedRecoNeutrinoVertexPosition.GetZ() << ")" << std::endl;

    std::cout << "NuReco hits: " << containedHits + uncontainedHits << std::endl;
    std::cout << "containmentFraction: " << containmentFraction << std::endl;
    //std::cout << "everyPfoEndpointContained: " << everyPfoEndpointContained << std::endl;
    //std::cout << "everyPfoSlidingFitEndpointContained: " << everyPfoSlidingFitEndpointContained << std::endl;
    //std::cout << "recoNeutrinoVertexContained: " << recoNeutrinoVertexContained << std::endl;

    //MC equivalent definitions
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const PfoList *pPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::MCRelationMap mcPrimaryMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);

    LArMCParticleHelper::CaloHitToMCMap hitToMCMap;
    LArMCParticleHelper::MCContributionMap mcToTrueHitListMap; 
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcPrimaryMap, hitToMCMap, mcToTrueHitListMap);

    pandora::MCParticleVector trueNeutrinos;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);

    const pandora::CartesianVector trueNeutrinoVertexPosition(trueNeutrinos.front()->GetVertex());
    bool trueNeutrinoVertexContained(this->IsInFiducialVolume(trueNeutrinoVertexPosition));
    bool bothTrueNeutrinoEndpointsContained(this->IsInFiducialVolume(trueNeutrinos.front()->GetVertex()) && this->IsInFiducialVolume(trueNeutrinos.front()->GetEndpoint()) ? true : false);

    int trueNeutrinoHits(0), trueNeutrinoContainedHits(0);

    pandora::CaloHitList allPfoHits;
    LArPfoHelper::GetCaloHits(*pPfoList, TPC_3D, allPfoHits);

    for (const auto pCaloHit : allPfoHits)
    {
        const pandora::CaloHit* pParentCaloHit(static_cast<const CaloHit*>(pCaloHit->GetParentAddress()));

        if (hitToMCMap.find(pParentCaloHit) != hitToMCMap.end())
        {
            if (!LArMCParticleHelper::IsBeamNeutrinoFinalState(hitToMCMap.at(pParentCaloHit)))
                continue;

            ++trueNeutrinoHits;

            CartesianVector correctedPosition(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(pCaloHit->GetPositionVector()));

            if (this->IsInFiducialVolume(correctedPosition))
                ++trueNeutrinoContainedHits;
        }
    }

    float trueNeutrinoContainmentFraction(static_cast<float>(trueNeutrinoContainedHits)/(trueNeutrinoHits));

    std::cout << "trueNeutrinoHits: " << trueNeutrinoHits << std::endl;
    //std::cout << "trueNeutrinoContainedHits: " << trueNeutrinoContainedHits << std::endl;
    std::cout << "trueNeutrinoContainmentFraction: " << trueNeutrinoContainmentFraction << std::endl;
    //std::cout << "bothTrueNeutrinoEndpointsContained: " << bothTrueNeutrinoEndpointsContained << std::endl;
    //std::cout << "trueNeutrinoVertexContained: " << trueNeutrinoVertexContained << std::endl;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NuTrueContained", trueNeutrinoContainmentFraction > 0.82 ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NuTrueContainmentFraction", trueNeutrinoContainmentFraction));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NuTrueEndpointsContained", bothTrueNeutrinoEndpointsContained ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NuTrueVertexContained", trueNeutrinoVertexContained ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NuRecoContained", containmentFraction > 0.82 ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NuRecoContainmentFraction", containmentFraction));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NuRecoEndpointsContained", everyPfoEndpointContained ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NuRecoSlidingFitEndpointsContained", everyPfoSlidingFitEndpointContained ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NuRecoVertexContained", recoNeutrinoVertexContained ? 1 : 0));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::WriteEventVariables(PfoList &recoNeutrinoPrimaryDaughters) const
{
    float totalEventCharge(0.f);

    for (const auto pPfo : recoNeutrinoPrimaryDaughters)
        totalEventCharge += GetPfoCharge(pPfo);

    //Get shortest and longest PFO
    const pandora::ParticleFlowObject* pLongestPfo(GetLongestPfo(recoNeutrinoPrimaryDaughters));
    const pandora::ParticleFlowObject* pShortestPfo(GetShortestPfo(recoNeutrinoPrimaryDaughters));

    float openingAngle(GetPfoOpeningAngle(pLongestPfo, pShortestPfo));

    pandora::CartesianVector neutrinoMomentum(GetApproximateNeutrinoMomentum(recoNeutrinoPrimaryDaughters, pLongestPfo));

    //Write all information to tree
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "TotalEventCharge", totalEventCharge));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "OpeningAngle", openingAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NeutrinoMomentumX", neutrinoMomentum.GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NeutrinoMomentumY", neutrinoMomentum.GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", "NeutrinoMomentumZ", neutrinoMomentum.GetZ()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::WriteTopologicalVariables(const ParticleFlowObject *const pPfo, std::string variableNamePrefix) const
{
    float pfoCharge(GetPfoCharge(pPfo));

    CartesianVector xAxis(1.f, 0.f, 0.f), yAxis(0.f, 1.f, 0.f);
    float pfoAzimuthalAngle(GetAngleWithVector(pPfo, xAxis)), pfoPolarAngle(GetAngleWithVector(pPfo, yAxis));
    float pfoTheta(GetThetaBeamPfo(pPfo));

    float pfoLength(std::sqrt(LArPfoHelper::GetThreeDLengthSquared(pPfo)));
    float pfoTrackProbability(this->CalculateTrackProbability(pPfo));
    float pfoDaughterHitFraction(this->GetDaughterHitFraction(pPfo));

    pandora::ClusterList clusterListW;
    LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, clusterListW);
    float pfoLength2D(clusterListW.size() == 1 ? LArClusterHelper::GetLength(clusterListW.front()) : -1.f);
    
    float pfo3DTo2DRatio(pfoLength/pfoLength2D);

    const pandora::Vertex *const pVertex = LArPfoHelper::GetVertex(pPfo);
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    LArTrackStateVector trackStateVector;
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, pVertex, m_slidingFitWindow, slidingFitPitch, trackStateVector);

    pandora::CartesianVector lowZVector(0.f, 0.f, 0.f), highZVector(0.f, 0.f, 0.f);

    if (trackStateVector.size() > 0)
    {
        TrackState firstTrackState(*(trackStateVector.begin())), lastTrackState(trackStateVector.back());
        pandora::CartesianVector initialPosition(firstTrackState.GetPosition());
        pandora::CartesianVector endPosition(lastTrackState.GetPosition());

        lowZVector = (initialPosition.GetZ() < endPosition.GetZ() ? initialPosition : endPosition);
        highZVector = (initialPosition.GetZ() > endPosition.GetZ() ? initialPosition : endPosition);
    }

    pandora::CartesianVector correctedLowZVector(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(lowZVector)), correctedHighZVector(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(highZVector));

    float minX(std::min(lowZVector.GetX(), highZVector.GetX())), maxX(std::max(lowZVector.GetX(), highZVector.GetX()));
    float minY(std::min(lowZVector.GetY(), highZVector.GetY())), maxY(std::max(lowZVector.GetY(), highZVector.GetY()));
    float minZ(std::min(lowZVector.GetZ(), highZVector.GetZ())), maxZ(std::max(lowZVector.GetZ(), highZVector.GetZ()));

    float correctedMinX(std::min(correctedLowZVector.GetX(), correctedHighZVector.GetX())), correctedMaxX(std::max(correctedLowZVector.GetX(), correctedHighZVector.GetX()));
    float correctedMinY(std::min(correctedLowZVector.GetY(), correctedHighZVector.GetY())), correctedMaxY(std::max(correctedLowZVector.GetY(), correctedHighZVector.GetY()));
    float correctedMinZ(std::min(correctedLowZVector.GetZ(), correctedHighZVector.GetZ())), correctedMaxZ(std::max(correctedLowZVector.GetZ(), correctedHighZVector.GetZ()));
    
    pandora::PfoList downstreamPfos;
    LArPfoHelper::GetAllDownstreamPfos(pPfo, downstreamPfos);
    int numberDaughters(downstreamPfos.size());

    bool isContained(this->IsContained(pPfo));

    const pandora::MCParticle* pMCParticle(this->GetMainMCParticle(pPfo));
    bool mcIsContained(this->IsContained(pMCParticle));

    //Write all information to tree
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCForwards", (pMCParticle->GetVertex().GetZ() < pMCParticle->GetEndpoint().GetZ()) ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCDownwards", (pMCParticle->GetVertex().GetY() > pMCParticle->GetEndpoint().GetY()) ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCContained", mcIsContained? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "Contained", isContained ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "Charge", pfoCharge));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "Length", pfoLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "TrackProbability", pfoTrackProbability));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "PolarAngle", pfoPolarAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "AzimuthalAngle", pfoAzimuthalAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "Theta", pfoTheta));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "DaughterHitFraction", pfoDaughterHitFraction));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "BeginX", lowZVector.GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "BeginY", lowZVector.GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "BeginZ", lowZVector.GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "EndX", highZVector.GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "EndY", highZVector.GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "EndZ", highZVector.GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MinX", minX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MaxX", maxX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MinY", minY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MaxY", maxY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MinZ", minZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MaxZ", maxZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "CorrectedMinX", correctedMinX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "CorrectedMaxX", correctedMaxX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "CorrectedMinY", correctedMinY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "CorrectedMaxY", correctedMaxY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "CorrectedMinZ", correctedMinZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "CorrectedMaxZ", correctedMaxZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "NumberDaughters", numberDaughters));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "3DTo2DRatio", pfo3DTo2DRatio));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::IsContained(const ParticleFlowObject *const pPfo) const
{
    pandora::CaloHitList pfoHits;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, pfoHits);

    pandora::CartesianVector trackBeginpoint(0.f, 0.f, 0.f), trackEndpoint(0.f, 0.f, 0.f);

    float maxDistanceFromVertex(0.f);

    for (const auto pCaloHit1 : pfoHits)
    {
        pandora::CartesianVector correctedHitPosition1(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(pCaloHit1->GetPositionVector()));

        for (const auto pCaloHit2 : pfoHits)
        {
            pandora::CartesianVector correctedHitPosition2(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(pCaloHit2->GetPositionVector()));

            if ((correctedHitPosition1 - correctedHitPosition2).GetMagnitude() > maxDistanceFromVertex) 
            {
                maxDistanceFromVertex = (correctedHitPosition1 - correctedHitPosition2).GetMagnitude(); 
                trackBeginpoint = correctedHitPosition1; 
                trackEndpoint = correctedHitPosition2; 
            }
        }
    }

    if (trackBeginpoint.GetMagnitude() == 0 || pfoHits.size() == 0)
        return true;

    //std::cout << "trackBeginpoint : (" << trackBeginpoint.GetX() << ", " << trackBeginpoint.GetY() << ", " << trackBeginpoint.GetZ() << ")" << std::endl;
    //std::cout << "trackEndpoint : (" << trackEndpoint.GetX() << ", " << trackEndpoint.GetY() << ", " << trackEndpoint.GetZ() << ")" << std::endl;
    //std::cout << "this->IsInFiducialVolume(trackBeginpoint): " << this->IsInFiducialVolume(trackBeginpoint) << std::endl; 
    //std::cout << "this->IsInFiducialVolume(trackEndpoint): " << this->IsInFiducialVolume(trackEndpoint) << std::endl; 

    if (this->IsInFiducialVolume(trackBeginpoint) && this->IsInFiducialVolume(trackEndpoint))
        return true;
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::EndpointsContained(const ParticleFlowObject *const pPfo) const
{
    LArTrackStateVector trackStateVector;
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, LArPfoHelper::GetVertex(pPfo), m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()), trackStateVector);

    if (trackStateVector.size() < 2)
        return true;

    const pandora::CartesianVector trackBeginpoint(trackStateVector.front().GetPosition());
    const pandora::CartesianVector trackEndpoint(trackStateVector.back().GetPosition());

    //std::cout << "trackBeginpoint : (" << trackBeginpoint.GetX() << ", " << trackBeginpoint.GetY() << ", " << trackBeginpoint.GetZ() << ")" << std::endl;
    //std::cout << "trackEndpoint : (" << trackEndpoint.GetX() << ", " << trackEndpoint.GetY() << ", " << trackEndpoint.GetZ() << ")" << std::endl;
    //std::cout << "this->IsInFiducialVolume(trackBeginpoint): " << this->IsInFiducialVolume(trackBeginpoint) << std::endl; 
    //std::cout << "this->IsInFiducialVolume(trackEndpoint): " << this->IsInFiducialVolume(trackEndpoint) << std::endl; 


    if (this->IsInFiducialVolume(trackBeginpoint) && this->IsInFiducialVolume(trackEndpoint))
        return true;
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::IsContained(const MCParticle* const pMCParticle) const
{
    CartesianVector mcVertex(pMCParticle->GetVertex()), mcEndpoint(pMCParticle->GetEndpoint());

    if (this->IsInFiducialVolume(mcVertex) && this->IsInFiducialVolume(mcEndpoint))
        return true;
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::IsInFiducialVolume(pandora::CartesianVector positionVector) const
{
    if ((positionVector.GetX() > 12.0 && positionVector.GetX() < (256.35 - 12.0)) && (positionVector.GetY() > (-116.5 + 35.0) && positionVector.GetY() < (116.5 - 35.0)) && ((positionVector.GetZ() > 25.0 && positionVector.GetZ() < 675.0) || (positionVector.GetZ() > 775.0 && positionVector.GetZ() < (1036.8 - 85.0))))
        return true;
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EventValidationAlgorithm::GetDaughterHitFraction(const ParticleFlowObject *const pPfo) const
{
    pandora::PfoList thisPfoList, downstreamPfos;
    thisPfoList.insert(thisPfoList.begin(), pPfo);
    this->GetAllDownstreamPfos(thisPfoList, downstreamPfos);

    pandora::CaloHitList pfoHits, downstreamPfoHits;

    this->GetClusterHits(thisPfoList, pfoHits);
    this->GetClusterHits(downstreamPfos, downstreamPfoHits);

    float numberDownstreamHits(downstreamPfoHits.size()), numberPrimaryHits(pfoHits.size());
    float daughterHitFraction(numberDownstreamHits/numberPrimaryHits);

    return daughterHitFraction; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetAllDownstreamPfos(const PfoList &inputPfoList, PfoList &outputPfoList) const
{
    for (const ParticleFlowObject *const pPfo : inputPfoList)
        this->GetAllDownstreamPfos(pPfo, outputPfoList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetAllDownstreamPfos(const ParticleFlowObject *const pPfo, PfoList &outputPfoList) const
{
    if (outputPfoList.end() != std::find(outputPfoList.begin(), outputPfoList.end(), pPfo))
        return;

    LArPfoHelper::GetAllDownstreamPfos(pPfo->GetDaughterPfoList(), outputPfoList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::WriteDirectionFitVariables(const ParticleFlowObject *const pPfo, std::string variableNamePrefix) const
{
    try
    {
        TrackDirectionTool::DirectionFitObject fitResult = m_pTrackDirectionTool->GetPfoDirection(pPfo);
        TrackDirectionTool::FitParameters fitParameters(fitResult.GetFitParameters());
        TrackDirectionTool::HitChargeVector bestFitCharges(fitResult.GetDeltaChiSquaredPerHit() < 0 ? fitResult.GetForwardsFitCharges() : fitResult.GetBackwardsFitCharges());

        //fitResult.DrawFit(); 

        //Primitive particle ID variables
        float meanQW(0.f), meanHitCharge(0.f);

        for (const auto &fitCharge : bestFitCharges)
            meanQW += fitCharge.GetChargeOverWidth();

        for (const auto &hitCharge : fitResult.GetHitChargeVector())
            meanHitCharge += hitCharge.GetChargeOverWidth();

        meanQW /= bestFitCharges.size();
        meanHitCharge /= fitResult.GetHitChargeVector().size();

        float fitChargeRange(std::abs(bestFitCharges.front().GetChargeOverWidth() - bestFitCharges.back().GetChargeOverWidth()));
        float trackLength(std::abs(bestFitCharges.front().GetLongitudinalPosition() - bestFitCharges.back().GetLongitudinalPosition()));
        float fitChargeRangeOverLength(trackLength > 1 ? (fitChargeRange/trackLength) : -1.f);

        float angleXAxis(GetAngleXAxis(pPfo));
        float angleCorrectedFitChargeRangeOverLength(std::abs(fitChargeRangeOverLength * cos(angleXAxis)));

        float leftSumQW(bestFitCharges.at(0).GetChargeOverWidth() + bestFitCharges.at(1).GetChargeOverWidth() + bestFitCharges.at(2).GetChargeOverWidth()), rightSumQW(bestFitCharges.at(bestFitCharges.size() - 3).GetChargeOverWidth() + bestFitCharges.at(bestFitCharges.size() - 2).GetChargeOverWidth() + bestFitCharges.at(bestFitCharges.size() - 1).GetChargeOverWidth());

        float minSumQW(std::min(leftSumQW, rightSumQW)), maxSumQW(std::max(leftSumQW, rightSumQW));
        float sumQWRatio(maxSumQW/minSumQW);

        //Straightness ratio
        float pfoLength(std::sqrt(LArPfoHelper::GetThreeDLengthSquared(pPfo)));
        float straightnessRatio((fitResult.GetBeginpoint() - fitResult.GetEndpoint()).GetMagnitude() > 1 ? (pfoLength/(fitResult.GetBeginpoint() - fitResult.GetEndpoint()).GetMagnitude()) : -1.f);

        //Cosmic downwards probability from previous work
        float cosmicProbability(this->CalculateCosmicProbability(fitResult));
        //std::cout << "cosmicProbability: " << cosmicProbability << std::endl; 

        int intersectsYFace(this->IntersectsYFace(fitResult) ? 1 : 0); 

        /*
        if (variableNamePrefix == "LongestPfo")
        {
            std::cout << "DeltaChiSquaredPerHit: " << fitResult.GetDeltaChiSquaredPerHit() << std::endl; 
            std::cout << "Longest PFO UpDownDeltaChiSquaredPerHit: " << fitResult.GetUpDownDeltaChiSquaredPerHit() << std::endl;
        }
        */

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "DirectionFitSuccesful", 1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "FitParameterZero", fitParameters.GetParameterZero()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "FitParameterOne", fitParameters.GetParameterOne()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "FitParameterTwo", fitParameters.GetParameterTwo()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "UpDownDeltaChiSquaredPerHit", fitResult.GetUpDownDeltaChiSquaredPerHit()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "DeltaChiSquaredPerHit", fitResult.GetDeltaChiSquaredPerHit()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MinChiSquaredPerHit", fitResult.GetMinChiSquaredPerHit()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "SplitApplied", (fitResult.GetSplitObject().GetSplitApplied() ? 1 : 0)));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "FitChargeRangeOverLength", fitChargeRangeOverLength));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "AngleCorrectedFitChargeRangeOverLength", angleCorrectedFitChargeRangeOverLength));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "CosmicProbability", cosmicProbability));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MinSumQW", minSumQW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MaxSumQW", maxSumQW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "SumQWRatio", sumQWRatio));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MeanQW", meanQW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MeanHitCharge", meanHitCharge));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "StraightnessRatio", straightnessRatio));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "IntersectsYFace", intersectsYFace));
    }
    catch (...)
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "DirectionFitSuccesful", 0));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "FitParameterZero", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "FitParameterOne", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "FitParameterTwo", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "UpDownDeltaChiSquaredPerHit", -100.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "DeltaChiSquaredPerHit", -100.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MinChiSquaredPerHit", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "SplitApplied", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "FitChargeRangeOverLength", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "AngleCorrectedFitChargeRangeOverLength", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "CosmicProbability", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MinSumQW", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MaxSumQW", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "SumQWRatio", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MeanQW", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MeanHitCharge", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "StraightnessRatio", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "IntersectsYFace", -1));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::WriteCosmicVariables(const ParticleFlowObject *const pPfo, std::string variableNamePrefix) const
{
    const pandora::MCParticle* pMCParticle(this->GetMainMCParticle(pPfo));
    CartesianVector mcVertex(pMCParticle->GetVertex()), mcEndpoint(pMCParticle->GetEndpoint());
    CartesianVector mcLowYPosition(mcVertex.GetY() <= mcEndpoint.GetY() ? mcVertex : mcEndpoint);

    LArTrackStateVector trackStateVector;
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, LArPfoHelper::GetVertex(pPfo), m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()), trackStateVector);

    if (trackStateVector.size() == 0)
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCDeltaY", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCVertexY", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoDeltaY", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoVertexY", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCTrueUpwards", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCPDG", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "AbsMCPDG", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCCosmicRay", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCFiducialLowY", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoFiducialLowY", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCIntersectsYFace", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoIntersectsYFace", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCLowestTenCmTotalCharge", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoLowestTenCmTotalCharge", -1.f));
        return;
    }

    std::sort(trackStateVector.begin(), trackStateVector.end(), [](TrackState pTrackState1, TrackState pTrackState2) -> bool { return pTrackState1.GetPosition().GetY() < pTrackState2.GetPosition().GetY(); });

    pandora::CartesianVector recoLowYPosition(trackStateVector.front().GetPosition());
    pandora::CartesianVector recoHighYPosition(trackStateVector.back().GetPosition());

    /*
    pandora::CartesianVector recoLowYPosition(0.f, 1e6, 0.f);
    pandora::CartesianVector recoHighYPosition(0.f, 0.f, 0.f);

    pandora::CaloHitList pfoHits;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, pfoHits);

    for (const auto pCaloHit : pfoHits)
    {
        if (pCaloHit->GetPositionVector().GetY() < recoLowYPosition.GetY())
            recoLowYPosition = pCaloHit->GetPositionVector();

        if (pCaloHit->GetPositionVector().GetY() > recoHighYPosition.GetY())
            recoHighYPosition = pCaloHit->GetPositionVector();
    }
    */

    ClusterList clusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_3D, clusterList);
    ClusterVector clusterVector(clusterList.begin(), clusterList.end());

    if (clusterVector.size() == 0)
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCDeltaY", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCVertexY", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoDeltaY", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoVertexY", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCTrueUpwards", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCPDG", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "AbsMCPDG", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCCosmicRay", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCFiducialLowY", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoFiducialLowY", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCIntersectsYFace", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoIntersectsYFace", -1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCLowestTenCmTotalCharge", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoLowestTenCmTotalCharge", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCLowestTenCmTotalChargeWView", -1.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoLowestTenCmTotalChargeWView", -1.f));
        return;
    }

    std::sort(clusterVector.begin(), clusterVector.end(), [](const Cluster* const pCluster1, const Cluster* const pCluster2) -> bool { return LArClusterHelper::GetLength(pCluster1) > LArClusterHelper::GetLength(pCluster2); });

    /*
    std::vector<HitWithDistance> hitWithDistanceVector, cleanHitsWithDistance;
    this->FilterHitCollection(clusterVector.front(), 5, hitWithDistanceVector);

    if (hitWithDistanceVector.size() > 10)
        cleanHitsWithDistance.insert(cleanHitsWithDistance.begin(), hitWithDistanceVector.begin(), hitWithDistanceVector.begin() + 0.72 * hitWithDistanceVector.size());

    CaloHitVector cleanHits;

    for (const auto &hitWithDistance : cleanHitsWithDistance)
        cleanHits.push_back(hitWithDistance.m_calohit);
    */

    OrderedCaloHitList orderedCaloHitList(clusterVector.front()->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::MCRelationMap mcPrimaryMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);
    LArMCParticleHelper::CaloHitToMCMap hitToMCMap;
    LArMCParticleHelper::MCContributionMap mcToTrueHitListMap; 
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcPrimaryMap, hitToMCMap, mcToTrueHitListMap);

    float mcLowestTenCmTotalCharge(0.f), recoLowestTenCmTotalCharge(0.f);
    float mcLowestTenCmTotalChargeWView(0.f), recoLowestTenCmTotalChargeWView(0.f);

    for (const auto pCaloHit : caloHitList)
    {
        const pandora::CaloHit* pParentCaloHit(static_cast<const CaloHit*>(pCaloHit->GetParentAddress()));

        if (hitToMCMap.find(pParentCaloHit) != hitToMCMap.end() && (recoLowYPosition - pCaloHit->GetPositionVector()).GetMagnitude() <= 10.0)
            mcLowestTenCmTotalCharge += pCaloHit->GetInputEnergy();

        if ((recoLowYPosition - pCaloHit->GetPositionVector()).GetMagnitude() <= 10.0)
            recoLowestTenCmTotalCharge += pCaloHit->GetInputEnergy();

        if (!(static_cast<const CaloHit*>(pCaloHit->GetParentAddress())->GetHitType() == TPC_VIEW_W))
            continue;

        //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitInformation", "HitWidth", pParentCaloHit->GetCellSize1()));
        //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitInformation", "HitInputEnergy", pParentCaloHit->GetInputEnergy()));
        //PANDORA_MONITORING_API(FillTree(this->GetPandora(), "HitInformation"));

        if (hitToMCMap.at(pParentCaloHit) == pMCParticle && (recoLowYPosition - pCaloHit->GetPositionVector()).GetMagnitude() <= 10.0)
            mcLowestTenCmTotalChargeWView += pCaloHit->GetInputEnergy();

        if ((recoLowYPosition - pCaloHit->GetPositionVector()).GetMagnitude() <= 10.0)
            recoLowestTenCmTotalChargeWView += pCaloHit->GetInputEnergy();
    } 

    /*
    for (const auto hitWithDistance : hitWithDistanceVector)
    {
        const auto pCaloHit(hitWithDistance.m_calohit);

        if (!(static_cast<const CaloHit*>(pCaloHit->GetParentAddress())->GetHitType() == TPC_VIEW_W))
            continue;

        const pandora::CaloHit* pParentCaloHit(static_cast<const CaloHit*>(pCaloHit->GetParentAddress()));

        if ((recoLowYPosition - pCaloHit->GetPositionVector()).GetMagnitude() <= 10.0 && variableNamePrefix == "LongestPfo" && pMCParticle->GetParticleId() == 13)
        {
            const MCParticleWeightMap &hitMCParticleWeightMap(pParentCaloHit->GetMCParticleWeightMap());
            float weight(0.f);
        
            if (hitMCParticleWeightMap.find(pMCParticle) != hitMCParticleWeightMap.end())
                weight = hitMCParticleWeightMap.at(pMCParticle);

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BraggHitInformation", "MCWeight", weight));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BraggHitInformation", "MCBraggPeak", (mcVertex.GetY() > mcEndpoint.GetY()) ? 1 : 0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BraggHitInformation", "MCPure", weight >= 0.9 ? 1 : 0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BraggHitInformation", "HitNNDistance", hitWithDistance.m_distance));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BraggHitInformation", "MCFiducialLowY", this->IsInFiducialVolume(mcLowYPosition) ? 1 : 0));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), "BraggHitInformation"));
        }
    } 
    */

    /*
    for (const auto pCaloHit : cleanHits)
    {
        if (!(static_cast<const CaloHit*>(pCaloHit->GetParentAddress())->GetHitType() == TPC_VIEW_W))
            continue;

        if ((recoLowYPosition - pCaloHit->GetPositionVector()).GetMagnitude() <= 10.0)
            filteredLowestTenCmTotalCharge += pCaloHit->GetInputEnergy();
    } 

    if (variableNamePrefix == "LongestPfo")
    {
        std::cout << "Longest PFO MC PDG: " << pMCParticle->GetParticleId() << std::endl;
        std::cout << "Longest PFO IsCosmic: " << LArMCParticleHelper::IsCosmicRay(pMCParticle) << std::endl;
    }
    */

    std::cout << "recoLowestTenCmTotalCharge: " << recoLowestTenCmTotalCharge << std::endl;

    //Reconstructed DeltaY and Vertex Y quantities
    CartesianVector correctedRecoLowYPosition(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(recoLowYPosition));
    CartesianVector correctedRecoHighYPosition(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(recoHighYPosition));
    const auto pVertex(LArPfoHelper::GetVertex(pPfo));
    CartesianVector vertexPosition(pVertex->GetPosition());
    CartesianVector correctedRecoVertexPosition(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(vertexPosition));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCDeltaY", std::abs(mcEndpoint.GetY() - mcVertex.GetY())));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCVertexY", mcVertex.GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoDeltaY", std::abs(correctedRecoHighYPosition.GetY() - correctedRecoLowYPosition.GetY())));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoVertexY", correctedRecoVertexPosition.GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCTrueUpwards", (mcVertex.GetY() < mcEndpoint.GetY()) ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCPDG", pMCParticle->GetParticleId()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "AbsMCPDG", std::abs(pMCParticle->GetParticleId())));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCCosmicRay", LArMCParticleHelper::IsCosmicRay(pMCParticle) ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCFiducialLowY", this->IsInFiducialVolume(mcLowYPosition) ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoFiducialLowY", this->IsInFiducialVolume(correctedRecoLowYPosition) ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCIntersectsYFace", this->MCIntersectsYFace(pMCParticle) ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoIntersectsYFace", this->RecoIntersectsYFace(pPfo) ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCLowestTenCmTotalCharge", mcLowestTenCmTotalCharge));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoLowestTenCmTotalCharge", recoLowestTenCmTotalCharge));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCLowestTenCmTotalChargeWView", mcLowestTenCmTotalChargeWView));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "RecoLowestTenCmTotalChargeWView", recoLowestTenCmTotalChargeWView));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::FilterHitCollection(const pandora::Cluster* const pCluster, int nNeighboursToConsider, std::vector<HitWithDistance> &hitWithDistanceVector) const
{
    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);

    if (caloHitList.size() < 10)
        return;

    CaloHitVector caloHitVector;

    for (const auto pCaloHit : caloHitList)
    {
        if (!((pCaloHit->GetCellSize1() > 0.3 && pCaloHit->GetCellSize1() < 1.8) && (pCaloHit->GetInputEnergy() > 50.0 && pCaloHit->GetInputEnergy() < 600.0)))
            continue;

        caloHitVector.push_back(pCaloHit);
    }

    std::sort(caloHitVector.begin(), caloHitVector.end(), [](const pandora::CaloHit* const pCaloHit1, const pandora::CaloHit* const pCaloHit2){return pCaloHit1->GetInputEnergy() < pCaloHit2->GetInputEnergy();});
    float chargeRange(caloHitVector.back()->GetInputEnergy() - caloHitVector.front()->GetInputEnergy());
    float trackLength(LArClusterHelper::GetLength(pCluster));

    for (const auto pCaloHit1 : caloHitVector)
    {    
        std::vector<float> distancesToNN;
    
        for (const auto pCaloHit2 : caloHitVector)
        {    
            if (pCaloHit1 == pCaloHit2)
                continue;

            float chargeDistance((trackLength/chargeRange) * std::abs((pCaloHit1->GetInputEnergy()) - (pCaloHit2->GetInputEnergy())));
            float lengthDistance((pCaloHit1->GetPositionVector() - pCaloHit2->GetPositionVector()).GetMagnitude());
            float distanceToNN(std::sqrt(chargeDistance*chargeDistance + lengthDistance*lengthDistance));

            distancesToNN.push_back(distanceToNN);
        }    

        std::sort(distancesToNN.begin(), distancesToNN.end());
        float nearestNeighboursDistanceSum(std::accumulate(distancesToNN.begin(), distancesToNN.begin() + nNeighboursToConsider, 0.f));
        HitWithDistance hitWithDistance(pCaloHit1, nearestNeighboursDistanceSum);
        hitWithDistanceVector.push_back(hitWithDistance);
    }    

    std::sort(hitWithDistanceVector.begin(), hitWithDistanceVector.end(), [](HitWithDistance &hitWithDistance1, HitWithDistance &hitWithDistance2) {return hitWithDistance1.m_distance < hitWithDistance2.m_distance;} );
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::MCIntersectsYFace(const pandora::MCParticle* pMCParticle) const
{
    const pandora::CartesianVector initialPosition(pMCParticle->GetVertex());
    const pandora::CartesianVector endPosition(pMCParticle->GetEndpoint());

    const pandora::CartesianVector lowYVector(initialPosition.GetY() < endPosition.GetY() ? initialPosition : endPosition);
    const pandora::CartesianVector highYVector(initialPosition.GetY() > endPosition.GetY() ? initialPosition : endPosition);

    float xExtent(highYVector.GetX() - lowYVector.GetX()), yExtent(highYVector.GetY() - lowYVector.GetY()), zExtent(highYVector.GetZ() - lowYVector.GetZ());
    float xSlope(xExtent/yExtent), zSlope(zExtent/yExtent);
    float yDistanceToTravel(116.5 - lowYVector.GetY());
    CartesianVector yFaceIntersection(lowYVector.GetX() + xSlope*yDistanceToTravel, 116.5, lowYVector.GetZ() + zSlope*yDistanceToTravel);

    if (yFaceIntersection.GetX() > 0.0 && yFaceIntersection.GetX() < 256.35 && yFaceIntersection.GetZ() > 0.0 && yFaceIntersection.GetZ() < 1036.8)
        return true;
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::RecoIntersectsYFace(const ParticleFlowObject *const pPfo) const
{
    LArTrackStateVector trackStateVector;
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, LArPfoHelper::GetVertex(pPfo), m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()), trackStateVector);

    std::sort(trackStateVector.begin(), trackStateVector.end(), [](TrackState pTrackState1, TrackState pTrackState2) -> bool { return pTrackState1.GetPosition().GetY() < pTrackState2.GetPosition().GetY(); });

    const pandora::CartesianVector initialPosition(trackStateVector.front().GetPosition());
    const pandora::CartesianVector endPosition(trackStateVector.back().GetPosition());

    pandora::CartesianVector lowYVector(initialPosition.GetY() < endPosition.GetY() ? initialPosition : endPosition);
    pandora::CartesianVector highYVector(initialPosition.GetY() > endPosition.GetY() ? initialPosition : endPosition);
    
    const pandora::CartesianVector correctedLowYVector(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(lowYVector));
    const pandora::CartesianVector correctedHighYVector(LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(highYVector));

    float xExtent(correctedHighYVector.GetX() - correctedLowYVector.GetX()), yExtent(correctedHighYVector.GetY() - correctedLowYVector.GetY()), zExtent(correctedHighYVector.GetZ() - correctedLowYVector.GetZ());
    float xSlope(xExtent/yExtent), zSlope(zExtent/yExtent);
    float yDistanceToTravel(116.5 - correctedLowYVector.GetY());
    CartesianVector yFaceIntersection(correctedLowYVector.GetX() + xSlope*yDistanceToTravel, 116.5, correctedLowYVector.GetZ() + zSlope*yDistanceToTravel);

    if (yFaceIntersection.GetX() > 0.0 && yFaceIntersection.GetX() < 256.35 && yFaceIntersection.GetZ() > 0.0 && yFaceIntersection.GetZ() < 1036.8)
        return true;
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::WritePIDVariables(const ParticleFlowObject *const pPfo, std::string variableNamePrefix) const
{
    pandora::CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, caloHitList);

    pandora::CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
    std::sort(caloHitVector.begin(), caloHitVector.end(), [](const CaloHit* const pCaloHit1, const CaloHit* const pCaloHit2) -> bool { return pCaloHit1->GetPositionVector().GetZ() < pCaloHit2->GetPositionVector().GetZ(); }); 

    float meanHitCharge(0.f), minHitCharge(1e6), maxHitCharge(0.f);

    for (const auto pCaloHit : caloHitVector)
    {
        meanHitCharge += pCaloHit->GetInputEnergy();

        if (pCaloHit->GetInputEnergy() < minHitCharge)
            minHitCharge = pCaloHit->GetInputEnergy();

        if (pCaloHit->GetInputEnergy() > maxHitCharge)
            maxHitCharge = pCaloHit->GetInputEnergy();
    }

    float pfoLength(sqrt(LArPfoHelper::GetThreeDLengthSquared(pPfo)));
    float hitChargeRangeOverLength(pfoLength > 1.0 ? (maxHitCharge - minHitCharge)/pfoLength : 0.f);

    float leftSumCharge(-1.f), rightSumCharge(-1.f);

    if (caloHitVector.size() >= 3)
    {
        leftSumCharge = (caloHitVector.at(0)->GetInputEnergy() + caloHitVector.at(1)->GetInputEnergy() + caloHitVector.at(2)->GetInputEnergy());
        rightSumCharge = (caloHitVector.at(caloHitVector.size() - 3)->GetInputEnergy() + caloHitVector.at(caloHitVector.size() - 2)->GetInputEnergy() + caloHitVector.at(caloHitVector.size() - 1)->GetInputEnergy());
    }

    float minSumCharge(std::min(leftSumCharge, rightSumCharge)), maxSumCharge(std::max(leftSumCharge, rightSumCharge));
    float sumChargeRatio(maxSumCharge/minSumCharge);

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "HitChargeRangeOverLength", hitChargeRangeOverLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MinSumCharge", minSumCharge));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MaxSumCharge", maxSumCharge));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "SumChargeRatio", sumChargeRatio));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MeanHitCharge", meanHitCharge));

    try
    {
        DirectionFittingThreeDTool::DirectionFitObject fitResult3D = m_pDirectionFittingThreeDTool->GetPfoDirection(pPfo);
        //fitResult3D.DrawFit();

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "FitSuccesful3D", 1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "FitMass3D", fitResult3D.GetFitParameters().GetParameterOne()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCPDG3D", fitResult3D.GetMCParent()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "Contained3D", fitResult3D.GetContained() ? 1 : 0));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MinChiSquaredPerHit3D", fitResult3D.GetMinChiSquaredPerHit()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "NumberHits3D", fitResult3D.GetNHits()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "FitStatus3D", fitResult3D.GetFitStatus()));
    }
    catch (...)
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "FitSuccesful3D", 0));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "FitMass3D", 0.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MCPDG3D", 0));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "Contained3D", 0));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "MinChiSquaredPerHit3D", 0.f));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "NumberHits3D", 0));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventSelection", variableNamePrefix + "FitStatus3D", -1));
    }

    
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EventValidationAlgorithm::CalculateCosmicProbability(TrackDirectionTool::DirectionFitObject &directionFit) const
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

    if (!this->IntersectsYFace(directionFit))
        probability = 0.0;

    return probability;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::IntersectsYFace(TrackDirectionTool::DirectionFitObject &fitResult) const
{
    const pandora::CartesianVector initialPosition(fitResult.GetBeginpoint());
    const pandora::CartesianVector endPosition(fitResult.GetEndpoint());

    const pandora::CartesianVector lowYVector(initialPosition.GetY() < endPosition.GetY() ? initialPosition : endPosition);
    const pandora::CartesianVector highYVector(initialPosition.GetY() > endPosition.GetY() ? initialPosition : endPosition);

    float xExtent(highYVector.GetX() - lowYVector.GetX()), yExtent(highYVector.GetY() - lowYVector.GetY()), zExtent(highYVector.GetZ() - lowYVector.GetZ());
    float xSlope(xExtent/yExtent), zSlope(zExtent/yExtent);
    float yDistanceToTravel(116.5 - lowYVector.GetY());
    CartesianVector yFaceIntersection(lowYVector.GetX() + xSlope*yDistanceToTravel, 116.5, lowYVector.GetZ() + zSlope*yDistanceToTravel);

    if (yFaceIntersection.GetX() > 0.0&& yFaceIntersection.GetX() < 256.35 && yFaceIntersection.GetZ() > 0.0 && yFaceIntersection.GetZ() < 1036.8)
        return true;
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EventValidationAlgorithm::CalculateTrackProbability(const ParticleFlowObject *const pPfo) const
{
    const PropertiesMap &properties(pPfo->GetPropertiesMap());
    const PropertiesMap::const_iterator iter(properties.find("TrackScore"));

    if (iter != properties.end())
        return iter->second; 
    else
        return -1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int EventValidationAlgorithm::GetInteractionType() const
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

    if (mcPrimaryList.size() == 0)
        return -1;

    const LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(mcPrimaryList));
    //return LArInteractionTypeHelper::ToString(interactionType);
    return static_cast<int>(static_cast<unsigned int>(interactionType));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetNumberTracksAndShowers(pandora::PfoList pfoList, int &nTracks, int &nShowers) const
{
    for (const auto pPfo : pfoList)
    {
        if (LArPfoHelper::IsTrack(pPfo))
            ++nTracks;

        if (LArPfoHelper::IsShower(pPfo))
            ++nShowers;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EventValidationAlgorithm::GetTotalEventCharge(const pandora::CaloHitList *const pCaloHitList) const
{
    float totalEventCharge(0.f);
    
    for (const auto pCaloHit : *pCaloHitList)
        totalEventCharge += pCaloHit->GetInputEnergy();

    return totalEventCharge; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::ParticleFlowObject* EventValidationAlgorithm::GetLongestPfo(pandora::PfoList pfoList) const
{
    float longestPfoLength(0.f);
    const pandora::ParticleFlowObject* pLongestPfo = NULL;

    for (const auto pPfo : pfoList)
    {
        if (!LArPfoHelper::IsThreeD(pPfo))
            continue;

        if (LArPfoHelper::GetThreeDLengthSquared(pPfo) > longestPfoLength)
        {
            pLongestPfo = pPfo;
            longestPfoLength = LArPfoHelper::GetThreeDLengthSquared(pPfo);
        }
    }

    if (pLongestPfo == NULL)
    {
        throw STATUS_CODE_NOT_FOUND;
    }

    return pLongestPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::ParticleFlowObject* EventValidationAlgorithm::GetShortestPfo(pandora::PfoList pfoList) const
{
    float shortestPfoLength(1e6);
    const pandora::ParticleFlowObject* pShortestPfo = NULL;

    for (const auto pPfo : pfoList)
    {
        if (!LArPfoHelper::IsThreeD(pPfo))
            continue;

        if (LArPfoHelper::GetThreeDLengthSquared(pPfo) < shortestPfoLength)
        {
            pShortestPfo = pPfo;
            shortestPfoLength = LArPfoHelper::GetThreeDLengthSquared(pPfo);
        }
    }

    if (pShortestPfo == NULL)
    {
        throw STATUS_CODE_NOT_FOUND;
    }

    return pShortestPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EventValidationAlgorithm::GetPfoCharge(const pandora::ParticleFlowObject* pPfo) const
{
    float pfoTotalCharge(0.f);

    pandora::PfoList pfoList;
    pfoList.insert(pfoList.begin(), pPfo);

    pandora::CaloHitList pfoHitList;
    this->GetClusterHits(pfoList, pfoHitList);

    for (const auto pCaloHit : pfoHitList)
        pfoTotalCharge += pCaloHit->GetInputEnergy();

    return pfoTotalCharge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int EventValidationAlgorithm::CountPfoHits(const pandora::ParticleFlowObject* pPfo) const
{
    pandora::PfoList pfoList;
    pfoList.insert(pfoList.begin(), pPfo);

    pandora::CaloHitList pfoHitList;
    this->GetClusterHits(pfoList, pfoHitList);

    return pfoHitList.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EventValidationAlgorithm::GetThetaBeamPfo(const pandora::ParticleFlowObject* pPfo) const
{
    LArTrackStateVector trackStateVector; 
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, LArPfoHelper::GetVertex(pPfo), m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()), trackStateVector);

    if (trackStateVector.size() == 0)
        return -1.f;

    const pandora::CartesianVector zAxis(0.f, 0.f, 1.f);
    const pandora::CartesianVector vertexDirection(trackStateVector.front().GetDirection());
    return vertexDirection.GetOpeningAngle(zAxis); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EventValidationAlgorithm::GetAngleXAxis(const pandora::ParticleFlowObject* pPfo) const
{
    LArTrackStateVector trackStateVector; 
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, LArPfoHelper::GetVertex(pPfo), m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()), trackStateVector);

    if (trackStateVector.size() == 0)
        return -1.f;

    const pandora::CartesianVector xAxis(1.f, 0.f, 0.f);
    const pandora::CartesianVector vertexDirection(trackStateVector.front().GetDirection());
    return vertexDirection.GetOpeningAngle(xAxis); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EventValidationAlgorithm::GetAngleWithVector(const pandora::ParticleFlowObject* pPfo, CartesianVector &axisVector) const
{
    LArTrackStateVector trackStateVector; 
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, LArPfoHelper::GetVertex(pPfo), m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()), trackStateVector);

    if (trackStateVector.size() == 0)
        return -1.f;

    const pandora::CartesianVector vertexDirection(trackStateVector.front().GetDirection());
    return vertexDirection.GetOpeningAngle(axisVector); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EventValidationAlgorithm::GetPfoOpeningAngle(const pandora::ParticleFlowObject* pPfo1, const pandora::ParticleFlowObject* pPfo2) const
{
    try
    {
        LArTrackStateVector trackStateVector1; 
        LArPfoHelper::GetSlidingFitTrajectory(pPfo1, LArPfoHelper::GetVertex(pPfo1), m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()), trackStateVector1);

        LArTrackStateVector trackStateVector2; 
        LArPfoHelper::GetSlidingFitTrajectory(pPfo2, LArPfoHelper::GetVertex(pPfo2), m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()), trackStateVector2);

        if (trackStateVector1.size() == 0 || trackStateVector2.size() == 0)
            return 0.f;

        const pandora::CartesianVector vertexDirection1(trackStateVector1.front().GetDirection());
        const pandora::CartesianVector vertexDirection2(trackStateVector2.front().GetDirection());

        return vertexDirection1.GetOpeningAngle(vertexDirection2); 
    }
    catch (...)
    {
        return -1.f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector EventValidationAlgorithm::GetApproximateNeutrinoMomentum(pandora::PfoList pfoList, const pandora::ParticleFlowObject* pLongestPfo) const
{
    const float muonMass(105.7), protonMass(938.3);
    CartesianVector muonMomentum(GetApproximatePfoMomentum(pLongestPfo, muonMass)); 

    pandora::CartesianVector neutrinoMomentum = muonMomentum;

    for (const auto pPfo : pfoList)
    {
        if (pPfo == pLongestPfo)
            continue;

        try
        {
            pandora::CartesianVector protonMomentum(GetApproximatePfoMomentum(pPfo, protonMass)); 
            neutrinoMomentum += protonMomentum; 
        }
        catch (...)
        {
            continue;
        }
    }
   
    return neutrinoMomentum.GetUnitVector(); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector EventValidationAlgorithm::GetApproximatePfoMomentum(const pandora::ParticleFlowObject* pPfo, const float &particleMass) const
{
    LArTrackStateVector trackStateVector; 
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, LArPfoHelper::GetVertex(pPfo), m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()), trackStateVector);

    if (trackStateVector.size() == 0)
    {
        CartesianVector dummyVector(1.f, 1.f, 1.f);
        return dummyVector;
    }

    const pandora::CartesianVector particleDirection(trackStateVector.front().GetDirection().GetUnitVector());
    float particleCharge(GetPfoCharge(pPfo));

    const float adcToMeV(0.0000236 * 197 * (1.0/0.62));
    const float momentumNorm(std::sqrt(adcToMeV * adcToMeV * particleCharge * particleCharge + 2 * adcToMeV * particleCharge * particleMass));

    pandora::CartesianVector particleMomentum(particleDirection.GetX() * momentumNorm, particleDirection.GetY() * momentumNorm, particleDirection.GetZ() * momentumNorm);
    return particleMomentum;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    LArSpaceChargeHelper::Configure("/usera/jjd49/pandora_direction/PandoraPFA/LArContent-origin/vertex_direction/larpandoracontent/LArDirection/SCEoffsets_MicroBooNE_E273.root");

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseTrueNeutrinosOnly", m_useTrueNeutrinosOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TestBeamMode", m_testBeamMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectInputHits", m_selectInputHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitSharingFraction", m_minHitSharingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPhotonPropagation", m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintAllToScreen", m_printAllToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintMatchingToScreen", m_printMatchingToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseSmallPrimaries", m_useSmallPrimaries));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinSharedHits", m_matchingMinSharedHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinCompleteness", m_matchingMinCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinPurity", m_matchingMinPurity));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ViewEvent", m_viewEvent));

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

    //3D direction tool
    AlgorithmTool *pDirectionFittingThreeDTool(nullptr);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "TrackDirection3D", pDirectionFittingThreeDTool));

    if (!(this->m_pDirectionFittingThreeDTool = dynamic_cast<DirectionFittingThreeDTool*>(pDirectionFittingThreeDTool)))
        throw STATUS_CODE_FAILURE;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
