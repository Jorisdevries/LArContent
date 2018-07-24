/**
 *  @file   LArDirection/EventSelectionAlgorithm.cc
 * 
 *  @brief  Implementation of an event selection algorithm 
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "EventSelectionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

EventSelectionAlgorithm::EventSelectionAlgorithm() :
    m_writeToTree(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventSelectionAlgorithm::~EventSelectionAlgorithm()
{
    if (m_writeToTree)
    {
        //PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
        //PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "TrueNeutrinos", m_fileName.c_str(), "UPDATE"));
        //PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ReconstructedNeutrinos", m_fileName.c_str(), "UPDATE"));
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "EventDescription", m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventSelectionAlgorithm::Run()
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    this->WriteEventDescription(pMCParticleList, pPfoList, pCaloHitList);

    //this->WriteVariables(pPfoList, pCaloHitList);
    //this->WriteNeutrinoInformation(pMCParticleList, pPfoList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSelectionAlgorithm::WriteEventDescription(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList, const pandora::CaloHitList* pCaloHitList) const
{
    LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, LArMCParticleHelper::PrimaryParameters(),
        LArMCParticleHelper::IsBeamNeutrinoFinalState, nuMCParticlesToGoodHitsMap);

    int trueInteractionType(this->GetInteractionType(nuMCParticlesToGoodHitsMap));

    pandora::MCParticleVector trueNeutrinos;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);

    if (trueNeutrinos.size() != 1) return;

    int neutrinoNuanceCode(LArMCParticleHelper::GetNuanceCode(trueNeutrinos.front()));

    pandora::PfoList trueNeutrinoPfos(this->GetTrueNeutrinoAssociatedPfos(pPfoList, nuMCParticlesToGoodHitsMap));  
    int nTrueNeutrinoAssociatedTracks(this->CountTrackPfos(trueNeutrinoPfos));

    PfoList neutrinoPfos;
    LArPfoHelper::GetRecoNeutrinos(pPfoList, neutrinoPfos);

    if (neutrinoPfos.size() != 1) return;

    pandora::PfoList recoNeutrinoPrimaryDaughters(this->GetPrimaryDaughters(neutrinoPfos, nuMCParticlesToGoodHitsMap));
    int nRecoNeutrinoAssociatedTracks(this->CountTrackPfos(recoNeutrinoPrimaryDaughters));

    int interactionTypeResemblance(this->GetRecoNeutrinoInteractionTypeResemblance(recoNeutrinoPrimaryDaughters, nuMCParticlesToGoodHitsMap));
    int chosenSliceContainsNeutrino(trueInteractionType == interactionTypeResemblance ? 1 : 0);
    int nSliceCosmicRays(this->GetNumberCosmicRaysChosenSlice(neutrinoPfos));

    pandora::PfoList allConnectedPfos;
    LArPfoHelper::GetAllConnectedPfos(neutrinoPfos, allConnectedPfos);

    pandora::CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(allConnectedPfos, TPC_3D, caloHitList);
    int chosenSliceNumberHits(caloHitList.size());

    std::cout << "chosenSliceNumberHits: " << chosenSliceNumberHits << std::endl;
    std::cout << "nSliceCosmicRays: " << nSliceCosmicRays << std::endl;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventDescription", "TrueInteractionType", trueInteractionType));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventDescription", "NeutrinoNuanceCode", neutrinoNuanceCode));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventDescription", "TrueNeutrinoNumberAssociatedParticles", (int)trueNeutrinoPfos.size()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventDescription", "TrueNeutrinoNumberAssociatedTracks", nTrueNeutrinoAssociatedTracks));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventDescription", "RecoNeutrinoNumberAssociatedParticles", (int)recoNeutrinoPrimaryDaughters.size()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventDescription", "RecoNeutrinoNumberAssociatedTracks", nRecoNeutrinoAssociatedTracks));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventDescription", "ChosenSliceInteractionType", interactionTypeResemblance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventDescription", "ChosenSliceContainsTrueNeutrino", chosenSliceContainsNeutrino));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventDescription", "ChosenSliceNumberCosmicRays", nSliceCosmicRays));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventDescription", "ChosenSliceNumberHits", chosenSliceNumberHits));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "EventDescription"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

int EventSelectionAlgorithm::GetInteractionType(LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap) const
{
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

pandora::PfoList EventSelectionAlgorithm::GetTrueNeutrinoAssociatedPfos(const pandora::PfoList* pPfoList, LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap) const
{
    PfoList neutrinoPfos, trueNeutrinoPfos;
    LArPfoHelper::GetRecoNeutrinos(pPfoList, neutrinoPfos);

    for (const auto pPfo : *pPfoList)
    {
        if (pPfo == neutrinoPfos.front() || !LArPfoHelper::IsFinalState(pPfo))
            continue;

        const pandora::MCParticle* pMCParticle(this->GetMainMCParticle(pPfo));
        const pandora::MCParticle* pPrimaryMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);
        const pandora::MCParticle* pParentMCParticle = LArMCParticleHelper::GetParentMCParticle(pMCParticle);

        bool mcParticleReconstructable = (nuMCParticlesToGoodHitsMap.find(pPrimaryMCParticle) != nuMCParticlesToGoodHitsMap.end()); 

        if (mcParticleReconstructable && LArMCParticleHelper::IsNeutrino(pParentMCParticle))
            trueNeutrinoPfos.insert(trueNeutrinoPfos.begin(), pPfo);
    }

    return trueNeutrinoPfos;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int EventSelectionAlgorithm::CountTrackPfos(pandora::PfoList &pfoList) const
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

pandora::PfoList EventSelectionAlgorithm::GetPrimaryDaughters(pandora::PfoList &neutrinoPfos, LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap) const
{
    pandora::PfoList allConnectedPfos, primaryDaughters;
    LArPfoHelper::GetAllConnectedPfos(neutrinoPfos, allConnectedPfos);

    for (const auto pPfo : allConnectedPfos)
    {
        if (pPfo == neutrinoPfos.front() || !LArPfoHelper::IsFinalState(pPfo))
            continue;

        const pandora::MCParticle* pMCParticle(this->GetMainMCParticle(pPfo));
        const pandora::MCParticle* pPrimaryMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);

        bool mcParticleReconstructable = (nuMCParticlesToGoodHitsMap.find(pPrimaryMCParticle) != nuMCParticlesToGoodHitsMap.end()); 
        
        if (mcParticleReconstructable)
            primaryDaughters.insert(primaryDaughters.begin(), pPfo);
    }

    return primaryDaughters;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int EventSelectionAlgorithm::GetRecoNeutrinoInteractionTypeResemblance(pandora::PfoList &primaryDaughters, LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap) const
{
    //Takes the PFOs in the reconstructed neutrino, their main MC particles, select the MC particles that are neutrino-induced, and construct an interaction type
    pandora::MCParticleList mcPrimaryList;

    for (const auto pPfo : primaryDaughters)
    {
        const pandora::MCParticle* pMCParticle(this->GetMainMCParticle(pPfo));
        const pandora::MCParticle* pPrimaryMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);
        const pandora::MCParticle* pParentMCParticle = LArMCParticleHelper::GetParentMCParticle(pMCParticle);

        bool mcParticleInPrimariesList = (std::find(mcPrimaryList.begin(), mcPrimaryList.end(), pPrimaryMCParticle) != mcPrimaryList.end());
        bool mcParticleReconstructable = (nuMCParticlesToGoodHitsMap.find(pPrimaryMCParticle) != nuMCParticlesToGoodHitsMap.end()); 

        if (LArMCParticleHelper::IsNeutrino(pParentMCParticle) && !mcParticleInPrimariesList && mcParticleReconstructable)
            mcPrimaryList.insert(mcPrimaryList.begin(), pPrimaryMCParticle);   
    }

    if (mcPrimaryList.size() == 0)
        return -1;

    LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(mcPrimaryList));
    //std::string interactionTypeString(LArInteractionTypeHelper::ToString(interactionType));
    return static_cast<int>(static_cast<unsigned int>(interactionType));
}

//------------------------------------------------------------------------------------------------------------------------------------------

int EventSelectionAlgorithm::GetNumberCosmicRaysChosenSlice(pandora::PfoList &neutrinoPfos) const
{
    int nCosmicRays(0);

    pandora::PfoList allConnectedPfos;
    LArPfoHelper::GetAllConnectedPfos(neutrinoPfos, allConnectedPfos);

    for (const auto pPfo : allConnectedPfos)
    {
        if (pPfo == neutrinoPfos.front())
            continue;

        if (!pPfo->GetParentPfoList().empty() && pPfo->GetParentPfoList().front() != neutrinoPfos.front())
            continue;

        const pandora::MCParticle* pMCParticle(this->GetMainMCParticle(pPfo));
        const pandora::MCParticle* pParentMCParticle = LArMCParticleHelper::GetParentMCParticle(pMCParticle);

        if (LArMCParticleHelper::IsCosmicRay(pParentMCParticle))
            ++nCosmicRays; 
    }

    return nCosmicRays;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void EventSelectionAlgorithm::WriteNeutrinoInformation(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList) const
{
    this->WriteTrueNeutrinoInformation(pMCParticleList, pPfoList);
    this->WriteReconstructedNeutrinoInformation(pMCParticleList, pPfoList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSelectionAlgorithm::WriteTrueNeutrinoInformation(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList) const
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, LArMCParticleHelper::PrimaryParameters(), LArMCParticleHelper::IsBeamNeutrinoFinalState, nuMCParticlesToGoodHitsMap);

    pandora::MCParticleVector trueNeutrinos;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);

    //Only consider events with one neutrino
    if (trueNeutrinos.size() != 1)
        return;

    int nuanceCode(LArMCParticleHelper::GetNuanceCode(trueNeutrinos.front()));

    PfoList neutrinoPfos;
    LArPfoHelper::GetRecoNeutrinos(pPfoList, neutrinoPfos);

    int nDaughterParticles(0), nDaughterTracks(0);
    int nReconstructableDaughterParticles(0), nReconstructableDaughterTracks(0);

    for (const auto pPfo : *pPfoList)
    {
        if (pPfo == neutrinoPfos.front() || !LArPfoHelper::IsFinalState(pPfo))
            continue;

        const pandora::MCParticle* pMCParticle(this->GetMainMCParticle(pPfo));
        const pandora::MCParticle* pParentMCParticle = LArMCParticleHelper::GetParentMCParticle(pMCParticle);

        if (LArMCParticleHelper::IsNeutrino(pParentMCParticle))
        {
            ++nDaughterParticles;
        
            if (LArPfoHelper::IsTrack(pPfo))
                ++nDaughterTracks;

            bool mcParticleReconstructable = (nuMCParticlesToGoodHitsMap.find(pMCParticle) != nuMCParticlesToGoodHitsMap.end()); 

            if (mcParticleReconstructable)
            {
                ++nReconstructableDaughterParticles;
            
                if (LArPfoHelper::IsTrack(pPfo))
                    ++nReconstructableDaughterTracks;
            }
        }
    }


    std::cout << "nuanceCode: " << nuanceCode << std::endl;
    std::cout << "True neutrino nDaughterParticles: " << nDaughterParticles << std::endl;
    std::cout << "True neutrino nDaughterTracks: " << nDaughterTracks << std::endl;
    std::cout << "True neutrino nReconstructableDaughterParticles: " << nReconstructableDaughterParticles << std::endl;
    std::cout << "True neutrino nReconstructableDaughterTracks: " << nReconstructableDaughterTracks << std::endl;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrueNeutrinos", "NuanceCode", nuanceCode));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrueNeutrinos", "nDaughterParticles", nDaughterParticles));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrueNeutrinos", "nDaughterTracks", nDaughterTracks));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrueNeutrinos", "nReconstructableDaughterParticles", nReconstructableDaughterParticles));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrueNeutrinos", "nReconstructableDaughterTracks", nReconstructableDaughterTracks));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "TrueNeutrinos"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSelectionAlgorithm::WriteReconstructedNeutrinoInformation(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList) const
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, LArMCParticleHelper::PrimaryParameters(), LArMCParticleHelper::IsBeamNeutrinoFinalState, nuMCParticlesToGoodHitsMap);

    PfoList neutrinoPfos, allConnectedPfos;
    LArPfoHelper::GetRecoNeutrinos(pPfoList, neutrinoPfos);
    LArPfoHelper::GetAllConnectedPfos(neutrinoPfos, allConnectedPfos);

    if (neutrinoPfos.size() != 1)
        return;

    pandora::MCParticleList mcPrimaryList;

    int nDaughterParticles(0), nDaughterTracks(0);
    int nReconstructableDaughterParticles(0), nReconstructableDaughterTracks(0);
    int nHits(0), nReconstructableHits(0);

    for (const auto pPfo : allConnectedPfos)
    {
        if (pPfo == neutrinoPfos.front() || !LArPfoHelper::IsFinalState(pPfo))
            continue;

        ++nDaughterParticles;
    
        if (LArPfoHelper::IsTrack(pPfo))
            ++nDaughterTracks;

        const pandora::MCParticle* pMCParticle(this->GetMainMCParticle(pPfo));
        const pandora::MCParticle* pPrimaryMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);

        bool mcParticleInPrimariesList = (std::find(mcPrimaryList.begin(), mcPrimaryList.end(), pPrimaryMCParticle) != mcPrimaryList.end());
        bool mcParticleReconstructable = (nuMCParticlesToGoodHitsMap.find(pPrimaryMCParticle) != nuMCParticlesToGoodHitsMap.end()); 

        pandora::CaloHitList pfoHits;
        LArPfoHelper::GetCaloHits(pPfo, TPC_3D, pfoHits);
        nHits += pfoHits.size();

        if (!mcParticleInPrimariesList && mcParticleReconstructable)
        {
            mcPrimaryList.insert(mcPrimaryList.begin(), pPrimaryMCParticle);   
            ++nReconstructableDaughterParticles;

            pandora::CaloHitList reconstructablePfoHits;
            LArPfoHelper::GetCaloHits(pPfo, TPC_3D, reconstructablePfoHits);
            nReconstructableHits += reconstructablePfoHits.size();

            if (LArPfoHelper::IsTrack(pPfo))
                ++nReconstructableDaughterTracks;
        }
    }

    LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(mcPrimaryList));
    std::string interactionTypeString(LArInteractionTypeHelper::ToString(interactionType));

    std::cout << "interactionTypeString: " << interactionTypeString << std::endl;
    std::cout << "Reco neutrino nDaughterParticles: " << nDaughterParticles << std::endl;
    std::cout << "Reco neutrino nDaughterTracks: " << nDaughterTracks << std::endl;
    std::cout << "Reco neutrino nReconstructableDaughterParticles: " << nReconstructableDaughterParticles << std::endl;
    std::cout << "Reco neutrino nReconstructableDaughterTracks: " << nReconstructableDaughterTracks << std::endl;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ReconstructedNeutrinos", "InteractionType", static_cast<int>(static_cast<unsigned int>(interactionType))));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ReconstructedNeutrinos", "nDaughterParticles", nDaughterParticles));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ReconstructedNeutrinos", "nDaughterTracks", nDaughterTracks));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ReconstructedNeutrinos", "nReconstructableDaughterParticles", nReconstructableDaughterParticles));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ReconstructedNeutrinos", "nReconstructableDaughterTracks", nReconstructableDaughterTracks));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ReconstructedNeutrinos", "nHits", nHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ReconstructedNeutrinos", "nReconstructableHits", nReconstructableHits));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ReconstructedNeutrinos"));

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* EventSelectionAlgorithm::GetMainMCParticle(const ParticleFlowObject *const pPfo)
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

void EventSelectionAlgorithm::WriteVariables(const pandora::PfoList* pPfoList, const pandora::CaloHitList *pCaloHitList) const
{
    PfoList allConnectedPfos;
    LArPfoHelper::GetAllConnectedPfos(*pPfoList, allConnectedPfos);

    //Get interaction type
    int interactionType(GetInteractionType());
    std::cout << "Interaction type: " << interactionType << std::endl;

    //Get number of tracks and showers
    int nTracks(0), nShowers(0);
    GetNumberTracksAndShowers(allConnectedPfos, nTracks, nShowers);

    //Get total event charge
    const float totalEventCharge(GetTotalEventCharge(pCaloHitList));

    //Get shortest and longest PFO
    const pandora::ParticleFlowObject* pLongestPfo(GetLongestPfo(allConnectedPfos));
    const pandora::ParticleFlowObject* pShortestPfo(GetShortestPfo(allConnectedPfos));

    //Longest and shortest pfo total charges
    const float shortestPfoCharge(GetPfoCharge(pShortestPfo));
    const float longestPfoCharge(GetPfoCharge(pLongestPfo));

    //longest pfo opening angle with beamline
    const float thetaBeamPfo(GetThetaBeamPfo(pLongestPfo));

    //longest and shortest pfo opening angle
    const float phiPfoOpeningAngle(GetPfoOpeningAngle(pLongestPfo, pShortestPfo));

    //longest and shortest pfo lengths 
    const float shortestPfoLength(LArPfoHelper::GetThreeDLengthSquared(pShortestPfo)), longestPfoLength(LArPfoHelper::GetThreeDLengthSquared(pLongestPfo));

    //longitudinal and transverse momenta
    pandora::CartesianVector neutrinoMomentum(GetApproximateNeutrinoMomentum(pPfoList, pLongestPfo));

    //dE/dx fit variables
    TrackDirectionTool::DirectionFitObject fitResult = m_pTrackDirectionTool->GetPfoDirection(pLongestPfo);
    TrackDirectionTool::FitParameters fitParameters(fitResult.GetFitParameters());

    //Write all information to tree
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "InteractionType", interactionType));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "IsSingleMuon", (interactionType == 0 ? 1 : 0)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "IsMuonProton", (interactionType == 1 ? 1 : 0)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NumberTracks", nTracks));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NumberShowers", nShowers));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TotalEventCharge", totalEventCharge));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "LongestPfoCharge", longestPfoCharge));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ShortestPfoCharge", shortestPfoCharge));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "LongestPfoLength", longestPfoLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ShortestPfoLength", shortestPfoLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Theta", thetaBeamPfo));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Phi", phiPfoOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NeutrinoMomentumX", neutrinoMomentum.GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NeutrinoMomentumY", neutrinoMomentum.GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NeutrinoMomentumZ", neutrinoMomentum.GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "FitParameterZero", fitParameters.GetParameterZero()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "FitParameterOne", fitParameters.GetParameterOne()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "FitParameterTwo", fitParameters.GetParameterTwo()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "FitParameterThree", fitParameters.GetParameterThree()));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));

}

//------------------------------------------------------------------------------------------------------------------------------------------

int EventSelectionAlgorithm::GetInteractionType() const
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

void EventSelectionAlgorithm::GetNumberTracksAndShowers(pandora::PfoList pfoList, int &nTracks, int &nShowers) const
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

float EventSelectionAlgorithm::GetTotalEventCharge(const pandora::CaloHitList *const pCaloHitList) const
{
    float totalEventCharge(0.f);
    
    for (const auto pCaloHit : *pCaloHitList)
        totalEventCharge += pCaloHit->GetInputEnergy();

    return totalEventCharge; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::ParticleFlowObject* EventSelectionAlgorithm::GetLongestPfo(pandora::PfoList pfoList) const
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

const pandora::ParticleFlowObject* EventSelectionAlgorithm::GetShortestPfo(pandora::PfoList pfoList) const
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

float EventSelectionAlgorithm::GetPfoCharge(const pandora::ParticleFlowObject* pPfo) const
{
    float pfoTotalCharge(0.f);

    pandora::CaloHitList pfoHitList;
    pandora::HitType hitType(TPC_VIEW_W);
    LArPfoHelper::GetCaloHits(pPfo, hitType, pfoHitList);

    for (const auto pCaloHit : pfoHitList)
        pfoTotalCharge += pCaloHit->GetInputEnergy();

    return pfoTotalCharge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EventSelectionAlgorithm::GetThetaBeamPfo(const pandora::ParticleFlowObject* pPfo) const
{
    LArTrackStateVector trackStateVector; 
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, LArPfoHelper::GetVertex(pPfo), 20, LArGeometryHelper::GetWireZPitch(this->GetPandora()), trackStateVector);

    const pandora::CartesianVector zAxis(0.f, 0.f, 1.f);
    const pandora::CartesianVector vertexDirection(trackStateVector.front().GetDirection());
    return vertexDirection.GetOpeningAngle(zAxis); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EventSelectionAlgorithm::GetPfoOpeningAngle(const pandora::ParticleFlowObject* pPfo1, const pandora::ParticleFlowObject* pPfo2) const
{
    LArTrackStateVector trackStateVector1; 
    LArPfoHelper::GetSlidingFitTrajectory(pPfo1, LArPfoHelper::GetVertex(pPfo1), 20, LArGeometryHelper::GetWireZPitch(this->GetPandora()), trackStateVector1);

    LArTrackStateVector trackStateVector2; 
    LArPfoHelper::GetSlidingFitTrajectory(pPfo2, LArPfoHelper::GetVertex(pPfo2), 20, LArGeometryHelper::GetWireZPitch(this->GetPandora()), trackStateVector2);

    if (trackStateVector1.size() == 0 || trackStateVector2.size() == 0)
        return 0.f;

    const pandora::CartesianVector vertexDirection1(trackStateVector1.front().GetDirection());
    const pandora::CartesianVector vertexDirection2(trackStateVector2.front().GetDirection());

    return vertexDirection1.GetOpeningAngle(vertexDirection2); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector EventSelectionAlgorithm::GetApproximateNeutrinoMomentum(const pandora::PfoList* pPfoList, const pandora::ParticleFlowObject* pLongestPfo) const
{
    const float muonMass(105.7), protonMass(938.3);
    CartesianVector muonMomentum(GetApproximatePfoMomentum(pLongestPfo, muonMass)); 

    pandora::CartesianVector neutrinoMomentum = muonMomentum;

    for (const auto pPfo : *pPfoList)
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

pandora::CartesianVector EventSelectionAlgorithm::GetApproximatePfoMomentum(const pandora::ParticleFlowObject* pPfo, const float &particleMass) const
{
    LArTrackStateVector trackStateVector; 
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, LArPfoHelper::GetVertex(pPfo), 20, LArGeometryHelper::GetWireZPitch(this->GetPandora()), trackStateVector);

    const pandora::CartesianVector particleDirection(trackStateVector.front().GetDirection().GetUnitVector());
    float particleCharge(GetPfoCharge(pPfo));

    const float adcToMeV(0.0000236 * 197 * (1.0/0.62));
    const float momentumNorm(std::sqrt(adcToMeV * adcToMeV * particleCharge * particleCharge + 2 * adcToMeV * particleCharge * particleMass));

    pandora::CartesianVector particleMomentum(particleDirection.GetX() * momentumNorm, particleDirection.GetY() * momentumNorm, particleDirection.GetZ() * momentumNorm);
    return particleMomentum;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
    "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
    "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
    "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
    "WriteToTree", m_writeToTree));
    
    if (m_writeToTree)
    {    
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));
    }   

    AlgorithmTool *pAlgorithmTool(nullptr);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "TrackDirection", pAlgorithmTool));

    if (!(m_pTrackDirectionTool = dynamic_cast<TrackDirectionTool *>(pAlgorithmTool)))
        throw STATUS_CODE_FAILURE;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} //namespace lar_content
