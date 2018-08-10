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
    m_writeToTree(false),
    m_viewEvent(false),
    m_fileIdentifier(0),
    m_eventNumber(-1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventSelectionAlgorithm::~EventSelectionAlgorithm()
{
    if (m_writeToTree)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
        //PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "EventDescription", m_fileName.c_str(), "UPDATE"));
        //PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "TrueNeutrinos", m_fileName.c_str(), "UPDATE"));
        //PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ReconstructedNeutrinos", m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventSelectionAlgorithm::Run()
{
    ++m_eventNumber;

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    this->WriteEventDescription(pMCParticleList, pPfoList, pCaloHitList);
    this->WriteVariables(pPfoList);
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));

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

    pandora::PfoList trueNeutrinoPfos(this->GetTrueNeutrinoAssociatedPfos(pPfoList));  
    int nTrueNeutrinoAssociatedTracks(this->CountTrackPfos(trueNeutrinoPfos));

    int firstNeutrinoPfoMCPDG(-1), secondNeutrinoPfoMCPDG(-1);
    int firstNeutrinoPfoParentMCPDG(-1), secondNeutrinoPfoParentMCPDG(-1);
   
    try
    {
        if (trueNeutrinoPfos.size() > 0)
        {
            firstNeutrinoPfoMCPDG = LArMCParticleHelper::GetMainMCParticle(trueNeutrinoPfos.front())->GetParticleId();
            firstNeutrinoPfoParentMCPDG = (LArMCParticleHelper::GetMainMCParticle(trueNeutrinoPfos.front()))->GetParentList().front()->GetParticleId();
        }

        if (trueNeutrinoPfos.size() > 1)
        {
            secondNeutrinoPfoMCPDG = LArMCParticleHelper::GetMainMCParticle(*(std::next(trueNeutrinoPfos.begin())))->GetParticleId();
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
    int nRecoNeutrinoAssociatedTracks(this->CountTrackPfos(recoNeutrinoPrimaryDaughters));

    int interactionTypeResemblance(this->GetInteractionTypeResemblance(recoNeutrinoPrimaryDaughters, nuMCParticlesToGoodHitsMap));

    int chosenSliceContainsNeutrino(trueInteractionType == interactionTypeResemblance ? 1 : 0);
    int nSliceCosmicRays(this->GetNumberCosmicRaysChosenSlice(recoNeutrinoPrimaryDaughters));

    pandora::PfoList allConnectedPfos;
    LArPfoHelper::GetAllConnectedPfos(neutrinoPfos, allConnectedPfos);

    pandora::CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(allConnectedPfos, TPC_3D, caloHitList);
    int chosenSliceNumberHits(caloHitList.size());

    std::vector<int> trueNeutrinoHitsVector(this->GetNeutrinoInducedHits(pMCParticleList, pPfoList, recoNeutrinoPrimaryDaughters, pCaloHitList));
 
    int nTrueNeutrinoInducedHits(trueNeutrinoHitsVector.at(0)), nPfoClustersTrueNeutrinoHits(trueNeutrinoHitsVector.at(1)), nPfoIsolatedTrueNeutrinoHits(trueNeutrinoHitsVector.at(2)), nPfoAllTrueNeutrinoHits(trueNeutrinoHitsVector.at(3)), nRecoNuClustersTrueNeutrinoHits(trueNeutrinoHitsVector.at(4)),  nRecoNuIsolatedTrueNeutrinoHits(trueNeutrinoHitsVector.at(5)), nRecoNuAllTrueNeutrinoHits(trueNeutrinoHitsVector.at(6));

    bool taggingFailure(this->IsTaggingFailure(pMCParticleList, pPfoList, pCaloHitList, nRecoNuAllTrueNeutrinoHits));

    std::cout << "Reconstructed neutrino particle multiplicity: " << recoNeutrinoPrimaryDaughters.size() << std::endl;
    std::cout << "True interaction type: " << LArInteractionTypeHelper::ToString(static_cast<LArInteractionTypeHelper::InteractionType>(trueInteractionType)) << std::endl;
    std::cout << "Number of cosmic rays: " << nSliceCosmicRays << std::endl;
    std::cout << "Tagging failure: " << taggingFailure << std::endl;
    std::cout << "firstNeutrinoPfoMCPDG: " << firstNeutrinoPfoMCPDG << std::endl;
    std::cout << "firstNeutrinoPfoParentMCPDG: " << firstNeutrinoPfoParentMCPDG << std::endl;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TrueInteractionType", trueInteractionType));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NeutrinoNuanceCode", neutrinoNuanceCode));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "FirstNeutrinoPfoMCPDG", firstNeutrinoPfoMCPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecondNeutrinoPfoMCPDG", secondNeutrinoPfoMCPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "FirstNeutrinoPfoParentMCPDG", firstNeutrinoPfoParentMCPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecondNeutrinoPfoParentMCPDG", secondNeutrinoPfoParentMCPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TrueNeutrinoNumberAssociatedParticles", (int)trueNeutrinoPfos.size()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TrueNeutrinoNumberAssociatedTracks", nTrueNeutrinoAssociatedTracks));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TrueNeutrinoNumberInducedHits", nTrueNeutrinoInducedHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PfoClustersTrueNeutrinoHits", nPfoClustersTrueNeutrinoHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PfoIsolatedTrueNeutrinoHits", nPfoIsolatedTrueNeutrinoHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PfoAllTrueNeutrinoHits", nPfoAllTrueNeutrinoHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "RecoNuClustersTrueNeutrinoHits", nRecoNuClustersTrueNeutrinoHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "RecoNuIsolatedTrueNeutrinoHits", nRecoNuIsolatedTrueNeutrinoHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "RecoNuAllTrueNeutrinoHits", nRecoNuAllTrueNeutrinoHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "RecoNeutrinoNumberAssociatedParticles", (int)recoNeutrinoPrimaryDaughters.size()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "RecoNeutrinoNumberAssociatedTracks", nRecoNeutrinoAssociatedTracks));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TaggingFailure", taggingFailure ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChosenSliceInteractionType", interactionTypeResemblance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChosenSliceContainsTrueNeutrino", chosenSliceContainsNeutrino));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChosenSliceNumberCosmicRays", nSliceCosmicRays));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChosenSliceNumberHits", chosenSliceNumberHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "FileIdentifier", m_fileIdentifier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "EventNumber", m_eventNumber));
    //PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
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
    return static_cast<int>(static_cast<unsigned int>(interactionType));
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<int> EventSelectionAlgorithm::GetNeutrinoInducedHits(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList, pandora::PfoList &recoNeutrinoPrimaryDaughters, const pandora::CaloHitList* pCaloHitList) const
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

    std::cout << "------------------------------------" << std::endl;
    std::cout << "nTrueNeutrinoInducedHits: " << nTrueNeutrinoInducedHits << std::endl;
    //std::cout << "nPfoClustersTrueNeutrinoHits: " << nPfoClustersTrueNeutrinoHits << std::endl;
    //std::cout << "nPfoIsolatedTrueNeutrinoHits: " << nPfoIsolatedTrueNeutrinoHits << std::endl;
    std::cout << "nPfoAllTrueNeutrinoHits: " << nPfoAllTrueNeutrinoHits << std::endl;
    //std::cout << "nRecoNuClustersTrueNeutrinoHits: " << nRecoNuClustersTrueNeutrinoHits << std::endl;
    //std::cout << "nRecoNuIsolatedTrueNeutrinoHits: " << nRecoNuIsolatedTrueNeutrinoHits << std::endl;
    std::cout << "nRecoNuAllTrueNeutrinoHits: " << nRecoNuAllTrueNeutrinoHits << std::endl;
    std::cout << "------------------------------------" << std::endl;


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

void EventSelectionAlgorithm::GetClusterHits(pandora::PfoList &pfoList, pandora::CaloHitList &caloHitList) const
{
    LArPfoHelper::GetCaloHits(pfoList, TPC_VIEW_U, caloHitList);
    LArPfoHelper::GetCaloHits(pfoList, TPC_VIEW_V, caloHitList);
    LArPfoHelper::GetCaloHits(pfoList, TPC_VIEW_W, caloHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSelectionAlgorithm::GetIsolatedHits(pandora::PfoList &pfoList, pandora::CaloHitList &caloHitList) const
{
    LArPfoHelper::GetIsolatedCaloHits(pfoList, TPC_VIEW_U, caloHitList);
    LArPfoHelper::GetIsolatedCaloHits(pfoList, TPC_VIEW_V, caloHitList);
    LArPfoHelper::GetIsolatedCaloHits(pfoList, TPC_VIEW_W, caloHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSelectionAlgorithm::MergeHitLists(pandora::CaloHitList &clusterHitList, pandora::CaloHitList &isolatedHitList, pandora::CaloHitList &combinedHitList) const
{
    combinedHitList = clusterHitList;
    pandora::CaloHitList isolatedCaloHitsCopy = isolatedHitList;
    combinedHitList.splice(combinedHitList.begin(), isolatedCaloHitsCopy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CaloHitList EventSelectionAlgorithm::GetTrueNeutrinoHits(pandora::CaloHitList &hitList, LArMCParticleHelper::CaloHitToMCMap &hitToMCMap) const
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

bool EventSelectionAlgorithm::IsTaggingFailure(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList, const pandora::CaloHitList* pCaloHitList, int nRecoNuAllTrueNeutrinoHits) const 
{
    LArMCParticleHelper::MCRelationMap mcPrimaryMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);

    LArMCParticleHelper::CaloHitToMCMap hitToMCMap;
    LArMCParticleHelper::MCContributionMap mcToTrueHitListMap; 
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcPrimaryMap, hitToMCMap, mcToTrueHitListMap);

    bool taggingFailure(false);

    std::cout << "------------------------------------" << std::endl;

    for (const auto pPfo : *pPfoList)
    {
        pandora::PfoList thisPfo;
        thisPfo.insert(thisPfo.begin(), pPfo);

        if (LArPfoHelper::IsNeutrino(pPfo))
            std::cout << "PFO hypothesis: NEUTRINO" << std::endl;
        else
            std::cout << "PFO hypothesis: COSMIC RAY" << std::endl;

        pandora::CaloHitList thisPfoAllHits, thisPfoClusterHits, thisPfoIsolatedHits;

        this->GetClusterHits(thisPfo, thisPfoClusterHits);
        this->GetIsolatedHits(thisPfo, thisPfoIsolatedHits);
        this->MergeHitLists(thisPfoClusterHits, thisPfoIsolatedHits, thisPfoAllHits);
        
        int nThisPfoClustersTrueNeutrinoHits(this->GetTrueNeutrinoHits(thisPfoClusterHits, hitToMCMap).size()), nThisPfoIsolatedTrueNeutrinoHits(this->GetTrueNeutrinoHits(thisPfoIsolatedHits, hitToMCMap).size()), nThisPfoAllTrueNeutrinoHits(this->GetTrueNeutrinoHits(thisPfoAllHits, hitToMCMap).size());

        if (!LArPfoHelper::IsNeutrino(pPfo) && nThisPfoAllTrueNeutrinoHits > nRecoNuAllTrueNeutrinoHits)
            taggingFailure = true;

        std::cout << "This PFO total number true nu hits: " << nThisPfoClustersTrueNeutrinoHits << std::endl;
        std::cout << "This PFO number isolated true nu hits: " << nThisPfoIsolatedTrueNeutrinoHits << std::endl;
        std::cout << "This PFO number cluster true nu hits: " << nThisPfoAllTrueNeutrinoHits << std::endl;
        std::cout << "********************************************" << std::endl;
    } 

    return taggingFailure;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::PfoList EventSelectionAlgorithm::GetTrueNeutrinoAssociatedPfos(const pandora::PfoList* pPfoList) const
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

pandora::PfoList EventSelectionAlgorithm::GetPrimaryDaughters(pandora::PfoList &neutrinoPfos) const
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

int EventSelectionAlgorithm::GetInteractionTypeResemblance(pandora::PfoList &primaryDaughters, LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap) const
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

int EventSelectionAlgorithm::GetNumberCosmicRaysChosenSlice(pandora::PfoList &primaryDaughters) const
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

void EventSelectionAlgorithm::WriteVariables(const pandora::PfoList* pPfoList) const
{
    PfoList neutrinoPfos;
    LArPfoHelper::GetRecoNeutrinos(pPfoList, neutrinoPfos);

    if (neutrinoPfos.size() != 1) return;

    pandora::PfoList recoNeutrinoPrimaryDaughters(this->GetPrimaryDaughters(neutrinoPfos));

    //Get interaction type
    int interactionType(GetInteractionType());
    std::cout << "Interaction type: " << interactionType << std::endl;

    //Get number of tracks and showers
    int nTracks(0), nShowers(0);
    GetNumberTracksAndShowers(recoNeutrinoPrimaryDaughters, nTracks, nShowers);

    //Get total event charge
    float totalEventCharge(0.f);

    for (const auto pPfo : recoNeutrinoPrimaryDaughters)
        totalEventCharge += GetPfoCharge(pPfo);

    //Get shortest and longest PFO
    const pandora::ParticleFlowObject* pLongestPfo(GetLongestPfo(recoNeutrinoPrimaryDaughters));
    const pandora::ParticleFlowObject* pShortestPfo(GetShortestPfo(recoNeutrinoPrimaryDaughters));

    //Longest and shortest pfo total charges
    float shortestPfoCharge(GetPfoCharge(pShortestPfo));
    float longestPfoCharge(GetPfoCharge(pLongestPfo));

    //longest pfo opening angle with beamline
    float thetaBeamPfo(GetThetaBeamPfo(pLongestPfo));

    //longest and shortest pfo opening angle
    float phiPfoOpeningAngle(GetPfoOpeningAngle(pLongestPfo, pShortestPfo));

    //longest and shortest pfo lengths 
    float shortestPfoLength(std::sqrt(LArPfoHelper::GetThreeDLengthSquared(pShortestPfo))), longestPfoLength(std::sqrt(LArPfoHelper::GetThreeDLengthSquared(pLongestPfo)));

    //longitudinal and transverse momenta
    pandora::CartesianVector neutrinoMomentum(GetApproximateNeutrinoMomentum(recoNeutrinoPrimaryDaughters, pLongestPfo));

    //dE/dx fit variables
    TrackDirectionTool::DirectionFitObject fitResult = m_pTrackDirectionTool->GetPfoDirection(pLongestPfo);
    TrackDirectionTool::FitParameters fitParameters(fitResult.GetFitParameters());

    //Write all information to tree
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "InteractionType", interactionType));
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "IsSingleMuon", (interactionType == 0 ? 1 : 0)));
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "IsMuonProton", (interactionType == 1 ? 1 : 0)));
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NumberTracks", nTracks));
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NumberShowers", nShowers));
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
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "FileIdentifier", m_fileIdentifier));
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "EventNumber", m_eventNumber));
    //PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));

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

pandora::CartesianVector EventSelectionAlgorithm::GetApproximateNeutrinoMomentum(pandora::PfoList pfoList, const pandora::ParticleFlowObject* pLongestPfo) const
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

pandora::CartesianVector EventSelectionAlgorithm::GetApproximatePfoMomentum(const pandora::ParticleFlowObject* pPfo, const float &particleMass) const
{
    LArTrackStateVector trackStateVector; 
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, LArPfoHelper::GetVertex(pPfo), 20, LArGeometryHelper::GetWireZPitch(this->GetPandora()), trackStateVector);

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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
    "ViewEvent", m_viewEvent));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
    "FileIdentifier", m_fileIdentifier));

    AlgorithmTool *pAlgorithmTool(nullptr);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "TrackDirection", pAlgorithmTool));

    if (!(m_pTrackDirectionTool = dynamic_cast<TrackDirectionTool *>(pAlgorithmTool)))
        throw STATUS_CODE_FAILURE;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} //namespace lar_content
