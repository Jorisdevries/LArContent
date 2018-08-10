/**
 *  @file   ExampleContent/include/ExampleAlgorithms/DirectionAnalysisAlgorithm.h
 * 
 *  @brief  Header file for the access lists algorithm class.
 * 
 *  $Log: $
 */

#ifndef EVENT_SELECTION_ALGORITHM_H
#define EVENT_SELECTION_ALGORITHM_H

#include "Pandora/Algorithm.h"
#include "larpandoracontent/LArDirection/TrackDirectionTool.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief EventSelectionAlgorithm class
 */

class EventSelectionAlgorithm : public pandora::Algorithm
{
    public:
        /**
         *  @brief Factory class that creates algorithm instances
         */

        class Factory : public pandora::AlgorithmFactory
        {
            public:
            pandora::Algorithm* CreateAlgorithm() const;
        };

        EventSelectionAlgorithm();

        ~EventSelectionAlgorithm();
    
    protected:

    private:
        pandora::StatusCode Run();

        void WriteEventDescription(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList, const pandora::CaloHitList* pCaloHitList) const;

        int GetInteractionType(LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap) const;

        std::vector<int> GetNeutrinoInducedHits(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList, pandora::PfoList &recoNeutrinoPrimaryDaughters, const pandora::CaloHitList* pCaloHitList) const;

        void GetClusterHits(pandora::PfoList &pfoList, pandora::CaloHitList &caloHitList) const;

        void GetIsolatedHits(pandora::PfoList &pfoList, pandora::CaloHitList &caloHitList) const;

        void MergeHitLists(pandora::CaloHitList &clusterHitList, pandora::CaloHitList &isolatedHitList, pandora::CaloHitList &combinedHitList) const;

        pandora::CaloHitList GetTrueNeutrinoHits(pandora::CaloHitList &hitList, LArMCParticleHelper::CaloHitToMCMap &hitToMCMap) const;

        bool IsTaggingFailure(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList, const pandora::CaloHitList* pCaloHitList, int nRecoNuAllTrueNeutrinoHits) const;

        pandora::PfoList GetTrueNeutrinoAssociatedPfos(const pandora::PfoList* pPfoList) const;

        int CountTrackPfos(pandora::PfoList &pfoList) const;

        pandora::PfoList GetPrimaryDaughters(pandora::PfoList &neutrinoPfos) const;

        int GetInteractionTypeResemblance(pandora::PfoList &primaryDaughters, LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap) const;

        int GetNumberCosmicRaysChosenSlice(pandora::PfoList &primaryDaughters) const;

        //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        void WriteNeutrinoInformation(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList) const;

        void WriteTrueNeutrinoInformation(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList) const;

        void WriteReconstructedNeutrinoInformation(const pandora::MCParticleList *pMCParticleList, const pandora::PfoList* pPfoList) const;

        static const pandora::MCParticle *GetMainMCParticle(const pandora::ParticleFlowObject *const pPfo);

        //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
        void WriteVariables(const pandora::PfoList* pPfoList) const;

        int GetInteractionType() const;

        void GetNumberTracksAndShowers(pandora::PfoList pfoList, int &nTracks, int &nShowers) const;

        float GetTotalEventCharge(const pandora::CaloHitList *const pCaloHitList) const;

        const pandora::ParticleFlowObject* GetLongestPfo(pandora::PfoList pfoList) const;

        const pandora::ParticleFlowObject* GetShortestPfo(pandora::PfoList pfoList) const;

        float GetPfoCharge(const pandora::ParticleFlowObject* pPfo) const;

        float GetThetaBeamPfo(const pandora::ParticleFlowObject* pPfo) const;

        float GetPfoOpeningAngle(const pandora::ParticleFlowObject* pPfo1, const pandora::ParticleFlowObject* pPfo2) const;

        pandora::CartesianVector GetApproximateNeutrinoMomentum(pandora::PfoList pfoList, const pandora::ParticleFlowObject* pLongestPfo) const; 

        pandora::CartesianVector GetApproximatePfoMomentum(const pandora::ParticleFlowObject* pLongestPfo, const float &particleMass) const; 

        pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

        std::string         m_mcParticleListName;
        std::string         m_caloHitListName;
        std::string         m_pfoListName;
        
        bool                m_writeToTree;
        std::string         m_treeName;
        std::string         m_fileName;

        bool                m_viewEvent;

        int                     m_fileIdentifier;
        int                     m_eventNumber;

        TrackDirectionTool  *m_pTrackDirectionTool;
};

inline pandora::Algorithm* EventSelectionAlgorithm::Factory::CreateAlgorithm() const
{
    return new EventSelectionAlgorithm();
}

} //namespace lar_content

#endif //#ifndef EVENT_SELECTION_ALGORITHM_H 
