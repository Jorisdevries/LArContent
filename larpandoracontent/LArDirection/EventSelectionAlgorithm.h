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

        void GetNumberTracksAndShowers(const pandora::PfoList* const pPfoList, int &nTracks, int &nShowers) const;

        float GetTotalEventCharge(const pandora::CaloHitList *const pCaloHitList) const;

        const pandora::ParticleFlowObject* GetLongestPfo(const pandora::PfoList* pPfoList) const;

        const pandora::ParticleFlowObject* GetShortestPfo(const pandora::PfoList* pPfoList) const;

        float GetPfoCharge(const pandora::ParticleFlowObject* pPfo) const;

        float GetThetaBeamPfo(const pandora::ParticleFlowObject* pPfo) const;

        float GetPfoOpeningAngle(const pandora::ParticleFlowObject* pPfo1, const pandora::ParticleFlowObject* pPfo2) const;

        pandora::CartesianVector GetApproximateNeutrinoMomentum(const pandora::PfoList* pPfoList, const pandora::ParticleFlowObject* pLongestPfo, pandora::CartesianVector neutrinoMomentum) const; 

        pandora::CartesianVector GetApproximatePfoMomentum(const pandora::ParticleFlowObject* pLongestPfo, const float &particleMass) const; 

        pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

        std::string         m_mcParticleListName;
        std::string         m_inputHitListName;
        std::string         m_pfoListName;
        
        bool                m_writeToTree;
        std::string         m_treeName;
        std::string         m_fileName;

        TrackDirectionTool  *m_pTrackDirectionTool;
};

inline pandora::Algorithm* EventSelectionAlgorithm::Factory::CreateAlgorithm() const
{
    return new EventSelectionAlgorithm();
}

} //namespace lar_content

#endif //#ifndef EVENT_SELECTION_ALGORITHM_H 
