/**
 *  @file   larpandoracontent/LArHelpers/LArInteractionTypeHelper.cc
 *
 *  @brief  Implementation of the interaction type helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include <regex>

using namespace pandora;

namespace lar_content
{

LArInteractionTypeHelper::InteractionType LArInteractionTypeHelper::GetInteractionType(const MCParticleList &mcPrimaryList)
{
    if (mcPrimaryList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    unsigned int nNonNeutrons(0), nMuons(0), nElectrons(0), nProtons(0), nPiPlus(0), nPiMinus(0), nPhotons(0), nKaonPlus(0), nKaonMinus(0);

    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        if (2112 != pMCPrimary->GetParticleId()) ++nNonNeutrons;
        if (13 == std::fabs(pMCPrimary->GetParticleId())) ++nMuons;
        if (11 == std::fabs(pMCPrimary->GetParticleId())) ++nElectrons;
        else if (2212 == std::fabs(pMCPrimary->GetParticleId())) ++nProtons;
        else if (22 == pMCPrimary->GetParticleId()) ++nPhotons;
        else if (211 == pMCPrimary->GetParticleId()) ++nPiPlus;
        else if (-211 == pMCPrimary->GetParticleId()) ++nPiMinus;
        else if (321 == pMCPrimary->GetParticleId()) ++nKaonPlus;
        else if (-321 == pMCPrimary->GetParticleId()) ++nKaonMinus;
    }

    if ((1 == mcPrimaryList.size()) && LArMCParticleHelper::IsCosmicRay(mcPrimaryList.front()))
    {
        if (1 == nMuons) return COSMIC_RAY_MU;
        if (1 == nProtons) return COSMIC_RAY_P;
        if (1 == nElectrons) return COSMIC_RAY_E;
        if (1 == nPhotons) return COSMIC_RAY_PHOTON;
        else return COSMIC_RAY_OTHER;
    }

    if ((1 == mcPrimaryList.size()) && LArMCParticleHelper::IsBeamParticle(mcPrimaryList.front()))
    {
        if (1 == nMuons) return BEAM_PARTICLE_MU;
        if (1 == nProtons) return BEAM_PARTICLE_P;
        if (1 == nElectrons) return BEAM_PARTICLE_E;
        if (1 == nPhotons) return BEAM_PARTICLE_PHOTON;
        if (1 == nPiPlus) return BEAM_PARTICLE_PI_PLUS;
        if (1 == nPiMinus) return BEAM_PARTICLE_PI_MINUS;
        if (1 == nKaonPlus) return BEAM_PARTICLE_KAON_PLUS;
        if (1 == nKaonMinus) return BEAM_PARTICLE_KAON_MINUS;
        else return BEAM_PARTICLE_OTHER;
    }

    const MCParticle *pMCNeutrino(nullptr);

    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        if (!LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) || (pMCNeutrino && (pMCNeutrino != LArMCParticleHelper::GetParentMCParticle(pMCPrimary))))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pMCNeutrino = LArMCParticleHelper::GetParentMCParticle(pMCPrimary);
    }

    if (!pMCNeutrino)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const int nuNuanceCode(LArMCParticleHelper::GetNuanceCode(pMCNeutrino));

    if (1001 == nuNuanceCode)
    {
        if ((1 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons)) return CCQEL_MU;
        if ((2 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons)) return CCQEL_MU_P;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons)) return CCQEL_MU_P_P;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons)) return CCQEL_MU_P_P_P;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons)) return CCQEL_MU_P_P_P_P;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons)) return CCQEL_MU_P_P_P_P_P;

        if ((1 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons)) return CCQEL_E;
        if ((2 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons)) return CCQEL_E_P;
        if ((3 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons)) return CCQEL_E_P_P;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons)) return CCQEL_E_P_P_P;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons)) return CCQEL_E_P_P_P_P;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons)) return CCQEL_E_P_P_P_P_P;
    }

    if (1002 == nuNuanceCode)
    {
        if ((1 == nNonNeutrons) && (1 == nProtons)) return NCQEL_P;
        if ((2 == nNonNeutrons) && (2 == nProtons)) return NCQEL_P_P;
        if ((3 == nNonNeutrons) && (3 == nProtons)) return NCQEL_P_P_P;
        if ((4 == nNonNeutrons) && (4 == nProtons)) return NCQEL_P_P_P_P;
        if ((5 == nNonNeutrons) && (5 == nProtons)) return NCQEL_P_P_P_P_P;
    }

    if ((nuNuanceCode >= 1003) && (nuNuanceCode <= 1005))
    {
        if ((1 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons)) return CCRES_MU;
        if ((2 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons)) return CCRES_MU_P;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons)) return CCRES_MU_P_P;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons)) return CCRES_MU_P_P_P;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons)) return CCRES_MU_P_P_P_P;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons)) return CCRES_MU_P_P_P_P_P;

        if ((2 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (1 == nPiPlus)) return CCRES_MU_PIPLUS;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_PIPLUS;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_P_PIPLUS;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_P_P_PIPLUS;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_P_P_P_PIPLUS;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_P_P_P_P_PIPLUS;

        if ((2 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (1 == nPhotons)) return CCRES_MU_PHOTON;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_PHOTON;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_P_PHOTON;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_P_P_PHOTON;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_P_P_P_PHOTON;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_P_P_P_P_PHOTON;

        if ((3 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (2 == nPhotons)) return CCRES_MU_PIZERO;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_PIZERO;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_P_PIZERO;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_P_P_PIZERO;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_P_P_P_PIZERO;
        if ((8 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_P_P_P_P_PIZERO;

        if ((1 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons)) return CCRES_E;
        if ((2 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons)) return CCRES_E_P;
        if ((3 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons)) return CCRES_E_P_P;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons)) return CCRES_E_P_P_P;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons)) return CCRES_E_P_P_P_P;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons)) return CCRES_E_P_P_P_P_P;

        if ((2 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons) && (1 == nPiPlus)) return CCRES_E_PIPLUS;
        if ((3 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_PIPLUS;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_P_PIPLUS;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_P_P_PIPLUS;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_P_P_P_PIPLUS;
        if ((7 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_P_P_P_P_PIPLUS;

        if ((2 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons) && (1 == nPhotons)) return CCRES_E_PHOTON;
        if ((3 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons) && (1 == nPhotons)) return CCRES_E_P_PHOTON;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons) && (1 == nPhotons)) return CCRES_E_P_P_PHOTON;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons) && (1 == nPhotons)) return CCRES_E_P_P_P_PHOTON;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons) && (1 == nPhotons)) return CCRES_E_P_P_P_P_PHOTON;
        if ((7 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons) && (1 == nPhotons)) return CCRES_E_P_P_P_P_P_PHOTON;

        if ((3 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons) && (2 == nPhotons)) return CCRES_E_PIZERO;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons) && (2 == nPhotons)) return CCRES_E_P_PIZERO;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons) && (2 == nPhotons)) return CCRES_E_P_P_PIZERO;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons) && (2 == nPhotons)) return CCRES_E_P_P_P_PIZERO;
        if ((7 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons) && (2 == nPhotons)) return CCRES_E_P_P_P_P_PIZERO;
        if ((8 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons) && (2 == nPhotons)) return CCRES_E_P_P_P_P_P_PIZERO;
    }

    if ((nuNuanceCode >= 1006) && (nuNuanceCode <= 1009))
    {
        if ((1 == nNonNeutrons) && (1 == nProtons)) return NCRES_P;
        if ((2 == nNonNeutrons) && (2 == nProtons)) return NCRES_P_P;
        if ((3 == nNonNeutrons) && (3 == nProtons)) return NCRES_P_P_P;
        if ((4 == nNonNeutrons) && (4 == nProtons)) return NCRES_P_P_P_P;
        if ((5 == nNonNeutrons) && (5 == nProtons)) return NCRES_P_P_P_P_P;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPiPlus)) return NCRES_PIPLUS;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPiPlus)) return NCRES_P_PIPLUS;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPiPlus)) return NCRES_P_P_PIPLUS;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPiPlus)) return NCRES_P_P_P_PIPLUS;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPiPlus)) return NCRES_P_P_P_P_PIPLUS;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPiPlus)) return NCRES_P_P_P_P_P_PIPLUS;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPiMinus)) return NCRES_PIMINUS;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPiMinus)) return NCRES_P_PIMINUS;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPiMinus)) return NCRES_P_P_PIMINUS;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPiMinus)) return NCRES_P_P_P_PIMINUS;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPiMinus)) return NCRES_P_P_P_P_PIMINUS;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPiMinus)) return NCRES_P_P_P_P_P_PIMINUS;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPhotons)) return NCRES_PHOTON;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPhotons)) return NCRES_P_PHOTON;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPhotons)) return NCRES_P_P_PHOTON;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPhotons)) return NCRES_P_P_P_PHOTON;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPhotons)) return NCRES_P_P_P_P_PHOTON;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPhotons)) return NCRES_P_P_P_P_P_PHOTON;

        if ((2 == nNonNeutrons) && (0 == nProtons) && (2 == nPhotons)) return NCRES_PIZERO;
        if ((3 == nNonNeutrons) && (1 == nProtons) && (2 == nPhotons)) return NCRES_P_PIZERO;
        if ((4 == nNonNeutrons) && (2 == nProtons) && (2 == nPhotons)) return NCRES_P_P_PIZERO;
        if ((5 == nNonNeutrons) && (3 == nProtons) && (2 == nPhotons)) return NCRES_P_P_P_PIZERO;
        if ((6 == nNonNeutrons) && (4 == nProtons) && (2 == nPhotons)) return NCRES_P_P_P_P_PIZERO;
        if ((7 == nNonNeutrons) && (5 == nProtons) && (2 == nPhotons)) return NCRES_P_P_P_P_P_PIZERO;
    }

    if (1091 == nuNuanceCode)
    {
        if ((1 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons)) return CCDIS_MU;
        if ((2 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons)) return CCDIS_MU_P;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons)) return CCDIS_MU_P_P;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons)) return CCDIS_MU_P_P_P;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons)) return CCDIS_MU_P_P_P_P;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons)) return CCDIS_MU_P_P_P_P_P;

        if ((2 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (1 == nPiPlus)) return CCDIS_MU_PIPLUS;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (1 == nPiPlus)) return CCDIS_MU_P_PIPLUS;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (1 == nPiPlus)) return CCDIS_MU_P_P_PIPLUS;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (1 == nPiPlus)) return CCDIS_MU_P_P_P_PIPLUS;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (1 == nPiPlus)) return CCDIS_MU_P_P_P_P_PIPLUS;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (1 == nPiPlus)) return CCDIS_MU_P_P_P_P_P_PIPLUS;

        if ((2 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (1 == nPhotons)) return CCDIS_MU_PHOTON;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (1 == nPhotons)) return CCDIS_MU_P_PHOTON;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (1 == nPhotons)) return CCDIS_MU_P_P_PHOTON;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (1 == nPhotons)) return CCDIS_MU_P_P_P_PHOTON;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (1 == nPhotons)) return CCDIS_MU_P_P_P_P_PHOTON;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (1 == nPhotons)) return CCDIS_MU_P_P_P_P_P_PHOTON;

        if ((3 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (2 == nPhotons)) return CCDIS_MU_PIZERO;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (2 == nPhotons)) return CCDIS_MU_P_PIZERO;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (2 == nPhotons)) return CCDIS_MU_P_P_PIZERO;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (2 == nPhotons)) return CCDIS_MU_P_P_P_PIZERO;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (2 == nPhotons)) return CCDIS_MU_P_P_P_P_PIZERO;
        if ((8 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (2 == nPhotons)) return CCDIS_MU_P_P_P_P_P_PIZERO;
    }

    if (1092 == nuNuanceCode)
    {
        if ((1 == nNonNeutrons) && (1 == nProtons)) return NCDIS_P;
        if ((2 == nNonNeutrons) && (2 == nProtons)) return NCDIS_P_P;
        if ((3 == nNonNeutrons) && (3 == nProtons)) return NCDIS_P_P_P;
        if ((4 == nNonNeutrons) && (4 == nProtons)) return NCDIS_P_P_P_P;
        if ((5 == nNonNeutrons) && (5 == nProtons)) return NCDIS_P_P_P_P_P;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPiPlus)) return NCDIS_PIPLUS;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPiPlus)) return NCDIS_P_PIPLUS;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPiPlus)) return NCDIS_P_P_PIPLUS;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPiPlus)) return NCDIS_P_P_P_PIPLUS;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPiPlus)) return NCDIS_P_P_P_P_PIPLUS;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPiPlus)) return NCDIS_P_P_P_P_P_PIPLUS;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPiMinus)) return NCDIS_PIMINUS;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPiMinus)) return NCDIS_P_PIMINUS;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPiMinus)) return NCDIS_P_P_PIMINUS;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPiMinus)) return NCDIS_P_P_P_PIMINUS;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPiMinus)) return NCDIS_P_P_P_P_PIMINUS;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPiMinus)) return NCDIS_P_P_P_P_P_PIMINUS;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPhotons)) return NCDIS_PHOTON;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPhotons)) return NCDIS_P_PHOTON;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPhotons)) return NCDIS_P_P_PHOTON;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPhotons)) return NCDIS_P_P_P_PHOTON;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPhotons)) return NCDIS_P_P_P_P_PHOTON;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPhotons)) return NCDIS_P_P_P_P_P_PHOTON;

        if ((2 == nNonNeutrons) && (0 == nProtons) && (2 == nPhotons)) return NCDIS_PIZERO;
        if ((3 == nNonNeutrons) && (1 == nProtons) && (2 == nPhotons)) return NCDIS_P_PIZERO;
        if ((4 == nNonNeutrons) && (2 == nProtons) && (2 == nPhotons)) return NCDIS_P_P_PIZERO;
        if ((5 == nNonNeutrons) && (3 == nProtons) && (2 == nPhotons)) return NCDIS_P_P_P_PIZERO;
        if ((6 == nNonNeutrons) && (4 == nProtons) && (2 == nPhotons)) return NCDIS_P_P_P_P_PIZERO;
        if ((7 == nNonNeutrons) && (5 == nProtons) && (2 == nPhotons)) return NCDIS_P_P_P_P_P_PIZERO;
    }

    if (1096 == nuNuanceCode) return NCCOH;
    if (1097 == nuNuanceCode) return CCCOH;

    //logic to deconstruct new nuance codes >= 5000 (if there are kaons we ignore this bit: kaons are not included in the interactionType enum anyway)
    if (nuNuanceCode >= 5000 && nKaonPlus == 0 && nKaonMinus == 0)
    {
        //deconstruct custom nuance code into CC/NC and mode
        int mutableNuNuanceCode(nuNuanceCode);
        std::vector<int> baseTenDeconstruction;

        while (mutableNuNuanceCode > 0)
        {
            baseTenDeconstruction.push_back(mutableNuNuanceCode % 10);
            mutableNuNuanceCode /= 10;
        }

        int CCNC(baseTenDeconstruction.at(2));
        int mode(10 * baseTenDeconstruction.at(1) + baseTenDeconstruction.at(0));
        std::string interactionTypeString("");

        //define CC/NC
        if (CCNC == 1)
            interactionTypeString.append("CC");
        else
            interactionTypeString.append("NC");
        
        //add interaction type, if not defined in modeToInteractionTypeMap bail out and return OTHER_INTERACTION
        std::map<int, std::string> modeToInteractionTypeMap = {{0, "QEL_"}, {1, "RES_"}, {2, "DIS_"}, {3, "COH_"}, {10, "MEC_"}};

        if (modeToInteractionTypeMap.find(mode) != modeToInteractionTypeMap.end())
            interactionTypeString.append(modeToInteractionTypeMap.at(mode));
        else
            return OTHER_INTERACTION;

        //add particles
        std::vector<std::pair<int, std::string>> particleCountsWithNames = {{nMuons, "MU_"}, {nElectrons, "E_"}, {nProtons, "P_"}, 
                                                                            {nPiPlus, "PIPLUS_"}, {nPiMinus, "PIMINUS_"}, {nPhotons, "PHOTON_"}};

        for (const auto &pair : particleCountsWithNames)
        {
            for (int i = 0; i < pair.first; ++i)
                interactionTypeString.append(pair.second);
        } 

        //remove trailing underscore
        interactionTypeString.erase(interactionTypeString.size() - 1);

        //map two photons to PIZERO, as before
        interactionTypeString = std::regex_replace(interactionTypeString, std::regex("PHOTON_PHOTON"), "PIZERO");
        return FromString(interactionTypeString); 
    }
    
    return OTHER_INTERACTION;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string LArInteractionTypeHelper::ToString(const InteractionType interactionType)
{
    switch (interactionType)
    {
    case CCQEL_MU: return "CCQEL_MU";
    case CCQEL_MU_P: return "CCQEL_MU_P";
    case CCQEL_MU_P_P: return "CCQEL_MU_P_P";
    case CCQEL_MU_P_P_P: return "CCQEL_MU_P_P_P";
    case CCQEL_MU_P_P_P_P: return "CCQEL_MU_P_P_P_P";
    case CCQEL_MU_P_P_P_P_P: return "CCQEL_MU_P_P_P_P_P";
    case CCQEL_E: return "CCQEL_E";
    case CCQEL_E_P: return "CCQEL_E_P";
    case CCQEL_E_P_P: return "CCQEL_E_P_P";
    case CCQEL_E_P_P_P: return "CCQEL_E_P_P_P";
    case CCQEL_E_P_P_P_P: return "CCQEL_E_P_P_P_P";
    case CCQEL_E_P_P_P_P_P: return "CCQEL_E_P_P_P_P_P";
    case NCQEL_P: return "NCQEL_P";
    case NCQEL_P_P: return "NCQEL_P_P";
    case NCQEL_P_P_P: return "NCQEL_P_P_P";
    case NCQEL_P_P_P_P: return "NCQEL_P_P_P_P";
    case NCQEL_P_P_P_P_P: return "NCQEL_P_P_P_P_P";
    case CCRES_MU: return "CCRES_MU";
    case CCRES_MU_P: return "CCRES_MU_P";
    case CCRES_MU_P_P: return "CCRES_MU_P_P";
    case CCRES_MU_P_P_P: return "CCRES_MU_P_P_P";
    case CCRES_MU_P_P_P_P: return "CCRES_MU_P_P_P_P";
    case CCRES_MU_P_P_P_P_P: return "CCRES_MU_P_P_P_P_P";
    case CCRES_MU_PIPLUS: return "CCRES_MU_PIPLUS";
    case CCRES_MU_P_PIPLUS: return "CCRES_MU_P_PIPLUS";
    case CCRES_MU_P_P_PIPLUS: return "CCRES_MU_P_P_PIPLUS";
    case CCRES_MU_P_P_P_PIPLUS: return "CCRES_MU_P_P_P_PIPLUS";
    case CCRES_MU_P_P_P_P_PIPLUS: return "CCRES_MU_P_P_P_P_PIPLUS";
    case CCRES_MU_P_P_P_P_P_PIPLUS: return "CCRES_MU_P_P_P_P_P_PIPLUS";
    case CCRES_MU_PHOTON: return "CCRES_MU_PHOTON";
    case CCRES_MU_P_PHOTON: return "CCRES_MU_P_PHOTON";
    case CCRES_MU_P_P_PHOTON: return "CCRES_MU_P_P_PHOTON";
    case CCRES_MU_P_P_P_PHOTON: return "CCRES_MU_P_P_P_PHOTON";
    case CCRES_MU_P_P_P_P_PHOTON: return "CCRES_MU_P_P_P_P_PHOTON";
    case CCRES_MU_P_P_P_P_P_PHOTON: return "CCRES_MU_P_P_P_P_P_PHOTON";
    case CCRES_MU_PIZERO: return "CCRES_MU_PIZERO";
    case CCRES_MU_P_PIZERO: return "CCRES_MU_P_PIZERO";
    case CCRES_MU_P_P_PIZERO: return "CCRES_MU_P_P_PIZERO";
    case CCRES_MU_P_P_P_PIZERO: return "CCRES_MU_P_P_P_PIZERO";
    case CCRES_MU_P_P_P_P_PIZERO: return "CCRES_MU_P_P_P_P_PIZERO";
    case CCRES_MU_P_P_P_P_P_PIZERO: return "CCRES_MU_P_P_P_P_P_PIZERO";
    case CCRES_E: return "CCRES_E";
    case CCRES_E_P: return "CCRES_E_P";
    case CCRES_E_P_P: return "CCRES_E_P_P";
    case CCRES_E_P_P_P: return "CCRES_E_P_P_P";
    case CCRES_E_P_P_P_P: return "CCRES_E_P_P_P_P";
    case CCRES_E_P_P_P_P_P: return "CCRES_E_P_P_P_P_P";
    case CCRES_E_PIPLUS: return "CCRES_E_PIPLUS";
    case CCRES_E_P_PIPLUS: return "CCRES_E_P_PIPLUS";
    case CCRES_E_P_P_PIPLUS: return "CCRES_E_P_P_PIPLUS";
    case CCRES_E_P_P_P_PIPLUS: return "CCRES_E_P_P_P_PIPLUS";
    case CCRES_E_P_P_P_P_PIPLUS: return "CCRES_E_P_P_P_P_PIPLUS";
    case CCRES_E_P_P_P_P_P_PIPLUS: return "CCRES_E_P_P_P_P_P_PIPLUS";
    case CCRES_E_PHOTON: return "CCRES_E_PHOTON";
    case CCRES_E_P_PHOTON: return "CCRES_E_P_PHOTON";
    case CCRES_E_P_P_PHOTON: return "CCRES_E_P_P_PHOTON";
    case CCRES_E_P_P_P_PHOTON: return "CCRES_E_P_P_P_PHOTON";
    case CCRES_E_P_P_P_P_PHOTON: return "CCRES_E_P_P_P_P_PHOTON";
    case CCRES_E_P_P_P_P_P_PHOTON: return "CCRES_E_P_P_P_P_P_PHOTON";
    case CCRES_E_PIZERO: return "CCRES_E_PIZERO";
    case CCRES_E_P_PIZERO: return "CCRES_E_P_PIZERO";
    case CCRES_E_P_P_PIZERO: return "CCRES_E_P_P_PIZERO";
    case CCRES_E_P_P_P_PIZERO: return "CCRES_E_P_P_P_PIZERO";
    case CCRES_E_P_P_P_P_PIZERO: return "CCRES_E_P_P_P_P_PIZERO";
    case CCRES_E_P_P_P_P_P_PIZERO: return "CCRES_E_P_P_P_P_P_PIZERO";
    case NCRES_P: return "NCRES_P";
    case NCRES_P_P: return "NCRES_P_P";
    case NCRES_P_P_P: return "NCRES_P_P_P";
    case NCRES_P_P_P_P: return "NCRES_P_P_P_P";
    case NCRES_P_P_P_P_P: return "NCRES_P_P_P_P_P";
    case NCRES_PIPLUS: return "NCRES_PIPLUS";
    case NCRES_P_PIPLUS: return "NCRES_P_PIPLUS";
    case NCRES_P_P_PIPLUS: return "NCRES_P_P_PIPLUS";
    case NCRES_P_P_P_PIPLUS: return "NCRES_P_P_P_PIPLUS";
    case NCRES_P_P_P_P_PIPLUS: return "NCRES_P_P_P_P_PIPLUS";
    case NCRES_P_P_P_P_P_PIPLUS: return "NCRES_P_P_P_P_P_PIPLUS";
    case NCRES_PIMINUS: return "NCRES_PIMINUS";
    case NCRES_P_PIMINUS: return "NCRES_P_PIMINUS";
    case NCRES_P_P_PIMINUS: return "NCRES_P_P_PIMINUS";
    case NCRES_P_P_P_PIMINUS: return "NCRES_P_P_P_PIMINUS";
    case NCRES_P_P_P_P_PIMINUS: return "NCRES_P_P_P_P_PIMINUS";
    case NCRES_P_P_P_P_P_PIMINUS: return "NCRES_P_P_P_P_P_PIMINUS";
    case NCRES_PHOTON: return "NCRES_PHOTON";
    case NCRES_P_PHOTON: return "NCRES_P_PHOTON";
    case NCRES_P_P_PHOTON: return "NCRES_P_P_PHOTON";
    case NCRES_P_P_P_PHOTON: return "NCRES_P_P_P_PHOTON";
    case NCRES_P_P_P_P_PHOTON: return "NCRES_P_P_P_P_PHOTON";
    case NCRES_P_P_P_P_P_PHOTON: return "NCRES_P_P_P_P_P_PHOTON";
    case NCRES_PIZERO: return "NCRES_PIZERO";
    case NCRES_P_PIZERO: return "NCRES_P_PIZERO";
    case NCRES_P_P_PIZERO: return "NCRES_P_P_PIZERO";
    case NCRES_P_P_P_PIZERO: return "NCRES_P_P_P_PIZERO";
    case NCRES_P_P_P_P_PIZERO: return "NCRES_P_P_P_P_PIZERO";
    case NCRES_P_P_P_P_P_PIZERO: return "NCRES_P_P_P_P_P_PIZERO";
    case CCDIS_MU: return "CCDIS_MU";
    case CCDIS_MU_P: return "CCDIS_MU_P";
    case CCDIS_MU_P_P: return "CCDIS_MU_P_P";
    case CCDIS_MU_P_P_P: return "CCDIS_MU_P_P_P";
    case CCDIS_MU_P_P_P_P: return "CCDIS_MU_P_P_P_P";
    case CCDIS_MU_P_P_P_P_P: return "CCDIS_MU_P_P_P_P_P";
    case CCDIS_MU_PIPLUS: return "CCDIS_MU_PIPLUS";
    case CCDIS_MU_P_PIPLUS: return "CCDIS_MU_P_PIPLUS";
    case CCDIS_MU_P_P_PIPLUS: return "CCDIS_MU_P_P_PIPLUS";
    case CCDIS_MU_P_P_P_PIPLUS: return "CCDIS_MU_P_P_P_PIPLUS";
    case CCDIS_MU_P_P_P_P_PIPLUS: return "CCDIS_MU_P_P_P_P_PIPLUS";
    case CCDIS_MU_P_P_P_P_P_PIPLUS: return "CCDIS_MU_P_P_P_P_P_PIPLUS";
    case CCDIS_MU_PHOTON: return "CCDIS_MU_PHOTON";
    case CCDIS_MU_P_PHOTON: return "CCDIS_MU_P_PHOTON";
    case CCDIS_MU_P_P_PHOTON: return "CCDIS_MU_P_P_PHOTON";
    case CCDIS_MU_P_P_P_PHOTON: return "CCDIS_MU_P_P_P_PHOTON";
    case CCDIS_MU_P_P_P_P_PHOTON: return "CCDIS_MU_P_P_P_P_PHOTON";
    case CCDIS_MU_P_P_P_P_P_PHOTON: return "CCDIS_MU_P_P_P_P_P_PHOTON";
    case CCDIS_MU_PIZERO: return "CCDIS_MU_PIZERO";
    case CCDIS_MU_P_PIZERO: return "CCDIS_MU_P_PIZERO";
    case CCDIS_MU_P_P_PIZERO: return "CCDIS_MU_P_P_PIZERO";
    case CCDIS_MU_P_P_P_PIZERO: return "CCDIS_MU_P_P_P_PIZERO";
    case CCDIS_MU_P_P_P_P_PIZERO: return "CCDIS_MU_P_P_P_P_PIZERO";
    case CCDIS_MU_P_P_P_P_P_PIZERO: return "CCDIS_MU_P_P_P_P_P_PIZERO";
    case NCDIS_P: return "NCDIS_P";
    case NCDIS_P_P: return "NCDIS_P_P";
    case NCDIS_P_P_P: return "NCDIS_P_P_P";
    case NCDIS_P_P_P_P: return "NCDIS_P_P_P_P";
    case NCDIS_P_P_P_P_P: return "NCDIS_P_P_P_P_P";
    case NCDIS_PIPLUS: return "NCDIS_PIPLUS";
    case NCDIS_P_PIPLUS: return "NCDIS_P_PIPLUS";
    case NCDIS_P_P_PIPLUS: return "NCDIS_P_P_PIPLUS";
    case NCDIS_P_P_P_PIPLUS: return "NCDIS_P_P_P_PIPLUS";
    case NCDIS_P_P_P_P_PIPLUS: return "NCDIS_P_P_P_P_PIPLUS";
    case NCDIS_P_P_P_P_P_PIPLUS: return "NCDIS_P_P_P_P_P_PIPLUS";
    case NCDIS_PIMINUS: return "NCDIS_PIMINUS";
    case NCDIS_P_PIMINUS: return "NCDIS_P_PIMINUS";
    case NCDIS_P_P_PIMINUS: return "NCDIS_P_P_PIMINUS";
    case NCDIS_P_P_P_PIMINUS: return "NCDIS_P_P_P_PIMINUS";
    case NCDIS_P_P_P_P_PIMINUS: return "NCDIS_P_P_P_P_PIMINUS";
    case NCDIS_P_P_P_P_P_PIMINUS: return "NCDIS_P_P_P_P_P_PIMINUS";
    case NCDIS_PHOTON: return "NCDIS_PHOTON";
    case NCDIS_P_PHOTON: return "NCDIS_P_PHOTON";
    case NCDIS_P_P_PHOTON: return "NCDIS_P_P_PHOTON";
    case NCDIS_P_P_P_PHOTON: return "NCDIS_P_P_P_PHOTON";
    case NCDIS_P_P_P_P_PHOTON: return "NCDIS_P_P_P_P_PHOTON";
    case NCDIS_P_P_P_P_P_PHOTON: return "NCDIS_P_P_P_P_P_PHOTON";
    case NCDIS_PIZERO: return "NCDIS_PIZERO";
    case NCDIS_P_PIZERO: return "NCDIS_P_PIZERO";
    case NCDIS_P_P_PIZERO: return "NCDIS_P_P_PIZERO";
    case NCDIS_P_P_P_PIZERO: return "NCDIS_P_P_P_PIZERO";
    case NCDIS_P_P_P_P_PIZERO: return "NCDIS_P_P_P_P_PIZERO";
    case NCDIS_P_P_P_P_P_PIZERO: return "NCDIS_P_P_P_P_P_PIZERO";
    case CCCOH: return "CCCOH";
    case NCCOH: return "NCCOH";
    case COSMIC_RAY_MU: return "COSMIC_RAY_MU";
    case COSMIC_RAY_P: return "COSMIC_RAY_P";
    case COSMIC_RAY_E: return "COSMIC_RAY_E";
    case COSMIC_RAY_PHOTON: return "COSMIC_RAY_PHOTON";
    case COSMIC_RAY_OTHER: return "COSMIC_RAY_OTHER";
    case BEAM_PARTICLE_MU: return "BEAM_PARTICLE_MU";
    case BEAM_PARTICLE_P: return "BEAM_PARTICLE_P";
    case BEAM_PARTICLE_E: return "BEAM_PARTICLE_E";
    case BEAM_PARTICLE_PHOTON: return "BEAM_PARTICLE_PHOTON";
    case BEAM_PARTICLE_PI_PLUS: return "BEAM_PARTICLE_PI_PLUS";
    case BEAM_PARTICLE_PI_MINUS: return "BEAM_PARTICLE_PI_MINUS";
    case BEAM_PARTICLE_KAON_PLUS: return "BEAM_PARTICLE_KAON_PLUS";
    case BEAM_PARTICLE_KAON_MINUS: return "BEAM_PARTICLE_KAON_MINUS";
    case BEAM_PARTICLE_OTHER: return "BEAM_PARTICLE_OTHER";
    case OTHER_INTERACTION: return "OTHER_INTERACTION";
    case ALL_INTERACTIONS: return "ALL_INTERACTIONS";
    case CCMEC_MU: return "CCMEC_MU";
    case CCMEC_MU_P: return "CCMEC_MU_P";
    case CCMEC_MU_P_P: return "CCMEC_MU_P_P";
    case CCMEC_MU_P_P_P: return "CCMEC_MU_P_P_P";
    case CCMEC_MU_P_P_P_P: return "CCMEC_MU_P_P_P_P";
    case CCMEC_MU_P_P_P_P_P: return "CCMEC_MU_P_P_P_P_P";
    case NCMEC_MU: return "NCMEC_MU";
    case NCMEC_MU_P: return "NCMEC_MU_P";
    case NCMEC_MU_P_P: return "NCMEC_MU_P_P";
    case NCMEC_MU_P_P_P: return "NCMEC_MU_P_P_P";
    case NCMEC_MU_P_P_P_P: return "NCMEC_MU_P_P_P_P";
    case NCMEC_MU_P_P_P_P_P: return "NCMEC_MU_P_P_P_P_P";
    case CCMEC_P: return "CCMEC_P";
    case NCMEC_P: return "NCMEC_P";
    case NCRES_MU: return "NCRES_MU";
    case NCRES_MU_P: return "NCRES_MU_P";
    default: return "UNKNOWN";
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArInteractionTypeHelper::InteractionType LArInteractionTypeHelper::FromString(std::string interactionTypeString)
{
    if (interactionTypeString == "CCQEL_MU") return CCQEL_MU;
    else if (interactionTypeString == "CCQEL_MU_P") return CCQEL_MU_P;
    else if (interactionTypeString == "CCQEL_MU_P_P") return CCQEL_MU_P_P;
    else if (interactionTypeString == "CCQEL_MU_P_P_P") return CCQEL_MU_P_P_P;
    else if (interactionTypeString == "CCQEL_MU_P_P_P_P") return CCQEL_MU_P_P_P_P;
    else if (interactionTypeString == "CCQEL_MU_P_P_P_P_P") return CCQEL_MU_P_P_P_P_P;
    else if (interactionTypeString == "CCQEL_E") return CCQEL_E;
    else if (interactionTypeString == "CCQEL_E_P") return CCQEL_E_P;
    else if (interactionTypeString == "CCQEL_E_P_P") return CCQEL_E_P_P;
    else if (interactionTypeString == "CCQEL_E_P_P_P") return CCQEL_E_P_P_P;
    else if (interactionTypeString == "CCQEL_E_P_P_P_P") return CCQEL_E_P_P_P_P;
    else if (interactionTypeString == "CCQEL_E_P_P_P_P_P") return CCQEL_E_P_P_P_P_P;
    else if (interactionTypeString == "NCQEL_P") return NCQEL_P;
    else if (interactionTypeString == "NCQEL_P_P") return NCQEL_P_P;
    else if (interactionTypeString == "NCQEL_P_P_P") return NCQEL_P_P_P;
    else if (interactionTypeString == "NCQEL_P_P_P_P") return NCQEL_P_P_P_P;
    else if (interactionTypeString == "NCQEL_P_P_P_P_P") return NCQEL_P_P_P_P_P;
    else if (interactionTypeString == "CCRES_MU") return CCRES_MU;
    else if (interactionTypeString == "CCRES_MU_P") return CCRES_MU_P;
    else if (interactionTypeString == "CCRES_MU_P_P") return CCRES_MU_P_P;
    else if (interactionTypeString == "CCRES_MU_P_P_P") return CCRES_MU_P_P_P;
    else if (interactionTypeString == "CCRES_MU_P_P_P_P") return CCRES_MU_P_P_P_P;
    else if (interactionTypeString == "CCRES_MU_P_P_P_P_P") return CCRES_MU_P_P_P_P_P;
    else if (interactionTypeString == "CCRES_MU_PIPLUS") return CCRES_MU_PIPLUS;
    else if (interactionTypeString == "CCRES_MU_P_PIPLUS") return CCRES_MU_P_PIPLUS;
    else if (interactionTypeString == "CCRES_MU_P_P_PIPLUS") return CCRES_MU_P_P_PIPLUS;
    else if (interactionTypeString == "CCRES_MU_P_P_P_PIPLUS") return CCRES_MU_P_P_P_PIPLUS;
    else if (interactionTypeString == "CCRES_MU_P_P_P_P_PIPLUS") return CCRES_MU_P_P_P_P_PIPLUS;
    else if (interactionTypeString == "CCRES_MU_P_P_P_P_P_PIPLUS") return CCRES_MU_P_P_P_P_P_PIPLUS;
    else if (interactionTypeString == "CCRES_MU_PHOTON") return CCRES_MU_PHOTON;
    else if (interactionTypeString == "CCRES_MU_P_PHOTON") return CCRES_MU_P_PHOTON;
    else if (interactionTypeString == "CCRES_MU_P_P_PHOTON") return CCRES_MU_P_P_PHOTON;
    else if (interactionTypeString == "CCRES_MU_P_P_P_PHOTON") return CCRES_MU_P_P_P_PHOTON;
    else if (interactionTypeString == "CCRES_MU_P_P_P_P_PHOTON") return CCRES_MU_P_P_P_P_PHOTON;
    else if (interactionTypeString == "CCRES_MU_P_P_P_P_P_PHOTON") return CCRES_MU_P_P_P_P_P_PHOTON;
    else if (interactionTypeString == "CCRES_MU_PIZERO") return CCRES_MU_PIZERO;
    else if (interactionTypeString == "CCRES_MU_P_PIZERO") return CCRES_MU_P_PIZERO;
    else if (interactionTypeString == "CCRES_MU_P_P_PIZERO") return CCRES_MU_P_P_PIZERO;
    else if (interactionTypeString == "CCRES_MU_P_P_P_PIZERO") return CCRES_MU_P_P_P_PIZERO;
    else if (interactionTypeString == "CCRES_MU_P_P_P_P_PIZERO") return CCRES_MU_P_P_P_P_PIZERO;
    else if (interactionTypeString == "CCRES_MU_P_P_P_P_P_PIZERO") return CCRES_MU_P_P_P_P_P_PIZERO;
    else if (interactionTypeString == "CCRES_E") return CCRES_E;
    else if (interactionTypeString == "CCRES_E_P") return CCRES_E_P;
    else if (interactionTypeString == "CCRES_E_P_P") return CCRES_E_P_P;
    else if (interactionTypeString == "CCRES_E_P_P_P") return CCRES_E_P_P_P;
    else if (interactionTypeString == "CCRES_E_P_P_P_P") return CCRES_E_P_P_P_P;
    else if (interactionTypeString == "CCRES_E_P_P_P_P_P") return CCRES_E_P_P_P_P_P;
    else if (interactionTypeString == "CCRES_E_PIPLUS") return CCRES_E_PIPLUS;
    else if (interactionTypeString == "CCRES_E_P_PIPLUS") return CCRES_E_P_PIPLUS;
    else if (interactionTypeString == "CCRES_E_P_P_PIPLUS") return CCRES_E_P_P_PIPLUS;
    else if (interactionTypeString == "CCRES_E_P_P_P_PIPLUS") return CCRES_E_P_P_P_PIPLUS;
    else if (interactionTypeString == "CCRES_E_P_P_P_P_PIPLUS") return CCRES_E_P_P_P_P_PIPLUS;
    else if (interactionTypeString == "CCRES_E_P_P_P_P_P_PIPLUS") return CCRES_E_P_P_P_P_P_PIPLUS;
    else if (interactionTypeString == "CCRES_E_PHOTON") return CCRES_E_PHOTON;
    else if (interactionTypeString == "CCRES_E_P_PHOTON") return CCRES_E_P_PHOTON;
    else if (interactionTypeString == "CCRES_E_P_P_PHOTON") return CCRES_E_P_P_PHOTON;
    else if (interactionTypeString == "CCRES_E_P_P_P_PHOTON") return CCRES_E_P_P_P_PHOTON;
    else if (interactionTypeString == "CCRES_E_P_P_P_P_PHOTON") return CCRES_E_P_P_P_P_PHOTON;
    else if (interactionTypeString == "CCRES_E_P_P_P_P_P_PHOTON") return CCRES_E_P_P_P_P_P_PHOTON;
    else if (interactionTypeString == "CCRES_E_PIZERO") return CCRES_E_PIZERO;
    else if (interactionTypeString == "CCRES_E_P_PIZERO") return CCRES_E_P_PIZERO;
    else if (interactionTypeString == "CCRES_E_P_P_PIZERO") return CCRES_E_P_P_PIZERO;
    else if (interactionTypeString == "CCRES_E_P_P_P_PIZERO") return CCRES_E_P_P_P_PIZERO;
    else if (interactionTypeString == "CCRES_E_P_P_P_P_PIZERO") return CCRES_E_P_P_P_P_PIZERO;
    else if (interactionTypeString == "CCRES_E_P_P_P_P_P_PIZERO") return CCRES_E_P_P_P_P_P_PIZERO;
    else if (interactionTypeString == "NCRES_P") return NCRES_P;
    else if (interactionTypeString == "NCRES_P_P") return NCRES_P_P;
    else if (interactionTypeString == "NCRES_P_P_P") return NCRES_P_P_P;
    else if (interactionTypeString == "NCRES_P_P_P_P") return NCRES_P_P_P_P;
    else if (interactionTypeString == "NCRES_P_P_P_P_P") return NCRES_P_P_P_P_P;
    else if (interactionTypeString == "NCRES_PIPLUS") return NCRES_PIPLUS;
    else if (interactionTypeString == "NCRES_P_PIPLUS") return NCRES_P_PIPLUS;
    else if (interactionTypeString == "NCRES_P_P_PIPLUS") return NCRES_P_P_PIPLUS;
    else if (interactionTypeString == "NCRES_P_P_P_PIPLUS") return NCRES_P_P_P_PIPLUS;
    else if (interactionTypeString == "NCRES_P_P_P_P_PIPLUS") return NCRES_P_P_P_P_PIPLUS;
    else if (interactionTypeString == "NCRES_P_P_P_P_P_PIPLUS") return NCRES_P_P_P_P_P_PIPLUS;
    else if (interactionTypeString == "NCRES_PIMINUS") return NCRES_PIMINUS;
    else if (interactionTypeString == "NCRES_P_PIMINUS") return NCRES_P_PIMINUS;
    else if (interactionTypeString == "NCRES_P_P_PIMINUS") return NCRES_P_P_PIMINUS;
    else if (interactionTypeString == "NCRES_P_P_P_PIMINUS") return NCRES_P_P_P_PIMINUS;
    else if (interactionTypeString == "NCRES_P_P_P_P_PIMINUS") return NCRES_P_P_P_P_PIMINUS;
    else if (interactionTypeString == "NCRES_P_P_P_P_P_PIMINUS") return NCRES_P_P_P_P_P_PIMINUS;
    else if (interactionTypeString == "NCRES_PHOTON") return NCRES_PHOTON;
    else if (interactionTypeString == "NCRES_P_PHOTON") return NCRES_P_PHOTON;
    else if (interactionTypeString == "NCRES_P_P_PHOTON") return NCRES_P_P_PHOTON;
    else if (interactionTypeString == "NCRES_P_P_P_PHOTON") return NCRES_P_P_P_PHOTON;
    else if (interactionTypeString == "NCRES_P_P_P_P_PHOTON") return NCRES_P_P_P_P_PHOTON;
    else if (interactionTypeString == "NCRES_P_P_P_P_P_PHOTON") return NCRES_P_P_P_P_P_PHOTON;
    else if (interactionTypeString == "NCRES_PIZERO") return NCRES_PIZERO;
    else if (interactionTypeString == "NCRES_P_PIZERO") return NCRES_P_PIZERO;
    else if (interactionTypeString == "NCRES_P_P_PIZERO") return NCRES_P_P_PIZERO;
    else if (interactionTypeString == "NCRES_P_P_P_PIZERO") return NCRES_P_P_P_PIZERO;
    else if (interactionTypeString == "NCRES_P_P_P_P_PIZERO") return NCRES_P_P_P_P_PIZERO;
    else if (interactionTypeString == "NCRES_P_P_P_P_P_PIZERO") return NCRES_P_P_P_P_P_PIZERO;
    else if (interactionTypeString == "CCDIS_MU") return CCDIS_MU;
    else if (interactionTypeString == "CCDIS_MU_P") return CCDIS_MU_P;
    else if (interactionTypeString == "CCDIS_MU_P_P") return CCDIS_MU_P_P;
    else if (interactionTypeString == "CCDIS_MU_P_P_P") return CCDIS_MU_P_P_P;
    else if (interactionTypeString == "CCDIS_MU_P_P_P_P") return CCDIS_MU_P_P_P_P;
    else if (interactionTypeString == "CCDIS_MU_P_P_P_P_P") return CCDIS_MU_P_P_P_P_P;
    else if (interactionTypeString == "CCDIS_MU_PIPLUS") return CCDIS_MU_PIPLUS;
    else if (interactionTypeString == "CCDIS_MU_P_PIPLUS") return CCDIS_MU_P_PIPLUS;
    else if (interactionTypeString == "CCDIS_MU_P_P_PIPLUS") return CCDIS_MU_P_P_PIPLUS;
    else if (interactionTypeString == "CCDIS_MU_P_P_P_PIPLUS") return CCDIS_MU_P_P_P_PIPLUS;
    else if (interactionTypeString == "CCDIS_MU_P_P_P_P_PIPLUS") return CCDIS_MU_P_P_P_P_PIPLUS;
    else if (interactionTypeString == "CCDIS_MU_P_P_P_P_P_PIPLUS") return CCDIS_MU_P_P_P_P_P_PIPLUS;
    else if (interactionTypeString == "CCDIS_MU_PHOTON") return CCDIS_MU_PHOTON;
    else if (interactionTypeString == "CCDIS_MU_P_PHOTON") return CCDIS_MU_P_PHOTON;
    else if (interactionTypeString == "CCDIS_MU_P_P_PHOTON") return CCDIS_MU_P_P_PHOTON;
    else if (interactionTypeString == "CCDIS_MU_P_P_P_PHOTON") return CCDIS_MU_P_P_P_PHOTON;
    else if (interactionTypeString == "CCDIS_MU_P_P_P_P_PHOTON") return CCDIS_MU_P_P_P_P_PHOTON;
    else if (interactionTypeString == "CCDIS_MU_P_P_P_P_P_PHOTON") return CCDIS_MU_P_P_P_P_P_PHOTON;
    else if (interactionTypeString == "CCDIS_MU_PIZERO") return CCDIS_MU_PIZERO;
    else if (interactionTypeString == "CCDIS_MU_P_PIZERO") return CCDIS_MU_P_PIZERO;
    else if (interactionTypeString == "CCDIS_MU_P_P_PIZERO") return CCDIS_MU_P_P_PIZERO;
    else if (interactionTypeString == "CCDIS_MU_P_P_P_PIZERO") return CCDIS_MU_P_P_P_PIZERO;
    else if (interactionTypeString == "CCDIS_MU_P_P_P_P_PIZERO") return CCDIS_MU_P_P_P_P_PIZERO;
    else if (interactionTypeString == "CCDIS_MU_P_P_P_P_P_PIZERO") return CCDIS_MU_P_P_P_P_P_PIZERO;
    else if (interactionTypeString == "NCDIS_P") return NCDIS_P;
    else if (interactionTypeString == "NCDIS_P_P") return NCDIS_P_P;
    else if (interactionTypeString == "NCDIS_P_P_P") return NCDIS_P_P_P;
    else if (interactionTypeString == "NCDIS_P_P_P_P") return NCDIS_P_P_P_P;
    else if (interactionTypeString == "NCDIS_P_P_P_P_P") return NCDIS_P_P_P_P_P;
    else if (interactionTypeString == "NCDIS_PIPLUS") return NCDIS_PIPLUS;
    else if (interactionTypeString == "NCDIS_P_PIPLUS") return NCDIS_P_PIPLUS;
    else if (interactionTypeString == "NCDIS_P_P_PIPLUS") return NCDIS_P_P_PIPLUS;
    else if (interactionTypeString == "NCDIS_P_P_P_PIPLUS") return NCDIS_P_P_P_PIPLUS;
    else if (interactionTypeString == "NCDIS_P_P_P_P_PIPLUS") return NCDIS_P_P_P_P_PIPLUS;
    else if (interactionTypeString == "NCDIS_P_P_P_P_P_PIPLUS") return NCDIS_P_P_P_P_P_PIPLUS;
    else if (interactionTypeString == "NCDIS_PIMINUS") return NCDIS_PIMINUS;
    else if (interactionTypeString == "NCDIS_P_PIMINUS") return NCDIS_P_PIMINUS;
    else if (interactionTypeString == "NCDIS_P_P_PIMINUS") return NCDIS_P_P_PIMINUS;
    else if (interactionTypeString == "NCDIS_P_P_P_PIMINUS") return NCDIS_P_P_P_PIMINUS;
    else if (interactionTypeString == "NCDIS_P_P_P_P_PIMINUS") return NCDIS_P_P_P_P_PIMINUS;
    else if (interactionTypeString == "NCDIS_P_P_P_P_P_PIMINUS") return NCDIS_P_P_P_P_P_PIMINUS;
    else if (interactionTypeString == "NCDIS_PHOTON") return NCDIS_PHOTON;
    else if (interactionTypeString == "NCDIS_P_PHOTON") return NCDIS_P_PHOTON;
    else if (interactionTypeString == "NCDIS_P_P_PHOTON") return NCDIS_P_P_PHOTON;
    else if (interactionTypeString == "NCDIS_P_P_P_PHOTON") return NCDIS_P_P_P_PHOTON;
    else if (interactionTypeString == "NCDIS_P_P_P_P_PHOTON") return NCDIS_P_P_P_P_PHOTON;
    else if (interactionTypeString == "NCDIS_P_P_P_P_P_PHOTON") return NCDIS_P_P_P_P_P_PHOTON;
    else if (interactionTypeString == "NCDIS_PIZERO") return NCDIS_PIZERO;
    else if (interactionTypeString == "NCDIS_P_PIZERO") return NCDIS_P_PIZERO;
    else if (interactionTypeString == "NCDIS_P_P_PIZERO") return NCDIS_P_P_PIZERO;
    else if (interactionTypeString == "NCDIS_P_P_P_PIZERO") return NCDIS_P_P_P_PIZERO;
    else if (interactionTypeString == "NCDIS_P_P_P_P_PIZERO") return NCDIS_P_P_P_P_PIZERO;
    else if (interactionTypeString == "NCDIS_P_P_P_P_P_PIZERO") return NCDIS_P_P_P_P_P_PIZERO;
    else if (interactionTypeString == "CCCOH") return CCCOH;
    else if (interactionTypeString == "NCCOH") return NCCOH;
    else if (interactionTypeString == "COSMIC_RAY_MU") return COSMIC_RAY_MU;
    else if (interactionTypeString == "COSMIC_RAY_P") return COSMIC_RAY_P;
    else if (interactionTypeString == "COSMIC_RAY_E") return COSMIC_RAY_E;
    else if (interactionTypeString == "COSMIC_RAY_PHOTON") return COSMIC_RAY_PHOTON;
    else if (interactionTypeString == "COSMIC_RAY_OTHER") return COSMIC_RAY_OTHER;
    else if (interactionTypeString == "BEAM_PARTICLE_MU") return BEAM_PARTICLE_MU;
    else if (interactionTypeString == "BEAM_PARTICLE_P") return BEAM_PARTICLE_P;
    else if (interactionTypeString == "BEAM_PARTICLE_E") return BEAM_PARTICLE_E;
    else if (interactionTypeString == "BEAM_PARTICLE_PHOTON") return BEAM_PARTICLE_PHOTON;
    else if (interactionTypeString == "BEAM_PARTICLE_PI_PLUS") return BEAM_PARTICLE_PI_PLUS;
    else if (interactionTypeString == "BEAM_PARTICLE_PI_MINUS") return BEAM_PARTICLE_PI_MINUS;
    else if (interactionTypeString == "BEAM_PARTICLE_KAON_PLUS") return BEAM_PARTICLE_KAON_PLUS;
    else if (interactionTypeString == "BEAM_PARTICLE_KAON_MINUS") return BEAM_PARTICLE_KAON_MINUS;
    else if (interactionTypeString == "BEAM_PARTICLE_OTHER") return BEAM_PARTICLE_OTHER;
    else if (interactionTypeString == "OTHER_INTERACTION") return OTHER_INTERACTION;
    else if (interactionTypeString == "ALL_INTERACTIONS") return ALL_INTERACTIONS;
    else if (interactionTypeString == "CCMEC_MU") return CCMEC_MU;
    else if (interactionTypeString == "CCMEC_MU_P") return CCMEC_MU_P;
    else if (interactionTypeString == "CCMEC_MU_P_P") return CCMEC_MU_P_P;
    else if (interactionTypeString == "CCMEC_MU_P_P_P") return CCMEC_MU_P_P_P;
    else if (interactionTypeString == "CCMEC_MU_P_P_P_P") return CCMEC_MU_P_P_P_P;
    else if (interactionTypeString == "CCMEC_MU_P_P_P_P_P") return CCMEC_MU_P_P_P_P_P;
    else if (interactionTypeString == "NCMEC_MU") return NCMEC_MU;
    else if (interactionTypeString == "NCMEC_MU_P") return NCMEC_MU_P;
    else if (interactionTypeString == "NCMEC_MU_P_P") return NCMEC_MU_P_P;
    else if (interactionTypeString == "NCMEC_MU_P_P_P") return NCMEC_MU_P_P_P;
    else if (interactionTypeString == "NCMEC_MU_P_P_P_P") return NCMEC_MU_P_P_P_P;
    else if (interactionTypeString == "NCMEC_MU_P_P_P_P_P") return NCMEC_MU_P_P_P_P_P;
    else if (interactionTypeString =="CCMEC_P") return CCMEC_P;
    else if (interactionTypeString =="NCMEC_P") return NCMEC_P;
    else if (interactionTypeString =="NCRES_MU") return NCRES_MU;
    else if (interactionTypeString =="NCRES_MU_P") return NCRES_MU_P;
    return OTHER_INTERACTION;
}

} // namespace lar_content
