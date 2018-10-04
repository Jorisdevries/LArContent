/**
 *  @file   larpandoracontent/LArVertex/DirectionFittingThreeDTool.cc
 *
 *  @brief  Implementation of the candidate vertex creation Tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "DirectionFittingThreeDTool.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include <ctime>

using namespace pandora;

//----------------------------------------------------------------------------------------------------------------------------------

//nasty global parameters necessary for TMinuit
lar_content::DirectionFittingThreeDTool::HitObjectVector* pMinuitVector3D = new lar_content::DirectionFittingThreeDTool::HitObjectVector;
float globalTrackLength3D(0.f), globalEndpointRange(10.f), globalMass(100.0);
bool globalFixedMass(false);
static lar_content::DirectionFittingThreeDTool::LookupTable globalLookupTable3D;

//----------------------------------------------------------------------------------------------------------------------------------

#include "MinuitThreeDFunctions.h"

namespace lar_content
{

DirectionFittingThreeDTool::DirectionFittingThreeDTool() :
    m_endpointRange(10.f),
    m_slidingFitWindow(5),
    m_minClusterCaloHits(20),
    m_minClusterLength(5.f),
    m_tableInitialEnergy(250.f),
    m_tableStepSize(0.01f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DirectionFittingThreeDTool::~DirectionFittingThreeDTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DirectionFittingThreeDTool::DirectionFitObject DirectionFittingThreeDTool::GetPfoDirection(const pandora::ParticleFlowObject *const pPfo, float particleMass)
{
    try 
    {
        globalLookupTable3D.SetInitialEnergy(m_tableInitialEnergy);
        globalLookupTable3D.SetBinWidth(m_tableStepSize);
        FillLookupTable3D(globalLookupTable3D, particleMass, 0.0);

        globalEndpointRange = m_endpointRange;

        if (particleMass != 100.0)
        {
            globalFixedMass = true;
            globalMass = particleMass;
        }

        HitObjectVector hitObjectVector, filteredHitObjectVector, endpointHitObjectVector;

        this->FillHitObjectVector(pPfo, hitObjectVector);
        this->FilterHitObjectVector(hitObjectVector, filteredHitObjectVector);
        this->SelectEndpoint(filteredHitObjectVector, endpointHitObjectVector);

        DirectionFitObject finalDirectionFitObject;
        this->FitHitObjectVector(endpointHitObjectVector, finalDirectionFitObject);
        this->SetMCInformation(pPfo, finalDirectionFitObject);

        std::cout << "------------------------------" << std::endl;
        std::cout << "MC PDG: " << finalDirectionFitObject.GetMCParent() << std::endl;
        std::cout << "Contained: " << finalDirectionFitObject.GetContained() << std::endl;
        std::cout << "Fit mass: " << finalDirectionFitObject.GetFitMass() << std::endl;
        std::cout << "Fit status: " << finalDirectionFitObject.GetFitStatus() << std::endl;
        std::cout << "MinChiSquaredPerHit: " << finalDirectionFitObject.GetMinChiSquaredPerHit() << std::endl;
        std::cout << "Number hits: " << finalDirectionFitObject.GetNHits() << std::endl;
        std::cout << "------------------------------" << std::endl;

        return finalDirectionFitObject;
    }

    catch (StatusCodeException &statusCodeException)
    {
        std::cout << "Failure." << std::endl;
        throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFittingThreeDTool::FillHitObjectVector(const pandora::ParticleFlowObject *const pPfo, HitObjectVector &hitObjectVector)
{
    ClusterList clusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_3D, clusterList);
    ClusterVector clusterVector(clusterList.begin(), clusterList.end());

    std::sort(clusterVector.begin(), clusterVector.end(), [](const Cluster* const pCluster1, const Cluster* const pCluster2) -> bool { return LArClusterHelper::GetLength(pCluster1) > LArClusterHelper::GetLength(pCluster2); });

    OrderedCaloHitList orderedCaloHitList(clusterVector.front()->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);

    const float wirePitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const ThreeDSlidingFitResult slidingFit(clusterList.front(), m_slidingFitWindow, wirePitch);

    float segmentLength(0.3);

    for (const auto pCaloHit : caloHitList)
    {
        if (!(static_cast<const CaloHit*>(pCaloHit->GetParentAddress())->GetHitType() == TPC_VIEW_W))
            continue;

        const CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
        float caloHitEnergy(pCaloHit->GetInputEnergy());
        float hitWidth(pCaloHit->GetCellSize1());

        caloHitEnergy *= 197.0;        //ADC to electron (gain)
        caloHitEnergy *= 23.6/1e6;     //Ionisation energy per electron in MeV
        caloHitEnergy /= 0.62;         //Recombination

        float rL(0.f);
        rL = slidingFit.GetLongitudinalDisplacement(caloHitPosition);

        HitObject hitObject(pCaloHit, rL, caloHitEnergy, hitWidth, segmentLength);
        hitObjectVector.push_back(hitObject);
    }

    std::sort(hitObjectVector.begin(), hitObjectVector.end(), SortByLongitudinalPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFittingThreeDTool::FilterHitObjectVector(HitObjectVector &hitObjectVector, HitObjectVector &filteredHitObjectVector)
{
    float endpointProtectionRange(0.00);

    filteredHitObjectVector.insert(filteredHitObjectVector.begin(), hitObjectVector.begin(),  hitObjectVector.begin() + endpointProtectionRange * hitObjectVector.size());
    filteredHitObjectVector.insert(filteredHitObjectVector.begin(), hitObjectVector.begin() + (1.0 - endpointProtectionRange) * hitObjectVector.size(), hitObjectVector.end());

    HitObjectVector innerHitObjectVector(hitObjectVector.begin() + endpointProtectionRange * hitObjectVector.size(), hitObjectVector.begin() + (1.0 - endpointProtectionRange) * hitObjectVector.size());

    int nNeighboursToConsider(5);
    this->SetNearestNeighbourValues(innerHitObjectVector, nNeighboursToConsider);

     std::sort(innerHitObjectVector.begin(), innerHitObjectVector.end(), SortByDistanceToNN);
     filteredHitObjectVector.insert(filteredHitObjectVector.begin(), innerHitObjectVector.begin(), innerHitObjectVector.begin() + 0.72 * innerHitObjectVector.size()); //percentage optimised by tests

     std::sort(filteredHitObjectVector.begin(), filteredHitObjectVector.end(), SortByLongitudinalPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFittingThreeDTool::DrawHitObjectVector(HitObjectVector &hitObjectVector)
{
    TGraphErrors *Hits = new TGraphErrors(hitObjectVector.size());

    int n(0);
    
    for (HitObject hitObject : hitObjectVector)
    {    
        Hits->SetPoint(n, hitObject.GetLongitudinalPosition(), hitObject.GetHitEnergy());
        ++n;
    }

    TCanvas *canvas = new TCanvas("canvas", "canvas", 900, 600);
    canvas->cd();

    Hits->GetXaxis()->SetLimits(hitObjectVector.front().GetLongitudinalPosition() - 0.05 * hitObjectVector.back().GetLongitudinalPosition(), hitObjectVector.back().GetLongitudinalPosition() + 0.05 * hitObjectVector.back().GetLongitudinalPosition());
    Hits->SetMarkerStyle(20);
    Hits->SetMarkerSize(0.5);
    Hits->SetMarkerColor(kBlack); 

    Hits->SetTitle(";Longitudinal Coordinate L_{3D} (cm); Hit Energy (MeV)"); //Bethe-Bloch Theory Fit
    Hits->Draw("AP");
    canvas->SaveAs("allhits_filtered.png");
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFittingThreeDTool::SelectEndpoint(HitObjectVector &hitObjectVector, HitObjectVector &filteredHitObjectVector)
{
    const float wirePitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    float leftTotalEnergy(0.f), rightTotalEnergy(0.f);
    float leftTotalHitWidth(0.f), rightTotalHitWidth(0.f);
    HitObjectVector leftHitObjectVector, rightHitObjectVector;
    
    for (const auto &hitObject : hitObjectVector)
    {
        if (hitObject.GetLongitudinalPosition() - hitObjectVector.front().GetLongitudinalPosition() <= m_endpointRange - wirePitch)
        {
            leftHitObjectVector.push_back(hitObject);
            leftTotalEnergy += hitObject.GetHitEnergy();
            leftTotalHitWidth += hitObject.GetHitWidth();
        }

        if (hitObjectVector.back().GetLongitudinalPosition() - hitObject.GetLongitudinalPosition() <= m_endpointRange - wirePitch)
        {
            rightHitObjectVector.push_back(hitObject);
            rightTotalEnergy += hitObject.GetHitEnergy(); 
            rightTotalHitWidth += hitObject.GetHitWidth(); 
        }
    }

    HitObjectVector selectedHitObjectVector(leftTotalEnergy > rightTotalEnergy ? leftHitObjectVector : rightHitObjectVector);
    float totalHitWidth(leftTotalEnergy > rightTotalEnergy ? leftTotalHitWidth : rightTotalHitWidth);
    float endpointLength(wirePitch + std::abs(selectedHitObjectVector.back().GetLongitudinalPosition() - selectedHitObjectVector.front().GetLongitudinalPosition()));

    for (const auto &hitObject : selectedHitObjectVector)
    {
        float longitudinalPosition(((m_endpointRange - endpointLength)/2) + hitObject.GetLongitudinalPosition() - std::min(selectedHitObjectVector.back().GetLongitudinalPosition(), selectedHitObjectVector.front().GetLongitudinalPosition()));

        float hitEnergy(hitObject.GetHitEnergy());
        float hitWidth(hitObject.GetHitWidth());
        float segmentLength((hitWidth/totalHitWidth) * endpointLength);

        HitObject newHitObject(hitObject.GetCaloHit(), longitudinalPosition, hitEnergy, hitWidth, segmentLength);
        filteredHitObjectVector.push_back(newHitObject);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFittingThreeDTool::SetNearestNeighbourValues(HitObjectVector &innerHitObjectVector, int &nNeighboursToConsider)
{
    float trackLength(innerHitObjectVector.back().GetLongitudinalPosition());

    std::sort(innerHitObjectVector.begin(), innerHitObjectVector.end(), SortByEnergy);
    float EnergyRange(innerHitObjectVector.back().GetHitEnergy() - innerHitObjectVector.front().GetHitEnergy());

    for (HitObject &hitObject1 : innerHitObjectVector)
    {
        std::vector<float> distancesToNN;

        for (HitObject &hitObject2 : innerHitObjectVector)
        {
            if (&hitObject1 == &hitObject2)
                continue;

            float EnergyDistance((trackLength/EnergyRange) * (std::abs(hitObject1.GetHitEnergy() - hitObject2.GetHitEnergy())));
            float Ldistance(std::abs(hitObject1.GetLongitudinalPosition() - hitObject2.GetLongitudinalPosition()));
            float distanceToNN(std::sqrt(EnergyDistance*EnergyDistance + Ldistance*Ldistance));

            distancesToNN.push_back(distanceToNN);
        }

        std::sort(distancesToNN.begin(), distancesToNN.end());
        float nearestNeighboursDistanceSum(std::accumulate(distancesToNN.begin(), distancesToNN.begin() + nNeighboursToConsider, 0.f));
        hitObject1.SetDistanceToNN(nearestNeighboursDistanceSum);
    }

    std::sort(innerHitObjectVector.begin(), innerHitObjectVector.end(), SortByLongitudinalPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFittingThreeDTool::FitHitObjectVector(HitObjectVector &hitObjectVector, DirectionFittingThreeDTool::DirectionFitObject &fitResult)
{
    LookupTable lookupTable = globalLookupTable3D;
    globalTrackLength3D = std::abs(hitObjectVector.back().GetLongitudinalPosition() - hitObjectVector.front().GetLongitudinalPosition());

    pMinuitVector3D->clear();
    pMinuitVector3D->insert(pMinuitVector3D->begin(), hitObjectVector.begin(), hitObjectVector.end());
    int fitStatus1, fitStatus2;

    //---------------------------------------------------------------------------------------------------
    //Forwards Fit

    int nParameters = 2;
    std::string parName[nParameters] = {"ENDENERGY", "MASS"};
    double vstart[nParameters] = {0.5, 1100.0};
    double step[nParameters] = {0.1, 1.0};
    double lowphysbound[nParameters] = {0.1, 900.0};
    double highphysbound[nParameters] = {1.0, 1300.0};

    if (globalFixedMass)
        nParameters = 1;

    int ierflg(0);

    TMinuit *pMinuit = new TMinuit(nParameters);
    pMinuit->SetPrintLevel(-1);
    pMinuit->SetFCN(GetForwardsChiSquared3D);

    for (int j = 0 ; j < nParameters ; ++j)
        pMinuit->mnparm(j, parName[j].c_str(), vstart[j], step[j], lowphysbound[j], highphysbound[j], ierflg);

    double arglist[2];
    arglist[0] = 1000;
    arglist[1] = 1;
    pMinuit->mnexcm("MIGRAD", arglist, 1, fitStatus1);

    double outpar[2], err[2];

    for (int k = 0; k < nParameters; k++)
        pMinuit->GetParameter(k, outpar[k], err[k]);

    delete pMinuit;

    //---------------------------------------------------------------------------------
    //Backwards Fit

    int nParameters2 = 2;
    std::string parName2[nParameters2] = {"ENDENERGY", "MASS"};
    double vstart2[nParameters2] = {0.5, 1100.0};
    double step2[nParameters2] = {0.1, 1.0};
    double lowphysbound2[nParameters2] = {0.1, 900.0};
    double highphysbound2[nParameters2] = {1.0, 1300.0};

    if (globalFixedMass)
        nParameters2 = 1;

    int ierflg2(0);

    TMinuit *pMinuit2 = new TMinuit(nParameters2);
    pMinuit2->SetPrintLevel(-1);
    pMinuit2->SetFCN(GetBackwardsChiSquared3D);

    for (int j = 0 ; j < nParameters2 ; ++j)
        pMinuit2->mnparm(j, parName2[j].c_str(), vstart2[j], step2[j], lowphysbound2[j], highphysbound2[j], ierflg2);

    double arglist2[2];
    arglist2[0] = 1000;
    arglist2[1] = 1;
    pMinuit2->mnexcm("MIGRAD", arglist2, 1, fitStatus2);

    double outpar2[2], err2[2];

    for (int k = 0; k < nParameters2; k++)
        pMinuit2->GetParameter(k, outpar2[k], err2[k]);

    delete pMinuit2;

    if (globalFixedMass)
    {
        outpar[1] = globalMass;
        outpar2[1] = globalMass;
    }

    //--------------------------------------------------------------------------

    float forwardsChiSquared(0.f), backwardsChiSquared(0.f);

    this->SetFinalFitValues(hitObjectVector, forwardsChiSquared, backwardsChiSquared, outpar, true);
    this->SetFinalFitValues(hitObjectVector, forwardsChiSquared, backwardsChiSquared, outpar2, false);

    fitResult.SetNHits(hitObjectVector.size());
    fitResult.SetForwardsChiSquared(forwardsChiSquared);
    fitResult.SetBackwardsChiSquared(backwardsChiSquared);

    //--------------------------------------------------------------------------

    if (forwardsChiSquared <= backwardsChiSquared)
    {
        float parameterZero(outpar[0]), parameterOne(outpar[1]), parameterTwo(0);
        FitParameters bestFitParameters(parameterZero, parameterOne, parameterTwo);

        fitResult.SetFitParameters(bestFitParameters);
        fitResult.SetFitStatus(fitStatus1);
        fitResult.SetFitMass(parameterOne);
    }
    else
    {
        float parameterZero(outpar2[0]), parameterOne(outpar2[1]), parameterTwo(0);
        FitParameters bestFitParameters(parameterZero, parameterOne, parameterTwo);

        fitResult.SetFitParameters(bestFitParameters);
        fitResult.SetFitStatus(fitStatus2);
        fitResult.SetFitMass(parameterOne);
    }

    globalTrackLength3D = 0.f;
    pMinuitVector3D->clear();
}

//---------------------------------------------------------------------------------------------------------------------------------

void DirectionFittingThreeDTool::SetFinalFitValues(HitObjectVector &hitObjectVector, float &forwardsChiSquared, float &backwardsChiSquared, double (&fitParameters)[2], bool forwards)
{
    double L_offset(fitParameters[0]);
    double M = fitParameters[1];

    FillLookupTable3D(globalLookupTable3D, M, L_offset); 

    for (auto &hitObject : hitObjectVector)
    {
        if (forwards)
        {
            double fitEnergy = GetForwardsFitEnergy(hitObject);
            double hitUncertainty(1.0);

            hitObject.SetForwardsFitEnergy(fitEnergy); 
            forwardsChiSquared += ((hitObject.GetHitEnergy() - fitEnergy) * (hitObject.GetHitEnergy() - fitEnergy) )/(hitUncertainty * hitUncertainty);
        }
        else
        {
            double fitEnergy = GetBackwardsFitEnergy(hitObject);
            double hitUncertainty(1.0);

            hitObject.SetBackwardsFitEnergy(fitEnergy); 
            backwardsChiSquared += ((hitObject.GetHitEnergy() - fitEnergy) * (hitObject.GetHitEnergy() - fitEnergy) )/(hitUncertainty * hitUncertainty);
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------

void DirectionFittingThreeDTool::SetMCInformation(const pandora::ParticleFlowObject *const pPfo, DirectionFittingThreeDTool::DirectionFitObject &fitResult)
{
    const pandora::MCParticle* pMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));
    fitResult.SetMCParent(pMCParticle->GetParticleId());
    fitResult.SetContained(this->IsParticleContained(pMCParticle));
}

//----------------------------------------------------------------------------------------------------------------------------------

bool DirectionFittingThreeDTool::IsParticleContained(const MCParticle* pMCParticle)
{
    const CartesianVector mcVertex(pMCParticle->GetVertex());
    const CartesianVector mcEndpoint(pMCParticle->GetEndpoint());

    const float eVx(256.35), eVy(233.), eVz(1036.8);
    const float xBorder(10.), yBorder(20.), zBorder(10.);

    if ((mcEndpoint.GetX() < (eVx - xBorder)) && (mcEndpoint.GetX() > xBorder) && (mcEndpoint.GetY() < (eVy / 2. - yBorder)) && (mcEndpoint.GetY() > (-eVy / 2. + yBorder)) && (mcEndpoint.GetZ() < (eVz - zBorder)) && (mcEndpoint.GetZ() > zBorder))
    {   
        if (!LArGeometryHelper::IsInGap(this->GetPandora(), mcVertex, TPC_VIEW_W, 0.05*((mcEndpoint - mcVertex).GetMagnitude())) && !LArGeometryHelper::IsInGap(this->GetPandora(), mcEndpoint, TPC_VIEW_W, 0.05*((mcEndpoint - mcVertex).GetMagnitude())))
            return true;
    }   

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionFittingThreeDTool::SortByLongitudinalPosition(HitObject &hitObject1, HitObject &hitObject2)
{
    return hitObject1.GetLongitudinalPosition() < hitObject2.GetLongitudinalPosition();
}

//----------------------------------------------------------------------------------------------------------------------------------

bool DirectionFittingThreeDTool::SortByDistanceToNN(HitObject &hitObject1, HitObject &hitObject2)
{
    return hitObject1.GetDistanceToNN() < hitObject2.GetDistanceToNN();
}

//----------------------------------------------------------------------------------------------------------------------------------

bool DirectionFittingThreeDTool::SortByEnergy(HitObject &hitEnergy1, HitObject &hitEnergy2)
{
    return hitEnergy1.GetHitEnergy() < hitEnergy2.GetHitEnergy();
}

//----------------------------------------------------------------------------------------------------------------------------------

StatusCode DirectionFittingThreeDTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    return STATUS_CODE_SUCCESS;
}

//----------------------------------------------------------------------------------------------------------------------------------

} // namespac lar_content
