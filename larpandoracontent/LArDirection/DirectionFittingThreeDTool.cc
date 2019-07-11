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
float globalTrackLength3D(0.f), globalEndpointRange(10.f);
static lar_content::DirectionFittingThreeDTool::LookupTable globalLookupTable3D;

//----------------------------------------------------------------------------------------------------------------------------------

#include "MinuitThreeDFunctions.h"

namespace lar_content
{

DirectionFittingThreeDTool::DirectionFittingThreeDTool() :
    m_endpointRange(10.f),
    m_slidingFitWindow(5),
    m_minClusterCaloHits(5),
    m_minClusterLength(1.f),
    m_tableInitialEnergy(250.f),
    m_tableStepSize(0.01f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DirectionFittingThreeDTool::~DirectionFittingThreeDTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DirectionFittingThreeDTool::DirectionFitObject DirectionFittingThreeDTool::GetPfoDirection(const pandora::ParticleFlowObject *const pPfo)
{
    try 
    {
        globalLookupTable3D.SetInitialEnergy(m_tableInitialEnergy);
        globalLookupTable3D.SetBinWidth(m_tableStepSize);
        FillLookupTable3D(globalLookupTable3D, 938.0, 0.0);
        globalEndpointRange = m_endpointRange;

        HitObjectVector hitObjectVector, filteredHitObjectVector, endpointHitObjectVector;

        this->FillHitObjectVector(pPfo, hitObjectVector);

        if (hitObjectVector.size() <= m_minClusterCaloHits || std::abs(hitObjectVector.back().GetLongitudinalPosition() - hitObjectVector.front().GetLongitudinalPosition()) < m_minClusterLength)
        {
            std::cout << "Not enough hits or PFO too short. Number of hits: " << hitObjectVector.size() << std::endl;
            throw STATUS_CODE_NOT_FOUND;
        }

        this->FilterHitObjectVector(hitObjectVector, filteredHitObjectVector);
        this->SelectEndpoint(filteredHitObjectVector, endpointHitObjectVector);

        globalTrackLength3D = std::abs(endpointHitObjectVector.back().GetLongitudinalPosition() - endpointHitObjectVector.front().GetLongitudinalPosition());

        DirectionFitObject finalDirectionFitObject;
        this->FitHitObjectVector(endpointHitObjectVector, finalDirectionFitObject);
        //this->SetMCInformation(pPfo, finalDirectionFitObject);

        globalTrackLength3D = 0.f;

        std::cout << "MC PDG: " << finalDirectionFitObject.GetMCParent() << std::endl;
        std::cout << "Contained: " << finalDirectionFitObject.GetContained() << std::endl;
        std::cout << "Fit mass: " << finalDirectionFitObject.GetFitMass() << std::endl;
        std::cout << "MinChiSquaredPerHit: " << finalDirectionFitObject.GetMinChiSquaredPerHit() << std::endl;
        std::cout << "---------------------" << std::endl;

        return finalDirectionFitObject;
    }

    catch (StatusCodeException &statusCodeException)
    {
        std::cout << "3D Fit Failure in tool." << std::endl;
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

    int numberWHits(0);

    for (const auto pCaloHit : caloHitList)
    {   
        if (!(static_cast<const CaloHit*>(pCaloHit->GetParentAddress())->GetHitType() == TPC_VIEW_W))
            continue;

        ++numberWHits;
    }

    const float wirePitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    int slidingFitWindowToUse(m_slidingFitWindow);

    for (int i = 1; i <= 10; ++i)
    {
        try 
        {
            const ThreeDSlidingFitResult testSlidingFit(clusterList.front(), slidingFitWindowToUse, wirePitch);
        }
        catch (...)
        {
            std::cout << "Can't fit sliding fit, retrying" << std::endl;
            slidingFitWindowToUse *= 2;
            continue;
        }

        break; 
    }

    const ThreeDSlidingFitResult slidingFit(clusterList.front(), slidingFitWindowToUse, wirePitch);

    for (const auto pCaloHit : caloHitList)
    {   
        if (numberWHits >= 5 && !(static_cast<const CaloHit*>(pCaloHit->GetParentAddress())->GetHitType() == TPC_VIEW_W))
            continue;

        const CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
        float caloHitEnergy(pCaloHit->GetInputEnergy());
        float hitWidth(pCaloHit->GetCellSize1());

        float rL(0.f);
        rL = slidingFit.GetLongitudinalDisplacement(caloHitPosition);

        CartesianVector localDirection(0.f, 0.f, 0.f), zAxis(0.f, 0.f, 1.f);
        slidingFit.GetGlobalFitDirection(rL, localDirection);
        float segmentLength(wirePitch/cos(localDirection.GetOpeningAngle(zAxis)));

        float dQdx(caloHitEnergy/segmentLength);
        dQdx *= GetDriftCoordinateCorrection(caloHitPosition);
        dQdx *= GetYZCoordinateCorrection(caloHitPosition);

        float dEdx(this->CalculateModBoxdEdx(dQdx));
    
        if (dEdx < 0.f)
        {
            std::cout << "ModBox failed." << std::endl;
            dEdx = 2.1;
            continue;
            //throw STATUS_CODE_NOT_FOUND; 
        }

        HitObject hitObject(pCaloHit, rL, caloHitEnergy, hitWidth, segmentLength, dQdx, dEdx);
        hitObjectVector.push_back(hitObject);
    }   

    std::sort(hitObjectVector.begin(), hitObjectVector.end(), SortByLongitudinalPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DirectionFittingThreeDTool::GetDriftCoordinateCorrection(const pandora::CartesianVector &positionVector)
{
    float x(positionVector.GetX());

    float parameter0(1.35970e-01);
    float parameter1(5.18310e-03);
    float parameter2(-3.45421e+00);
    float parameter3(1.37225e+00);
    float parameter4(-4.23881e-03);
    float parameter5(-6.83842e-01);

    return (parameter0 + parameter1 * (x - parameter2) + parameter3 * sin(parameter4 * x - parameter5));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DirectionFittingThreeDTool::GetYZCoordinateCorrection(const pandora::CartesianVector &positionVector)
{
    float y(positionVector.GetY());
    float z(positionVector.GetZ());

    if ((z > 0 && z < 400) && (y > (-120.0 + ((220.0/400.0) * z))) && (y < ((120.0/250.0) * z))) 
        return 1.3;
    else
        return 1.0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DirectionFittingThreeDTool::CalculateModBoxdEdx(float &dQdx)
{
    const float C_modbox(0.00507689);
    const float W_ion_modbox(23.6e-6);
    const float E_modbox(0.273);
    const float rho_modbox(1.38);
    const float beta_modbox(0.212);
    const float alpha_modbox(0.93);

    float exponential_component(exp((dQdx/C_modbox) * ((beta_modbox * W_ion_modbox)/(rho_modbox * E_modbox))));
    float dEdx((exponential_component - alpha_modbox)/(beta_modbox/(rho_modbox * E_modbox)));
    return dEdx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFittingThreeDTool::FilterHitObjectVector(HitObjectVector &hitObjectVector, HitObjectVector &filteredHitObjectVector)
{
    if (hitObjectVector.size() <= 5)
    {
        filteredHitObjectVector = hitObjectVector;  
        return;
    }

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
        Hits->SetPoint(n, hitObject.GetLongitudinalPosition(), hitObject.GetdEdx());
        ++n;
    }

    TCanvas *canvas = new TCanvas("canvas", "canvas", 900, 600);
    canvas->cd();

    Hits->GetXaxis()->SetLimits(hitObjectVector.front().GetLongitudinalPosition() - 0.05 * hitObjectVector.back().GetLongitudinalPosition(), hitObjectVector.back().GetLongitudinalPosition() + 0.05 * hitObjectVector.back().GetLongitudinalPosition());
    Hits->SetMarkerStyle(20);
    Hits->SetMarkerSize(0.5);
    Hits->SetMarkerColor(kBlack); 

    Hits->SetTitle(";Longitudinal Coordinate L_{3D} (cm); Hit dE/dx (MeV/cm)"); 
    Hits->Draw("AP");
    canvas->SaveAs("allhits_filtered.png");
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFittingThreeDTool::SelectEndpoint(HitObjectVector &hitObjectVector, HitObjectVector &filteredHitObjectVector)
{
    float leftTotaldEdx(0.f), rightTotaldEdx(0.f);
    HitObjectVector leftHitObjectVector, rightHitObjectVector;
    
    for (const auto &hitObject : hitObjectVector)
    {
        if (hitObject.GetLongitudinalPosition() - hitObjectVector.front().GetLongitudinalPosition() <= m_endpointRange)
        {
            leftHitObjectVector.push_back(hitObject);
            leftTotaldEdx += hitObject.GetdEdx();
        }

        if (hitObjectVector.back().GetLongitudinalPosition() - hitObject.GetLongitudinalPosition() <= m_endpointRange)
        {
            rightHitObjectVector.push_back(hitObject);
            rightTotaldEdx += hitObject.GetdEdx(); 
        }
    }

    HitObjectVector selectedHitObjectVector(leftTotaldEdx > rightTotaldEdx ? leftHitObjectVector : rightHitObjectVector);
    float longitudinalOffset(std::min(selectedHitObjectVector.front().GetLongitudinalPosition(), selectedHitObjectVector.back().GetLongitudinalPosition()));
    float offsetCorrection(0.5 * (std::ceil(selectedHitObjectVector.back().GetLongitudinalPosition()) - selectedHitObjectVector.back().GetLongitudinalPosition()));
    longitudinalOffset -= offsetCorrection;

    for (auto &hitObject : selectedHitObjectVector)
        hitObject.SetLongitudinalPosition(hitObject.GetLongitudinalPosition() - longitudinalOffset); 

    filteredHitObjectVector = selectedHitObjectVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionFittingThreeDTool::SetNearestNeighbourValues(HitObjectVector &innerHitObjectVector, int &nNeighboursToConsider)
{
    float trackLength(innerHitObjectVector.back().GetLongitudinalPosition());

    std::sort(innerHitObjectVector.begin(), innerHitObjectVector.end(), SortByEnergy);
    float EnergyRange(innerHitObjectVector.back().GetEnergy() - innerHitObjectVector.front().GetEnergy());

    for (HitObject &hitObject1 : innerHitObjectVector)
    {
        std::vector<float> distancesToNN;

        for (HitObject &hitObject2 : innerHitObjectVector)
        {
            if (&hitObject1 == &hitObject2)
                continue;

            float EnergyDistance((trackLength/EnergyRange) * (std::abs(hitObject1.GetEnergy() - hitObject2.GetEnergy())));
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
    int numberHits(hitObjectVector.size()), fitStatus(-1);
    float particleForwardsChiSquared(0.f), particleBackwardsChiSquared(0.f);
    FitParameters fitParameters;

    this->PerformFits(hitObjectVector, particleForwardsChiSquared, particleBackwardsChiSquared, fitParameters, fitStatus);

    DirectionFitObject finalDirectionFitObject(hitObjectVector, numberHits, particleForwardsChiSquared, particleBackwardsChiSquared, fitParameters, fitStatus);
    finalDirectionFitObject.SetFitParameters(fitParameters);

    fitResult = finalDirectionFitObject;
}

//----------------------------------------------------------------------------------------------------------------------------------

void DirectionFittingThreeDTool::PerformFits(HitObjectVector &hitObjectVector, float &forwardsChiSquared, float &backwardsChiSquared, FitParameters &fitParameters, int &fitStatus)
{
    LookupTable lookupTable = globalLookupTable3D;

    pMinuitVector3D->clear();
    pMinuitVector3D->insert(pMinuitVector3D->begin(), hitObjectVector.begin(), hitObjectVector.end());
    int fitStatus1(0), fitStatus2(0);

    //---------------------------------------------------------------------------------------------------
    //Forwards Fit

    int nParameters = 2;
    std::string parName[nParameters] = {"ENDENERGY", "EXTRADOF"};
    double vstart[nParameters] = {1.0, 1.0};
    double step[nParameters] = {0.001, 0.001};
    double lowphysbound[nParameters] = {0.01, 0.95};
    double highphysbound[nParameters] = {5.0, 1.05};

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

    /*
    //Refit in case of fit failure
    int refitCounter1(0);
    while (fitStatus1 != 0 && refitCounter1 < 3)
    {   
        vstart[1] = 15.0 + refitCounter1 * 250.0; //new start mass

        for (int j = 0 ; j < nParameters ; ++j)
            pMinuit->mnparm(j, parName[j].c_str(), vstart[j], step[j], lowphysbound[j], highphysbound[j], ierflg);

        pMinuit->mnexcm("MIGRAD", arglist, 1, fitStatus1);
        ++refitCounter1;
    }   
    */

    double outpar[2], err[2];

    for (int k = 0; k < nParameters; k++)
        pMinuit->GetParameter(k, outpar[k], err[k]);

    delete pMinuit;

    //---------------------------------------------------------------------------------
    //Backwards Fit

    int nParameters2 = 2;
    std::string parName2[nParameters2] = {"ENDENERGY", "EXTRADOF"};
    double vstart2[nParameters2] = {1.0, 1.0};
    double step2[nParameters2] = {0.001, 0.001};
    double lowphysbound2[nParameters2] = {0.01, 0.95};
    double highphysbound2[nParameters2] = {5.0, 1.05};

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
        
    /*

    //Refit in case of fit failure
    int refitCounter2(0);
    while (fitStatus2 != 0 && refitCounter2 < 3)
    {   
        vstart2[1] = 15.0 + refitCounter2 * 250.0; //new start mass

        for (int j = 0 ; j < nParameters2 ; ++j)
            pMinuit2->mnparm(j, parName2[j].c_str(), vstart2[j], step2[j], lowphysbound2[j], highphysbound2[j], ierflg2);

        pMinuit2->mnexcm("MIGRAD", arglist2, 1, fitStatus2);
        ++refitCounter2;
    }  
    */

    double outpar2[2], err2[2];

    for (int k = 0; k < nParameters2; k++)
        pMinuit2->GetParameter(k, outpar2[k], err2[k]);

    delete pMinuit2;

    //--------------------------------------------------------------------------

    this->SetFinalFitValues(hitObjectVector, forwardsChiSquared, backwardsChiSquared, outpar, true);
    this->SetFinalFitValues(hitObjectVector, forwardsChiSquared, backwardsChiSquared, outpar2, false);

    //--------------------------------------------------------------------------

    if (forwardsChiSquared <= backwardsChiSquared)
    {
        float parameterZero(outpar[0]), parameterOne(outpar[1]), parameterTwo(0);
        FitParameters bestFitParameters(parameterZero, parameterOne, parameterTwo);
        fitParameters = bestFitParameters;
        fitStatus = fitStatus1;
    }
    else
    {
        float parameterZero(outpar2[0]), parameterOne(outpar2[1]), parameterTwo(0);
        FitParameters bestFitParameters(parameterZero, parameterOne, parameterTwo);
        fitParameters = bestFitParameters;
        fitStatus = fitStatus2;
    }

    pMinuitVector3D->clear();
}

//---------------------------------------------------------------------------------------------------------------------------------

void DirectionFittingThreeDTool::SetFinalFitValues(HitObjectVector &hitObjectVector, float &forwardsChiSquared, float &backwardsChiSquared, double (&fitParameters)[2], bool forwards)
{
    double L_offset(fitParameters[0]);
    double M = 938.0;

    FillLookupTable3D(globalLookupTable3D, 938.0, L_offset); 

    for (auto &hitObject : hitObjectVector)
    {
        if (forwards)
        {
            try
            {
                int binNumber(std::floor(hitObject.GetLongitudinalPosition()/globalLookupTable3D.GetBinWidth()));
                double fitdEdx(fitParameters[1] * BetheBloch3D(globalLookupTable3D.GetMap().at(binNumber), M));
                double hitUncertainty(1.0);

                hitObject.SetForwardsFitdEdx(fitdEdx); 
                forwardsChiSquared += ((hitObject.GetdEdx() - fitdEdx) * (hitObject.GetdEdx() - fitdEdx) )/(hitUncertainty * hitUncertainty);
            }
            catch (...)
            {
                continue;
            }
        }
        else
        {
            try
            {
                int binNumber(std::floor((globalTrackLength3D - hitObject.GetLongitudinalPosition())/globalLookupTable3D.GetBinWidth()));
                double fitdEdx(fitParameters[1] * BetheBloch3D(globalLookupTable3D.GetMap().at(binNumber), M));
                double hitUncertainty(1.0);

                hitObject.SetBackwardsFitdEdx(fitdEdx); 
                backwardsChiSquared += ((hitObject.GetdEdx() - fitdEdx) * (hitObject.GetdEdx() - fitdEdx) )/(hitUncertainty * hitUncertainty);
            }
            catch (...)
            {
                continue;
            }
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
    return hitEnergy1.GetEnergy() < hitEnergy2.GetEnergy();
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
