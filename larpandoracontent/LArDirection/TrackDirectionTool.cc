/**
 *  @file   larpandoracontent/LArVertex/TrackDirectionTool.cc
 *
 *  @brief  Implementation of the candidate vertex creation Tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "TrackDirectionTool.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include <ctime>
#include <math.h>

using namespace pandora;

//----------------------------------------------------------------------------------------------------------------------------------

//nasty global parameters necessary for TMinuit
lar_content::TrackDirectionTool::HitChargeVector* pMinuitVector = new lar_content::TrackDirectionTool::HitChargeVector;

float globalTotalCharge(0.f), globalTrackLength(0.f), globalTotalHitWidth(0.f);
int globalTrueDirection(-1);

static lar_content::TrackDirectionTool::LookupTable globalLookupTable;

//----------------------------------------------------------------------------------------------------------------------------------

#include "ToolMinuitFunctions.h"

namespace lar_content
{

TrackDirectionTool::TrackDirectionTool() :
    m_slidingFitWindow(5),
    m_minClusterCaloHits(25),
    m_minClusterLength(5.f),
    m_targetParticlePDG(13),
    m_numberTrackEndHits(100000),
    m_endpointProtectionFraction(0.05),
    m_enableFragmentRemoval(true),
    m_enableBraggPeakFilter(true),
    m_enableSplitting(true),
    m_tableInitialEnergy(2000.f),
    m_tableStepSize(0.5f),
    m_writeTable(false),
    m_useMCInformation(false),
    m_lookupTableFileName("lookuptable.root"),
    m_probabilityFileName("probability.root"),
    m_treeName("lookuptable")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackDirectionTool::~TrackDirectionTool()
{
    if (m_writeTable)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "FilterTree", m_fileName.c_str(), "UPDATE"));
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "FilterHitTree", m_fileName.c_str(), "UPDATE"));
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "BraggPeakFilterTree", m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackDirectionTool::DirectionFitObject TrackDirectionTool::GetClusterDirection(const Cluster *const pTargetClusterW)
{
    try
    {
        if (pTargetClusterW->GetNCaloHits() < m_minClusterCaloHits || pTargetClusterW->GetOrderedCaloHitList().size() > 5000)
        {
            //std::cout << "Direction fit error: invalid cluster" << std::endl;
            throw STATUS_CODE_FAILURE;
        }

        //if (globalLookupTable.GetMap().empty())
            this->SetLookupTable(105.7);

        DirectionFitObject finalDirectionFitObject;

        this->AddToSlidingFitCache(pTargetClusterW);
        this->GetCalorimetricDirection(pTargetClusterW, finalDirectionFitObject);
        this->ComputeProbability(finalDirectionFitObject);
        this->SetEndpoints(finalDirectionFitObject, pTargetClusterW);
        //this->SetMCTruth(finalDirectionFitObject, pTargetClusterW);

        if (m_useMCInformation)
        {
            const auto pMCParticle(MCParticleHelper::GetMainMCParticle(pTargetClusterW));
            globalTrueDirection = (pMCParticle->GetVertex().GetZ() < pMCParticle->GetEndpoint().GetZ() ? 1 : 0); 
        }

        this->TidyUp();
        return finalDirectionFitObject;
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();
        throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackDirectionTool::DirectionFitObject TrackDirectionTool::GetPfoDirection(const pandora::ParticleFlowObject *const pPfo)
{
    try 
    {
        const pandora::Vertex *const pVertex = LArPfoHelper::GetVertex(pPfo);
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        LArTrackStateVector trackStateVector;
        LArPfoHelper::GetSlidingFitTrajectory(pPfo, pVertex, m_slidingFitWindow, slidingFitPitch, trackStateVector);

        if (trackStateVector.size() == 0)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        const Cluster *const pClusterW = GetTargetClusterFromPFO(pPfo, trackStateVector);

        DirectionFitObject finalDirectionFitObject = GetClusterDirection(pClusterW);
        this->ComputeProbability(finalDirectionFitObject);

        //If the PFO is 3D, then 3D endpoints should be set 
        if (LArPfoHelper::IsThreeD(pPfo))
            SetEndpoints(finalDirectionFitObject, trackStateVector);

        this->TidyUp();
        return finalDirectionFitObject;
    }

    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();
        throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackDirectionTool::DirectionFitObject TrackDirectionTool::GetClusterDirection(const Cluster *const pTargetClusterW, float massHypothesis)
{
    try
    {
        if (pTargetClusterW->GetNCaloHits() < m_minClusterCaloHits || pTargetClusterW->GetOrderedCaloHitList().size() > 5000)
        {
            //std::cout << "Direction fit error: invalid cluster" << std::endl;
            throw STATUS_CODE_FAILURE;
        }

        //if (globalLookupTable.GetMap().empty())
            this->SetLookupTable(massHypothesis);

        DirectionFitObject finalDirectionFitObject;

        this->AddToSlidingFitCache(pTargetClusterW);
        this->GetCalorimetricDirection(pTargetClusterW, finalDirectionFitObject);
        this->ComputeProbability(finalDirectionFitObject);
        this->SetEndpoints(finalDirectionFitObject, pTargetClusterW);
        //this->SetMCTruth(finalDirectionFitObject, pTargetClusterW);

        this->TidyUp();
        return finalDirectionFitObject;
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();
        throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackDirectionTool::DirectionFitObject TrackDirectionTool::GetPfoDirection(const pandora::ParticleFlowObject *const pPfo, float massHypothesis)
{
    try 
    {
        const pandora::Vertex *const pVertex = LArPfoHelper::GetVertex(pPfo);
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        LArTrackStateVector trackStateVector;
        LArPfoHelper::GetSlidingFitTrajectory(pPfo, pVertex, m_slidingFitWindow, slidingFitPitch, trackStateVector);

        if (trackStateVector.size() == 0)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        const Cluster *const pClusterW = GetTargetClusterFromPFO(pPfo, trackStateVector);

        DirectionFitObject finalDirectionFitObject = GetClusterDirection(pClusterW, massHypothesis);
        //this->ComputeProbability(finalDirectionFitObject);

        //If the PFO is 3D, then 3D endpoints should be set 
        if (LArPfoHelper::IsThreeD(pPfo))
            SetEndpoints(finalDirectionFitObject, trackStateVector);

        this->TidyUp();
        return finalDirectionFitObject;
    }

    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();
        throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::WriteLookupTableToTree(LookupTable &lookupTable)
{
    std::vector<int> mapVector1, reverseMapVector2;
    std::vector<double> mapVector2, reverseMapVector1;

    for (auto &element : lookupTable.GetMap())
    {
        mapVector1.push_back(element.first);
        mapVector2.push_back(element.second);
    }

    for (auto &element : lookupTable.GetReverseMap())
    {
        reverseMapVector1.push_back(element.first);
        reverseMapVector2.push_back(element.second);
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mapVector1", &mapVector1));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mapVector2", &mapVector2));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "reverseMapVector1", &reverseMapVector1));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "reverseMapVector2", &reverseMapVector2));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "binWidth", lookupTable.GetBinWidth()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "initialEnergy", lookupTable.GetInitialEnergy()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "maxRange", lookupTable.GetMaxRange()));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SetLookupTable(float particleMass)
{
    globalLookupTable.SetInitialEnergy(m_tableInitialEnergy);
    globalLookupTable.SetBinWidth(m_tableStepSize);

    ifstream inputFile(m_lookupTableFileName);

    if (inputFile) 
    {
        this->ReadLookupTableFromTree(globalLookupTable);
    
        if (m_writeTable)
        {
            FillLookupTable(globalLookupTable, particleMass);
            this->WriteLookupTableToTree(globalLookupTable);
        }
    }
    else
    {
        //std::cout << "WARNING: filling lookup table because lookuptable.root was not found. To create it, include <WriteTable>true</WriteTable> to the Pandora settings XML file." << std::endl;
        //
        FillLookupTable(globalLookupTable, particleMass);

        if (m_writeTable)
            this->WriteLookupTableToTree(globalLookupTable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster* TrackDirectionTool::GetTargetClusterFromPFO(const ParticleFlowObject* pPfo, const LArTrackStateVector &trackStateVector)
{
    //HitType hitType(TPC_VIEW_W);
    ClusterList clusterListW;
    LArPfoHelper::GetTwoDClusterList(pPfo, clusterListW);

    TrackState firstTrackState(*(trackStateVector.begin())), lastTrackState(trackStateVector.back());
    const pandora::CartesianVector initialPosition(firstTrackState.GetPosition());
    const pandora::CartesianVector endPosition(lastTrackState.GetPosition());

    const pandora::CartesianVector lowZVector(initialPosition.GetZ() < endPosition.GetZ() ? initialPosition : endPosition);
    const pandora::CartesianVector highZVector(initialPosition.GetZ() > endPosition.GetZ() ? initialPosition : endPosition);

    float currentEndpointDistance(std::numeric_limits<float>::max());
    ClusterList bestClusterList;

    for (const Cluster *const pCluster : clusterListW)
    {    
        CartesianVector innerCoordinate(0.f, 0.f, 0.f), outerCoordinate(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinates(pCluster, innerCoordinate, outerCoordinate);

        const pandora::CartesianVector lowZClusterVector(innerCoordinate.GetZ() < outerCoordinate.GetZ() ? innerCoordinate : outerCoordinate);
        const pandora::CartesianVector highZClusterVector(innerCoordinate.GetZ() > outerCoordinate.GetZ() ? innerCoordinate : outerCoordinate);

        /*
        std::cout << "Cluster inner coordinates: (" << innerCoordinate.GetX() << ", " << innerCoordinate.GetY() << ", " << innerCoordinate.GetZ() << ")" << std::endl;
        std::cout << "Cluster outer coordinates: (" << outerCoordinate.GetX() << ", " << outerCoordinate.GetY() << ", " << outerCoordinate.GetZ() << ")" << std::endl;
        std::cout << "Track low Z coordinates: (" << lowZVector.GetX() << ", " << lowZVector.GetY() << ", " << lowZVector.GetZ() << ")" << std::endl;
        std::cout << "Track high Z coordinates: (" << highZVector.GetX() << ", " << highZVector.GetY() << ", " << highZVector.GetZ() << ")" << std::endl;
        */

        if (innerCoordinate.GetY() != 0 || outerCoordinate.GetY() != 0) 
            continue;

        float endpointDistance(std::abs(lowZVector.GetZ() - lowZClusterVector.GetZ()));

        if (endpointDistance < currentEndpointDistance)
        {    
            currentEndpointDistance = endpointDistance; 
            bestClusterList.clear();
            bestClusterList.push_back(pCluster);
        }    
    } 

    if (bestClusterList.size() == 0)
    {
        std::cout << "ERROR: no W clusters could be extracted from the PFO!" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }

    const Cluster *const pCluster(*(bestClusterList.begin())); 
    return pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::ReadLookupTableFromTree(LookupTable &lookupTable)
{
    TFile *f = TFile::Open(m_lookupTableFileName.c_str());
    TTree *t1 = (TTree*)f->Get(m_treeName.c_str());

    std::vector<int> mapVector1, reverseMapVector2;
    std::vector<double> mapVector2, reverseMapVector1;

    std::vector<int> *pMapVector1 = 0;
    std::vector<int> *pReverseMapVector2 = 0;
    std::vector<double> *pMapVector2 = 0;
    std::vector<double> *pReverseMapVector1 = 0;

    TBranch *pBranchMapVector1 = 0;
    TBranch *pBranchMapVector2 = 0;
    TBranch *pBranchReverseMapVector1 = 0;
    TBranch *pBranchReverseMapVector2 = 0;

    t1->SetBranchAddress("mapVector1", &pMapVector1, &pBranchMapVector1);
    t1->SetBranchAddress("mapVector2", &pMapVector2, &pBranchMapVector2);
    t1->SetBranchAddress("reverseMapVector1", &pReverseMapVector1, &pBranchReverseMapVector1);
    t1->SetBranchAddress("reverseMapVector2", &pReverseMapVector2, &pBranchReverseMapVector2);

    const auto tEntry = t1->LoadTree(0);
    pBranchMapVector1->GetEntry(tEntry);
    pBranchMapVector2->GetEntry(tEntry);
    pBranchReverseMapVector1->GetEntry(tEntry);
    pBranchReverseMapVector2->GetEntry(tEntry);

    for (int j = 0; j < (int)pMapVector1->size(); ++j)
    {
        mapVector1.push_back(pMapVector1->at(j));
        mapVector2.push_back(pMapVector2->at(j));
        reverseMapVector1.push_back(pReverseMapVector1->at(j));
        reverseMapVector2.push_back(pReverseMapVector2->at(j));
    }

    std::map<int, double> map;
    std::map<double, int> reverseMap;
    double binWidth;
    double initialEnergy;
    double maxRange;

    t1->SetBranchAddress("binWidth", &binWidth);
    t1->SetBranchAddress("initialEnergy", &initialEnergy);
    t1->SetBranchAddress("maxRange", &maxRange);
    t1->GetEntry(0);

    for (int i = 0; i < mapVector1.size(); i++)
        map.insert(std::pair<int,double>(mapVector1.at(i), mapVector2.at(i)));

    for (int i = 0; i < mapVector1.size(); i++)
        reverseMap.insert(std::pair<double, int>(reverseMapVector1.at(i), reverseMapVector2.at(i)));


    lookupTable.SetMap(map);
    lookupTable.SetReverseMap(reverseMap);
    lookupTable.SetBinWidth(binWidth);
    lookupTable.SetInitialEnergy(initialEnergy);
    lookupTable.SetMaxRange(maxRange);

    f->Close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SetEndpoints(DirectionFitObject &fitResult, const Cluster *const pCluster)
{
    CartesianVector lowZVector(0.f, 0.f, 0.f), highZVector(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(pCluster, lowZVector, highZVector);

    fitResult.SetBeginpoint(lowZVector);
    fitResult.SetEndpoint(highZVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SetEndpoints(DirectionFitObject &fitResult, const LArTrackStateVector &trackStateVector)
{
    TrackState firstTrackState(*(trackStateVector.begin())), lastTrackState(trackStateVector.back());
    const pandora::CartesianVector initialPosition(firstTrackState.GetPosition());
    const pandora::CartesianVector endPosition(lastTrackState.GetPosition());

    const pandora::CartesianVector lowZVector(initialPosition.GetZ() < endPosition.GetZ() ? initialPosition : endPosition);
    const pandora::CartesianVector highZVector(initialPosition.GetZ() > endPosition.GetZ() ? initialPosition : endPosition);

    fitResult.SetBeginpoint(lowZVector);
    fitResult.SetEndpoint(highZVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SetMCTruth(DirectionFitObject &fitResult, const Cluster *const pCluster)
{
    CartesianVector mcEndpoint(MCParticleHelper::GetMainMCParticle(pCluster)->GetEndpoint());
    CartesianVector mcBeginpoint(MCParticleHelper::GetMainMCParticle(pCluster)->GetVertex());
    CartesianVector mcDirection((mcEndpoint - mcBeginpoint).GetUnitVector());

    CartesianVector recoBeginpoint(fitResult.GetBeginpoint());
    CartesianVector recoEndpoint(fitResult.GetEndpoint());

    if ((mcBeginpoint - recoBeginpoint).GetMagnitude() < (mcBeginpoint - recoEndpoint).GetMagnitude())
        fitResult.SetMCDirection(1);
    else
        fitResult.SetMCDirection(0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FillHitChargeVector(const Cluster *const pCluster, HitChargeVector &hitChargeVector)
{
    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);

    const TwoDSlidingFitResult &slidingFit(this->GetCachedSlidingFit(pCluster));

    for (CaloHitList::const_iterator hitIter = caloHitList.begin(), hitIterEnd = caloHitList.end(); hitIter != hitIterEnd; ++hitIter)
    {
        const CaloHit *const pCaloHit(*hitIter);
        const CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
        float hitWidth(pCaloHit->GetCellSize1());
        float caloHitEnergy(pCaloHit->GetInputEnergy());

        if (!((hitWidth > 0.3 && hitWidth < 1.8) && (caloHitEnergy > 50.0 && caloHitEnergy < 600.0)))
            continue;

        caloHitEnergy *= 187.6; //ADC to electron
        caloHitEnergy *= 23.6/1000000; //ionisation energy per electron in MeV
        caloHitEnergy /= 0.62;

        if (caloHitEnergy/hitWidth < 1.05)
            continue;

        float rL(0.f), rT(0.f);
        slidingFit.GetLocalPosition(caloHitPosition, rL, rT);
        if (rL == 0.)
            continue;

        float calibratedUncertainty(std::sqrt((0.00303236 * (caloHitEnergy) * (caloHitEnergy)) + (0.0038681 * (caloHitEnergy)))); //70%

        HitCharge hitCharge(pCaloHit, rL, hitWidth, caloHitEnergy, calibratedUncertainty);
        hitChargeVector.push_back(hitCharge);
    }

    std::sort(hitChargeVector.begin(), hitChargeVector.end(), SortHitChargeVectorByRL);

}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::TrackInnerFilter(HitChargeVector &hitChargeVector, HitChargeVector &filteredHitChargeVector)
{
    //Fill endpoint protected area into filtered vector and put all other hits in a separate vector
    float endpointProtectionRange(m_endpointProtectionFraction);

    //Smarter filtering
    /*
    HitChargeVector innerHitChargeVector;

    if (endpointProtectionRange != 0.f)
    {
        HitChargeVector leftEndpoint(hitChargeVector.begin(), std::next(hitChargeVector.begin(), endpointProtectionRange * hitChargeVector.size()));
        HitChargeVector rightEndpoint(std::next(hitChargeVector.begin(), (1.0 - endpointProtectionRange) * hitChargeVector.size()), hitChargeVector.end());

        float leftCharge(0.f), rightCharge(0.f);

        for (const auto &hitCharge : leftEndpoint)
            leftCharge += hitCharge.GetChargeOverWidth();

        for (const auto &hitCharge : rightEndpoint)
            rightCharge += hitCharge.GetChargeOverWidth();


        //protect the endpoint with the largest average charge
        if (rightCharge >= leftCharge)
        {
            filteredHitChargeVector.insert(filteredHitChargeVector.begin(), hitChargeVector.begin() + (1.0 - endpointProtectionRange) * hitChargeVector.size(), hitChargeVector.end());
            innerHitChargeVector.insert(innerHitChargeVector.begin(), hitChargeVector.begin(), std::next(hitChargeVector.begin(), (1.0 - endpointProtectionRange) * hitChargeVector.size()));
        }
        else
        {
            filteredHitChargeVector.insert(filteredHitChargeVector.begin(), hitChargeVector.begin(), std::next(hitChargeVector.begin(), endpointProtectionRange * hitChargeVector.size()));
            innerHitChargeVector.insert(innerHitChargeVector.begin(), std::next(hitChargeVector.begin(), endpointProtectionRange * hitChargeVector.size()), hitChargeVector.end());
        }
    }
    else
    {
        innerHitChargeVector = hitChargeVector;
    }
    */

    filteredHitChargeVector.insert(filteredHitChargeVector.begin(), hitChargeVector.begin(),  hitChargeVector.begin() + endpointProtectionRange * hitChargeVector.size());
    filteredHitChargeVector.insert(filteredHitChargeVector.begin(), hitChargeVector.begin() + (1.0 - endpointProtectionRange) * hitChargeVector.size(), hitChargeVector.end());
    HitChargeVector innerHitChargeVector(hitChargeVector.begin() + endpointProtectionRange * hitChargeVector.size(), hitChargeVector.begin() + (1.0 - endpointProtectionRange) * hitChargeVector.size());

    if (innerHitChargeVector.size() < 10)
    {
        filteredHitChargeVector = hitChargeVector;
        return;
    }

    int nNeighboursToConsider(5);
    this->SetNearestNeighbourValues(innerHitChargeVector, nNeighboursToConsider);

     std::sort(innerHitChargeVector.begin(), innerHitChargeVector.end(), SortByDistanceToNN);
     filteredHitChargeVector.insert(filteredHitChargeVector.begin(), innerHitChargeVector.begin(), innerHitChargeVector.begin() + 0.72 * innerHitChargeVector.size()); //lots of testing has been done to optimise percentage
     std::sort(filteredHitChargeVector.begin(), filteredHitChargeVector.end(), SortHitChargeVectorByRL);

    //-------------

    /*
    if (m_writeTable)
    {
        std::sort(innerHitChargeVector.begin(), innerHitChargeVector.end(), SortByDistanceToNN);

        for (int i = 0; i < innerHitChargeVector.size(); ++i)
        {
            auto hitCharge(innerHitChargeVector.at(i)); 
            const float fractionalPosition(static_cast<float>(i)/innerHitChargeVector.size());
            bool isPure(this->IsPureHit(hitCharge));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterHitTree", "FractionalPosition", fractionalPosition));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterHitTree", "IsPure", isPure ? 1 : 0));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), "FilterHitTree"));
        }

         std::sort(filteredHitChargeVector.begin(), filteredHitChargeVector.end(), SortHitChargeVectorByRL);
    }
    */

    int totalNumberHits(hitChargeVector.size()), nTotalHitsFiltered(0), nMuonHits(0), nNonMuonHits(0), nMuonHitsFiltered(0), nNonMuonHitsFiltered(0);

    for (HitCharge hitCharge : hitChargeVector)
    {   
        bool isPure(this->IsPureHit(hitCharge));
        
        if (isPure)
            nMuonHits++;
        else 
            nNonMuonHits++;
        
        int count(0);
        
        for (HitCharge &hitCharge2 : filteredHitChargeVector)
        {   
            if (hitCharge2.GetCharge() == hitCharge.GetCharge())
                count++;
        }
        
        if (count == 0)
        {   
            nTotalHitsFiltered++;
            
            if (isPure)
                nMuonHitsFiltered++;
            else
                nNonMuonHitsFiltered++;
        }
    }

    float fractionImpureHits(((float) nNonMuonHits) / (nNonMuonHits + nMuonHits));

    float fractionFilteredBelongingToMuon(((float)nMuonHitsFiltered)/nTotalHitsFiltered), fractionFilteredBelongingToNonMuon(((float)nNonMuonHitsFiltered)/nTotalHitsFiltered),
    fractionMuonHitsFiltered(((float)nMuonHitsFiltered)/nMuonHits), fractionNonMuonHitsFiltered(((float)nNonMuonHitsFiltered)/nNonMuonHits), fractionClusterHitsFiltered(((float)nTotalHitsFiltered)/totalNumberHits);

    if (m_writeTable)
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "totalNumberHits", totalNumberHits));
        //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "globalTrueTargetParticleEnergy", globalTrueTargetParticleEnergy));
        //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "trackLength", trackLength));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "fractionImpureHits", fractionImpureHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "nTotalHitsFiltered", nTotalHitsFiltered));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "nMuonHits", nMuonHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "nNonMuonHits", nNonMuonHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "nMuonHitsFiltered", nMuonHitsFiltered));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "nNonMuonHitsFiltered", nNonMuonHitsFiltered));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "fractionFilteredBelongingToMuon", fractionFilteredBelongingToMuon));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "fractionFilteredBelongingToNonMuon", fractionFilteredBelongingToNonMuon));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "fractionMuonHitsFiltered", fractionMuonHitsFiltered));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "fractionNonMuonHitsFiltered", fractionNonMuonHitsFiltered));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "fractionClusterHitsFiltered", fractionClusterHitsFiltered));
        //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "fileIdentifier", m_fileIdentifier));
        //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "FilterTree", "eventNumber", m_eventNumber));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), "FilterTree"));
    }    
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackDirectionTool::IsPureHit(HitCharge &hitCharge)
{
    float contributionThreshold(0.15);

    try  
    {    
        int nContributingMCParticles(0);
        const auto pCaloHit(hitCharge.GetCaloHit());
        MCParticleWeightMap mcParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

        for (MCParticleWeightMap::const_iterator mapIter = mcParticleWeightMap.begin(), mapIterEnd = mcParticleWeightMap.end(); mapIter != mapIterEnd; ++mapIter)
        {    
            float contribution(mapIter->second);

            if (contribution > contributionThreshold)
                nContributingMCParticles++;
        }    

        if (MCParticleHelper::GetMainMCParticle(pCaloHit)->GetParticleId() == 13 && nContributingMCParticles == 1)
            return true;
        else 
            return false;
    }    
    catch (...)
    {    
        std::cout << "Purity of hit uncertain." << std::endl;
        return true;
    }    
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SetNearestNeighbourValues(HitChargeVector &innerHitChargeVector, int &nNeighboursToConsider)
{
    float trackLength(0.f);
    this->GetTrackLength(innerHitChargeVector, trackLength);

    std::sort(innerHitChargeVector.begin(), innerHitChargeVector.end(), SortHitChargeVectorByChargeOverWidth);
    float ChargeOverWidthRange((*(std::prev(innerHitChargeVector.end(), 1))).GetChargeOverWidth() - (*innerHitChargeVector.begin()).GetChargeOverWidth());

    for (HitCharge &hitCharge1 : innerHitChargeVector)
    {
        std::vector<float> distancesToNN;

        for (HitCharge &hitCharge2 : innerHitChargeVector)
        {
            if (&hitCharge1 == &hitCharge2)
                continue;

            float ChargeOverWidthDistance((trackLength/ChargeOverWidthRange) * (std::abs(hitCharge1.GetChargeOverWidth() - hitCharge2.GetChargeOverWidth())));
            float Ldistance(std::abs(hitCharge1.GetLongitudinalPosition() - hitCharge2.GetLongitudinalPosition()));
            float distanceToNN(std::sqrt(ChargeOverWidthDistance*ChargeOverWidthDistance + Ldistance*Ldistance));

            distancesToNN.push_back(distanceToNN);
        }

        std::sort(distancesToNN.begin(), distancesToNN.end());
        float nearestNeighboursDistanceSum(std::accumulate(distancesToNN.begin(), distancesToNN.begin() + nNeighboursToConsider, 0.f));
        hitCharge1.SetDistanceToNN(nearestNeighboursDistanceSum);
    }

    std::sort(innerHitChargeVector.begin(), innerHitChargeVector.end(), SortHitChargeVectorByRL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FragmentRemoval(HitChargeVector &hitChargeVector, HitChargeVector &filteredHitChargeVector, float &splitPosition)
{
    float trackLength(0.f);
    this->GetTrackLength(hitChargeVector, trackLength);

    std::vector<JumpObject> jumpsVector;
    this->FindLargestJumps(hitChargeVector, jumpsVector);

    std::vector<JumpObject> peakJumps;
    this->FindPeakJumps(hitChargeVector, jumpsVector);

    this->AttemptFragmentRemoval(hitChargeVector, jumpsVector, filteredHitChargeVector, splitPosition);
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SimpleTrackEndFilter(HitChargeVector &hitChargeVector)
{
    float lowerBound(0.9), upperBound(10.f);

    while ((*(hitChargeVector.begin())).GetChargeOverWidth()/(*(std::next(hitChargeVector.begin(), 1))).GetChargeOverWidth() <= lowerBound || (*(hitChargeVector.begin())).GetChargeOverWidth()/(*(std::next(hitChargeVector.begin(), 1))).GetChargeOverWidth() >= upperBound)
        hitChargeVector.erase(hitChargeVector.begin());

    while ((*(std::prev(hitChargeVector.end(), 1))).GetChargeOverWidth()/(*(std::prev(hitChargeVector.end(), 2))).GetChargeOverWidth() <= lowerBound || (*(std::prev(hitChargeVector.end(), 1))).GetChargeOverWidth()/(*(std::prev(hitChargeVector.end(), 2))).GetChargeOverWidth() >= upperBound)
        hitChargeVector.pop_back();

    /*
    for (auto &hitCharge : hitChargeVector)
    {
        std::cout << hitCharge.GetLongitudinalPosition() << std::endl;
        std::cout << hitCharge.GetChargeOverWidth() << std::endl;
        std::cout << "----------------------" << std::endl;
    }
    */

    /*
    float lowerBound(0.9), upperBound(2.2);

    while ((*(hitChargeVector.begin())).GetChargeOverWidth()/(*(std::next(hitChargeVector.begin(), 1))).GetChargeOverWidth() <= lowerBound || (*(hitChargeVector.begin())).GetChargeOverWidth()/(*(std::next(hitChargeVector.begin(), 1))).GetChargeOverWidth() >= upperBound)
        hitChargeVector.erase(hitChargeVector.begin());

    while ((*(std::prev(hitChargeVector.end(), 1))).GetChargeOverWidth()/(*(std::prev(hitChargeVector.end(), 2))).GetChargeOverWidth() <= lowerBound || (*(std::prev(hitChargeVector.end(), 1))).GetChargeOverWidth()/(*(std::prev(hitChargeVector.end(), 2))).GetChargeOverWidth() >= upperBound)
        hitChargeVector.pop_back();

    //This piece of logic removes hits that have uncharacteristically high or low Q/w values (in tails of Q/w distribution)
    hitChargeVector.erase(
    std::remove_if(hitChargeVector.begin(), hitChargeVector.end(),
        [](HitCharge & hitCharge) { return hitCharge.m_intails; }),
    hitChargeVector.end());
    */

    /*
    //Get track length and Q over W span for last step
    float trackLength(0.f), minQoverW(1e6), maxQoverW(0.f);
    this->GetTrackLength(hitChargeVector, trackLength);

    if (trackLength == 0)
        return;

    for (HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetChargeOverWidth() < minQoverW)
            minQoverW = hitCharge.GetChargeOverWidth();

        if (hitCharge.GetChargeOverWidth() > maxQoverW)
            maxQoverW = hitCharge.GetChargeOverWidth();
    }

    //If there is only one hit in a wide charge range, remove it
    for (HitChargeVector::const_iterator iter = hitChargeVector.begin(); iter != hitChargeVector.end(); )
    {
        bool nearbyCharge(false);

        for (HitCharge &hitCharge : hitChargeVector)
        {
            if (std::abs(hitCharge.GetLongitudinalPosition() - (*iter).GetLongitudinalPosition()) <= 0.025 * trackLength)
                continue;

            if (std::abs(hitCharge.GetChargeOverWidth() - (*iter).GetChargeOverWidth()) <= 0.1 * (maxQoverW - minQoverW))
            {
                nearbyCharge = true;
                break;
            }
        }

        if (!nearbyCharge) 
            iter = hitChargeVector.erase(iter);
        else
            ++iter;
    }
    */
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::TrackEndFilter(HitChargeVector &hitChargeVector, DirectionFitObject &directionFitObject, DirectionFitObject &beforeDirectionFitObject)
{
    HitChargeVector filteredHitChargeVector(hitChargeVector);
    //int beforeNumberHits(hitChargeVector.size());

    float trackLength(0.f);
    this->GetTrackLength(hitChargeVector, trackLength);

    float bodyQoverW(0.f);
    this->GetAverageQoverWTrackBody(hitChargeVector, bodyQoverW);

    int nHitsToSkip(3), counterFromBeginning(-1), counterToEnd(hitChargeVector.size() - 1);
    float trackEndRange(0.05);

    float leftMean(0.f), rightMean(0.f);
    int leftHits(0), rightHits(0);

    for (auto &hitCharge : hitChargeVector)
    {
            if (hitCharge.GetLongitudinalPosition()/trackLength >= (1.0 - trackEndRange))
            {
                rightMean += hitCharge.GetChargeOverWidth();
                ++rightHits;
            }
            else if (hitCharge.GetLongitudinalPosition()/trackLength <= trackEndRange)
            {
                leftMean += hitCharge.GetChargeOverWidth();
                ++leftHits;
            }
    }

    rightMean /= rightHits;
    leftMean /= leftHits;
    

    for (HitChargeVector::const_iterator iter = filteredHitChargeVector.begin(); iter != filteredHitChargeVector.end(); )
    {
        //This counter exists so that the hit charge N hits over never points before begin() or after end(), hence the use of std::min below
        ++counterFromBeginning;
        --counterToEnd;

        HitCharge hitCharge(*iter), nextHitCharge(*std::next(iter, std::min(1, counterToEnd))), plusNHitCharge(*std::next(iter, std::min(nHitsToSkip, counterToEnd))), previousHitCharge(*std::prev(iter, std::min(1, counterFromBeginning))), minusNHitCharge(*std::prev(iter, std::min(nHitsToSkip, counterFromBeginning)));

        if (hitCharge.GetLongitudinalPosition()/trackLength <= trackEndRange || hitCharge.GetLongitudinalPosition()/trackLength >= (1.0 - trackEndRange))
        {
            float nearestRatio(std::max((hitCharge.GetChargeOverWidth()/previousHitCharge.GetChargeOverWidth()), (hitCharge.GetChargeOverWidth()/nextHitCharge.GetChargeOverWidth())));
            float plusMinusNRatio(std::max((hitCharge.GetChargeOverWidth()/minusNHitCharge.GetChargeOverWidth()), (hitCharge.GetChargeOverWidth()/plusNHitCharge.GetChargeOverWidth())));

            float relevantMean(hitCharge.GetLongitudinalPosition()/trackLength <= trackEndRange ? leftMean : rightMean);
            float distanceFromMeanCharge(std::abs(hitCharge.GetChargeOverWidth() - relevantMean));
            float distanceFromBodyQoverW(std::abs(hitCharge.GetChargeOverWidth() - bodyQoverW));

            int isBraggPeakHit(((hitCharge.GetLongitudinalPosition()/trackLength <= trackEndRange && globalTrueDirection == 0) || (hitCharge.GetLongitudinalPosition()/trackLength >= (1.0 - trackEndRange) && globalTrueDirection == 1)) ? 1 : 0);

            if (m_writeTable)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BraggPeakFilterTree", "IsBraggPeakHit", isBraggPeakHit));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BraggPeakFilterTree", "NearestNeighbourRatio", nearestRatio));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BraggPeakFilterTree", "ThirdNearestNeighbourRatio", plusMinusNRatio));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BraggPeakFilterTree", "DistanceFromAverageCharge", distanceFromMeanCharge));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BraggPeakFilterTree", "DistanceFromTrackBody", distanceFromBodyQoverW));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), "BraggPeakFilterTree"));
            }

            //if (distanceFromBodyQoverW >= 4.8 || std::abs(1.0 - nearestRatio) >= 0.4 || std::abs(1.0 - plusMinusNRatio) >= 0.7)
            //§if (distanceFromBodyQoverW >= 1.5)
            if (plusMinusNRatio >= 1.725)
                iter = filteredHitChargeVector.erase(iter);
            else
                ++iter;
        }
        else
        {
            ++iter;
        }
    }

    std::sort(filteredHitChargeVector.begin(), filteredHitChargeVector.end(), SortHitChargeVectorByRL);

    //*******************************************************************************************
    //Write split information

    DirectionFitObject afterDirectionFitObject;
    this->FitHitChargeVector(filteredHitChargeVector, afterDirectionFitObject);

    float tefBeforeChiSquaredPerHit(beforeDirectionFitObject.GetMinChiSquaredPerHit()), tefAfterChiSquaredPerHit(afterDirectionFitObject.GetMinChiSquaredPerHit());
    float chiSquaredPerHitChange(beforeDirectionFitObject.GetMinChiSquaredPerHit() - afterDirectionFitObject.GetMinChiSquaredPerHit());
    bool shouldApply(false);

    float fitParameter1(0.0175365), fitParameter2(1.98053), fitParameter3(-0.359743);

    if (chiSquaredPerHitChange >= fitParameter1 + exp(fitParameter2 + fitParameter3 * floor(static_cast<float>(filteredHitChargeVector.size())/50) ))
        shouldApply = true;

    //float N(beforeNumberHits);

    /*
    std::cout << "TEF chi squared per hit change: " << chiSquaredPerHitChange << std::endl;
    std::cout << "TEF before Forwards chi squared per hit: " << beforeDirectionFitObject.GetForwardsChiSquaredPerHit() << std::endl;
    std::cout << "TEF before Backwards chi squared per hit: " << beforeDirectionFitObject.GetBackwardsChiSquaredPerHit() << std::endl;
    */

    /*
    if (beforeNumberHits < 400 && chiSquaredPerHitChange < (6.0 - ((N/400) * 7.0)))
        shouldApply = false;
    if (beforeNumberHits >= 400 && chiSquaredPerHitChange < 1.0)
        shouldApply = false;
    */

    SplitObject tefObject;
    tefObject.SetSplitPosition(0.f);
    tefObject.SetBeforeMinChiSquaredPerHit(tefBeforeChiSquaredPerHit);
    tefObject.SetBeforeForwardsChiSquaredPerHit(beforeDirectionFitObject.GetForwardsChiSquaredPerHit());
    tefObject.SetBeforeBackwardsChiSquaredPerHit(beforeDirectionFitObject.GetBackwardsChiSquaredPerHit());

    tefObject.SetAfterMinChiSquaredPerHit(tefAfterChiSquaredPerHit);
    tefObject.SetMinChiSquaredPerHitChange(tefBeforeChiSquaredPerHit - tefAfterChiSquaredPerHit);
    tefObject.SetBeforeNHits(hitChargeVector.size());
    tefObject.SetAfterNHits(filteredHitChargeVector.size());
    tefObject.SetSplitApplied(shouldApply);
    tefObject.SetBeforeDeltaChiSquaredPerHit(beforeDirectionFitObject.GetDeltaChiSquaredPerHit());

    if (m_useMCInformation)
    {
        //bool splitCorrect((beforeDirectionFitObject.GetDirectionEstimate() != globalTrueDirection && afterDirectionFitObject.GetDirectionEstimate() == globalTrueDirection) ? true : false);

        int splitCorrect(-2);
        
        if (beforeDirectionFitObject.GetDirectionEstimate() == afterDirectionFitObject.GetDirectionEstimate()) 
            splitCorrect = -1;

        if (beforeDirectionFitObject.GetDirectionEstimate() == globalTrueDirection && afterDirectionFitObject.GetDirectionEstimate() != globalTrueDirection) 
            splitCorrect = 0;

        if (beforeDirectionFitObject.GetDirectionEstimate() != globalTrueDirection && afterDirectionFitObject.GetDirectionEstimate() == globalTrueDirection) 
            splitCorrect = 1;

        tefObject.SetSplitCorrect(splitCorrect);
        tefObject.SetBeforeForwardsChiSquaredPerHit(beforeDirectionFitObject.GetForwardsChiSquaredPerHit());
        tefObject.SetBeforeBackwardsChiSquaredPerHit(beforeDirectionFitObject.GetBackwardsChiSquaredPerHit());
    }

    directionFitObject.SetTEFObject(tefObject);

    if (shouldApply)
        hitChargeVector = filteredHitChargeVector;
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::Regularise(HitChargeVector &hitChargeVector, HitChargeVector &filteredHitChargeVector)
{
    std::sort(hitChargeVector.begin(), hitChargeVector.end(), SortHitChargeVectorByChargeOverWidth);

    int numberHitsToRemove(std::min(3, static_cast<int>((hitChargeVector.size()/2 - 1))));
    float chargeLowerBound(hitChargeVector.at((numberHitsToRemove - 1)).GetChargeOverWidth()), chargeUpperBound(hitChargeVector.at(hitChargeVector.size() - (numberHitsToRemove + 1)).GetChargeOverWidth());

    std::sort(hitChargeVector.begin(), hitChargeVector.end(), SortHitChargeVectorByRL);

    for (const auto &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetChargeOverWidth() == hitChargeVector.front().GetChargeOverWidth() || hitCharge.GetChargeOverWidth() == hitChargeVector.back().GetChargeOverWidth())
            continue;

        if (hitCharge.GetChargeOverWidth() > chargeLowerBound && hitCharge.GetChargeOverWidth() < chargeUpperBound)
            filteredHitChargeVector.push_back(hitCharge);
    }

    std::sort(filteredHitChargeVector.begin(), filteredHitChargeVector.end(), SortHitChargeVectorByRL);
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::AttemptFragmentRemoval(HitChargeVector &hitChargeVector, std::vector<JumpObject> &jumpsVector, HitChargeVector &filteredHitChargeVector, float &finalSplitPosition)
{
    DirectionFitObject beforeDirectionFitObject;
    this->FitHitChargeVector(hitChargeVector, beforeDirectionFitObject);

    float bestSplitPosition(0.f);

    HitChargeVector bestHitChargeVector;
    DirectionFitObject bestDirectionFitObject(beforeDirectionFitObject);

    for (JumpObject &jumpObject : jumpsVector)
    {
        float splitPosition(jumpObject.GetLongitudinalPosition());

        HitChargeVector smallHitCollection, largeHitCollection;
        this->SplitHitCollectionBySize(hitChargeVector, splitPosition, smallHitCollection, largeHitCollection);

        DirectionFitObject afterDirectionFitObject;
        this->FitHitChargeVector(largeHitCollection, afterDirectionFitObject);

        if (afterDirectionFitObject.GetMinChiSquaredPerHit() < bestDirectionFitObject.GetMinChiSquaredPerHit())
        {
            bestSplitPosition = splitPosition;
            bestHitChargeVector = largeHitCollection;
            bestDirectionFitObject = afterDirectionFitObject;
        }
    }

    finalSplitPosition = bestSplitPosition;
    filteredHitChargeVector = bestHitChargeVector;
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindLargestJumps(HitChargeVector &hitChargeVector, std::vector<JumpObject> &normalJumps)
{
    //HitChargeVector binnedHitChargeVector;
    //BinHitChargeVector(hitChargeVector, binnedHitChargeVector, 0.5);

    std::vector<HitChargeVector> bothVectors;
    bothVectors.push_back(hitChargeVector);
    //bothVectors.push_back(binnedHitChargeVector);

    for (HitChargeVector &vector : bothVectors)
    {
        int searchRange(0.05 * vector.size());

        for (int jumpRange = 1; jumpRange <= 5; jumpRange++)
        {
            for (int i = 0; i < searchRange; i++)
            {
                float binJump = (std::abs(vector.at(i).GetChargeOverWidth() - vector.at(i + jumpRange).GetChargeOverWidth()));
                float jumpPosition(vector.at(i + jumpRange).GetLongitudinalPosition());
                JumpObject jumpObject(jumpPosition, binJump);
                normalJumps.push_back(jumpObject);
            }

            for (int j = vector.size() - searchRange; j < vector.size() - jumpRange; j++)
            {
                float binJump = (std::abs(vector.at(j).GetChargeOverWidth() - vector.at(j + jumpRange).GetChargeOverWidth()));
                float jumpPosition(vector.at(j).GetLongitudinalPosition());
                JumpObject jumpObject(jumpPosition, binJump);

                normalJumps.push_back(jumpObject);
            }
        }
    }

    std::sort(normalJumps.begin(), normalJumps.end(), SortJumpVector);
    if (normalJumps.size() > 3)
        normalJumps.erase(normalJumps.begin() + 3, normalJumps.end());
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindPeakJumps(HitChargeVector &hitChargeVector, std::vector<JumpObject> &peakJumps)
{
    float jumpPosition(0.f), jumpValue(0.f), currentLargestQoverW(0.f);

    for (HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetChargeOverWidth() > currentLargestQoverW)
        {
            currentLargestQoverW = hitCharge.GetChargeOverWidth();
            jumpPosition = hitCharge.GetLongitudinalPosition();
        }
    }

    float jumpPosition1(jumpPosition - 0.5), jumpPosition2(jumpPosition + 0.5);
    JumpObject jumpObject(jumpPosition, jumpValue);
    JumpObject jumpObject1(jumpPosition1, jumpValue);
    JumpObject jumpObject2(jumpPosition2, jumpValue);
    peakJumps.push_back(jumpObject);
    peakJumps.push_back(jumpObject1);
    peakJumps.push_back(jumpObject2);
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindTrackEndJumps(HitChargeVector &hitChargeVector, std::vector<JumpObject> &trackEndJumps)
{
    float trackLength(0.f);
    this->GetTrackLength(hitChargeVector, trackLength);

    for (float edge = 0.01; edge <= 0.25; edge += 0.01)
    {
        float jumpPosition1(edge * trackLength), jumpPosition2((1.0 - edge) * trackLength), jumpValue(0.f);
        JumpObject jumpObject1(jumpPosition1, jumpValue);
        JumpObject jumpObject2(jumpPosition2, jumpValue);

        trackEndJumps.push_back(jumpObject1);
        trackEndJumps.push_back(jumpObject2);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::ParticleSplitting(const Cluster* pCluster, HitChargeVector &hitChargeVector, DirectionFitObject &backwardsDirectionFitObject, DirectionFitObject &forwardsDirectionFitObject, bool &splitApplied, SplitObject &splitObject)
{
    DirectionFitObject beforeDirectionFitObject;
    this->FitHitChargeVector(hitChargeVector, beforeDirectionFitObject);

    DirectionFitObject outputBackwardsDirectionFitObject(beforeDirectionFitObject), outputForwardsDirectionFitObject(beforeDirectionFitObject);
    backwardsDirectionFitObject = beforeDirectionFitObject;
    forwardsDirectionFitObject = beforeDirectionFitObject;

    float afterSplitChiSquared(beforeDirectionFitObject.GetMinChiSquaredPerHit()), bestSplitPosition(0.f);

    //std::vector<float> calorimetricSplitPositions;
    //this->CreateCalorimetricSplitHitVector(hitChargeVector, calorimetricSplitPositions);

    std::vector<JumpObject> jumpObjects;
    this->FindJumpSplit(hitChargeVector, jumpObjects);

    float bestJumpSize(0.f);

    for (auto &jumpObject : jumpObjects)
    {
        float splitPosition(jumpObject.GetLongitudinalPosition());

        HitChargeVector backwardsTestHitCollection, forwardsTestHitCollection;
        this->SplitHitCollectionByLeftRight(hitChargeVector, splitPosition, backwardsTestHitCollection, forwardsTestHitCollection);

        HitChargeVector backwardsTestFilteredHitCollection, forwardsFilteredTestHitCollection;
        this->TrackInnerFilter(backwardsTestHitCollection, backwardsTestFilteredHitCollection);
        this->TrackInnerFilter(forwardsTestHitCollection, forwardsFilteredTestHitCollection);

        DirectionFitObject backwardsTestDirectionFitObject, forwardsTestDirectionFitObject;
        this->FitHitChargeVector(backwardsTestFilteredHitCollection, forwardsFilteredTestHitCollection, backwardsTestDirectionFitObject, forwardsTestDirectionFitObject);

        float splitMinChiSquared((backwardsTestDirectionFitObject.GetNHits() > 0 ? backwardsTestDirectionFitObject.GetBackwardsChiSquared()/backwardsTestDirectionFitObject.GetNHits() : 0.f) + (forwardsTestDirectionFitObject.GetNHits() > 0 ? forwardsTestDirectionFitObject.GetForwardsChiSquared()/forwardsTestDirectionFitObject.GetNHits() : 0.f));

        //float kinkSize(0.f);
        //this->FindKinkSize(pTargetClusterW, splitPosition, kinkSize);

        if (splitMinChiSquared < afterSplitChiSquared)
        {
            bestJumpSize = jumpObject.GetJumpValue();
            afterSplitChiSquared = splitMinChiSquared;
            bestSplitPosition = splitPosition;
            outputBackwardsDirectionFitObject = backwardsTestDirectionFitObject;
            outputForwardsDirectionFitObject = forwardsTestDirectionFitObject;
        }
    }

    if (outputBackwardsDirectionFitObject.GetNHits() == 0 || outputForwardsDirectionFitObject.GetNHits() == 0)
        return;

    float chiSquaredPerHitChange(beforeDirectionFitObject.GetMinChiSquaredPerHit() - (outputBackwardsDirectionFitObject.GetBackwardsChiSquared()/outputBackwardsDirectionFitObject.GetNHits() + outputForwardsDirectionFitObject.GetForwardsChiSquared()/outputForwardsDirectionFitObject.GetNHits()));
    //int beforeNumberHits((int)beforeDirectionFitObject.GetHitChargeVector().size());
    //float N(beforeNumberHits);
    bool shouldApply(false);

    /*
    if (beforeNumberHits < 400 && chiSquaredPerHitChange < (5.65777 - (0.586666/50) * N))
        shouldApply = false;
    if (beforeNumberHits >= 400 && chiSquaredPerHitChange < 1.0)
        shouldApply = false;
    */

    float fitParameter1(0.384288), fitParameter2(1.10873), fitParameter3(-0.983149);

    if (bestJumpSize >= 0.6 && (chiSquaredPerHitChange >= fitParameter1 + exp(fitParameter2 + fitParameter3 * floor(static_cast<float>(hitChargeVector.size())/50) )))
        shouldApply = true;

    if (m_useMCInformation)
    {
        bool splitCorrect(this->IsClusterTwoParticles(pCluster, outputBackwardsDirectionFitObject.GetBackwardsFitCharges(), outputForwardsDirectionFitObject.GetForwardsFitCharges()));
        splitObject.SetSplitCorrect(splitCorrect);
        splitObject.SetBeforeForwardsChiSquaredPerHit(beforeDirectionFitObject.GetForwardsChiSquaredPerHit());
        splitObject.SetBeforeBackwardsChiSquaredPerHit(beforeDirectionFitObject.GetBackwardsChiSquaredPerHit());
    }

    /*
    std::cout << ">>> Splitting before chi squared: " << beforeDirectionFitObject.GetMinChiSquaredPerHit() << std::endl;
    std::cout << ">>> Splitting chi squared per hit change: " << chiSquaredPerHitChange << std::endl;
    std::cout << ">>> Splitting before Forwards chi squared per hit: " << beforeDirectionFitObject.GetForwardsChiSquaredPerHit() << std::endl;
    std::cout << ">>> Splitting before Backwards chi squared per hit: " << beforeDirectionFitObject.GetBackwardsChiSquaredPerHit() << std::endl;
    */

    splitObject.SetJumpSize(bestJumpSize);
    splitObject.SetSplitPosition(bestSplitPosition);
    splitObject.SetBeforeMinChiSquaredPerHit(beforeDirectionFitObject.GetMinChiSquaredPerHit());
    splitObject.SetBeforeForwardsChiSquaredPerHit(beforeDirectionFitObject.GetForwardsChiSquaredPerHit());
    splitObject.SetBeforeBackwardsChiSquaredPerHit(beforeDirectionFitObject.GetBackwardsChiSquaredPerHit());

    splitObject.SetAfterMinChiSquaredPerHit(beforeDirectionFitObject.GetMinChiSquaredPerHit() - chiSquaredPerHitChange);
    splitObject.SetMinChiSquaredPerHitChange(chiSquaredPerHitChange);
    splitObject.SetBeforeNHits(hitChargeVector.size());
    splitObject.SetSplitApplied(splitApplied);
    splitObject.SetBeforeDeltaChiSquaredPerHit(beforeDirectionFitObject.GetDeltaChiSquaredPerHit());

    if (shouldApply)
    {
        splitApplied = true;
        splitObject.SetSplitPosition(bestSplitPosition);
        backwardsDirectionFitObject = outputBackwardsDirectionFitObject;
        forwardsDirectionFitObject = outputForwardsDirectionFitObject;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindKinkSize(const Cluster* pCluster, float &splitPosition, float &kinkSize)
{
    try
    {
        const TwoDSlidingFitResult &slidingFit(this->GetCachedSlidingFit(pCluster));
        const LayerFitResultMap &layerFitResultMap(slidingFit.GetLayerFitResultMap());
        const int minLayer(layerFitResultMap.begin()->first), maxLayer(layerFitResultMap.rbegin()->first);

        const int nLayersHalfWindow(slidingFit.GetLayerFitHalfWindow());
        const int nLayersSpanned(1 + maxLayer - minLayer);

        if (nLayersSpanned <= 2 * nLayersHalfWindow)
            return;

        for (LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end(); iter != iterEnd; ++iter)
        {
            const int iLayer(iter->first);

            const float rL(slidingFit.GetL(iLayer));
            const float rL1(slidingFit.GetL(iLayer - nLayersHalfWindow));
            const float rL2(slidingFit.GetL(iLayer + nLayersHalfWindow));

            CartesianVector centralPosition(0.f,0.f,0.f), firstDirection(0.f,0.f,0.f), secondDirection(0.f,0.f,0.f);

            if ((STATUS_CODE_SUCCESS != slidingFit.GetGlobalFitPosition(rL, centralPosition)) ||
                (STATUS_CODE_SUCCESS != slidingFit.GetGlobalFitDirection(rL1, firstDirection)) ||
                (STATUS_CODE_SUCCESS != slidingFit.GetGlobalFitDirection(rL2, secondDirection)))
            {
                continue;
            }

            const float cosTheta(firstDirection.GetDotProduct(secondDirection));
            if (std::abs(splitPosition - rL) <= 3.0 && cosTheta > kinkSize && cosTheta < 1.0)
                kinkSize = (180.0 / 3.1415926535) * std::acos(cosTheta);
        }
    }
    catch (...)
    {
        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::CreateCalorimetricSplitHitVector(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions)
{
    std::vector<JumpObject> jumpObjects;
    this->FindJumpSplit(hitChargeVector, jumpObjects);

    //this->FindPlateauSplit(hitChargeVector, jumpObjects);

    for (auto &jumpObject: jumpObjects)
        splitPositions.push_back(jumpObject.GetLongitudinalPosition());

    /*
    float QoverWRange(0.f);
    this->GetQoverWRange(hitChargeVector, QoverWRange);

    //only try kink split if no jumps found
    auto it = find_if(jumpObjects.begin(), jumpObjects.end(), [&QoverWRange](JumpObject& obj) {return obj.GetJumpValue() > 0.1 * QoverWRange;});

    if (it == jumpObjects.end())
        this->FindKinkSplit(hitChargeVector, splitPositions);
    */

    //this->FindSamplingSplit(hitChargeVector, splitPositions);

    std::sort(splitPositions.begin(), splitPositions.end());
    splitPositions.erase(std::unique(splitPositions.begin(), splitPositions.end()), splitPositions.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SplitHitCollectionBySize(HitChargeVector &hitChargeVector, float &splitPosition, HitChargeVector &smallHitChargeVector, HitChargeVector &largeHitChargeVector)
{
    HitChargeVector leftHitCollection, rightHitCollection;

    for (const  HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetLongitudinalPosition() <= splitPosition)
            leftHitCollection.push_back(hitCharge);
        else
            rightHitCollection.push_back(hitCharge);
    }

    if (leftHitCollection.size() <= rightHitCollection.size())
    {
        smallHitChargeVector = leftHitCollection;
        largeHitChargeVector = rightHitCollection;
    }
    else
    {
        smallHitChargeVector = rightHitCollection;
        largeHitChargeVector = leftHitCollection;
    }
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SplitHitCollectionByLeftRight(HitChargeVector &hitChargeVector, float &splitPosition, HitChargeVector &leftHitChargeVector, HitChargeVector &rightHitChargeVector)
{
    for (const  HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetLongitudinalPosition() <= splitPosition)
            leftHitChargeVector.push_back(hitCharge);
        else
            rightHitChargeVector.push_back(hitCharge);
    }
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::GetTrackLength(HitChargeVector &hitChargeVector, float &trackLength)
{
    trackLength = 0.f;

    if (hitChargeVector.size() == 0)
        return;

    for (const auto hitCharge : hitChargeVector)
    {
        if (hitCharge.GetLongitudinalPosition() > trackLength)
            trackLength = hitCharge.GetLongitudinalPosition();
    }
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::GetAverageQoverWTrackBody(HitChargeVector &hitChargeVector, float &averageChargeTrackBody)
{
    //temp vector because I do not want to mess with the sorting of the original vector
    HitChargeVector tempHitChargeVector;

    for (HitCharge &hitCharge : hitChargeVector)
        tempHitChargeVector.push_back(hitCharge);

    std::sort(tempHitChargeVector.begin(), tempHitChargeVector.end(), SortHitChargeVectorByChargeOverWidth);

    int nEntries(0);

    for (int q = 0; q < tempHitChargeVector.size(); q++)
    {
        if (q <= 0.1 * tempHitChargeVector.size() || q >= 0.6 * tempHitChargeVector.size())
            continue;

        averageChargeTrackBody += tempHitChargeVector.at(q).GetChargeOverWidth();
        nEntries++;
    }

    averageChargeTrackBody /= nEntries;
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::GetQoverWRange(HitChargeVector &hitChargeVector, float &QoverWRange)
{
    //temp vector because I do not want to mess with the sorting of the original vector
    HitChargeVector tempHitChargeVector;

    for (HitCharge &hitCharge : hitChargeVector)
        tempHitChargeVector.push_back(hitCharge);

    std::sort(tempHitChargeVector.begin(), tempHitChargeVector.end(), SortHitChargeVectorByChargeOverWidth);

    float minQoverW(1e6), maxQoverW(0.f);

    for (HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetChargeOverWidth() > maxQoverW)
            maxQoverW = hitCharge.GetChargeOverWidth();

        if (hitCharge.GetChargeOverWidth() < minQoverW)
            minQoverW = hitCharge.GetChargeOverWidth();
    }

    QoverWRange = (maxQoverW - minQoverW);
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindKinkSplit(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions)
{
    HitChargeVector binnedHitChargeVector = hitChargeVector;
    //BinHitChargeVector(hitChargeVector, binnedHitChargeVector, 1.0);

    std::vector<JumpObject> kinkObjects;

    float minCharge(10000.f), maxCharge(0.f);

    for (HitCharge &hitCharge : binnedHitChargeVector)
    {
        if (hitCharge.GetCharge() < minCharge)
            minCharge = hitCharge.GetCharge();

        if (hitCharge.GetCharge() > maxCharge)
            maxCharge = hitCharge.GetCharge();
    }

    float fullChargeRange(maxCharge - minCharge);
    float chargeHalfWidth(0.1 * fullChargeRange);

    for (HitCharge &bin1 : binnedHitChargeVector)
    {
        HitChargeVector leftHitCollection, rightHitCollection;

        for (HitCharge &vector : binnedHitChargeVector)
        {
            if (vector.GetLongitudinalPosition() <= bin1.GetLongitudinalPosition())
                leftHitCollection.push_back(vector);
            else
                rightHitCollection.push_back(vector);
        }

        if (leftHitCollection.size() == 0 || rightHitCollection.size() == 0)
            continue;

        float bestLeftScore(0.f), bestLeftSlope(0.f);

        for (HitCharge &bin2 : leftHitCollection)
        {
            float chargeDifference(bin2.GetCharge() - bin1.GetCharge());
            float positionDifference(bin1.GetLongitudinalPosition() - bin2.GetLongitudinalPosition());
            float slope(chargeDifference/positionDifference);

            int nLeftHits(0);
            for (HitCharge &bin3 : leftHitCollection)
            {
                float lineValue(bin2.GetCharge() - ((bin3.GetLongitudinalPosition() - bin2.GetLongitudinalPosition()) * slope));
                if (std::abs(bin3.GetCharge() - lineValue) < chargeHalfWidth)
                    nLeftHits++;
            }

            float score((float)nLeftHits/leftHitCollection.size());

            if (score > bestLeftScore)
            {
                bestLeftScore = score;
                bestLeftSlope = slope;
            }
        }


        float bestRightScore(0.f), bestRightSlope(0.f);

        for (HitCharge &bin2 : rightHitCollection)
        {
            float chargeDifference(bin2.GetCharge() - bin1.GetCharge());
            float positionDifference(bin1.GetLongitudinalPosition() - bin2.GetLongitudinalPosition());
            float slope(chargeDifference/positionDifference);

            int nRightHits(0);
            for (HitCharge &bin3 : rightHitCollection)
            {
                float lineValue(bin2.GetCharge() - ((bin3.GetLongitudinalPosition() - bin2.GetLongitudinalPosition()) * slope));
                if (std::abs(bin3.GetCharge() - lineValue) < chargeHalfWidth)
                    nRightHits++;
            }

            float score((float)nRightHits/rightHitCollection.size());

            if (score > bestRightScore)
            {
                bestRightScore = score;
                bestRightSlope = slope;
            }
        }

        float kinkPosition(bin1.GetLongitudinalPosition());
        float totalScore(bestLeftScore + bestRightScore);

        CartesianVector leftSlopeVector(1.f, 0.f, bestLeftSlope);
        CartesianVector rightSlopeVector(1.f, 0.f, bestRightSlope);
        float openingAngle(leftSlopeVector.GetOpeningAngle(rightSlopeVector));

        JumpObject kinkObject(kinkPosition, totalScore, openingAngle);
        kinkObjects.push_back(kinkObject);
    }

    std::sort(kinkObjects.begin(), kinkObjects.end(), SortJumpVector);

    int cutOff(3), nAdded(0);
    float latestJumpPosition(0.f), latestJumpValue(0.f), range(3.f);

    for (int i = 0; i < kinkObjects.size(); i++)
    {
        if (nAdded >= cutOff)
            break;

        JumpObject kinkObject(kinkObjects.at(i));

        if (kinkObject.GetJumpValue() < 1.5)
            continue;

        if (nAdded == 0)
        {
            if (kinkObject.GetOpeningAngle() > 0.05)
                splitPositions.push_back(kinkObjects.at(i).GetLongitudinalPosition());

            latestJumpPosition = kinkObjects.at(i).GetLongitudinalPosition();
            latestJumpValue = kinkObjects.at(i).GetJumpValue();
            nAdded++;
        }
        else
        {
            if ((kinkObject.GetLongitudinalPosition() - latestJumpPosition) < range && kinkObject.GetJumpValue() > latestJumpValue)
            {
                latestJumpPosition = kinkObjects.at(i).GetLongitudinalPosition();
                latestJumpValue = kinkObjects.at(i).GetJumpValue();

                if (kinkObject.GetOpeningAngle() > 0.05)
                {
                    splitPositions.pop_back();
                    splitPositions.push_back(kinkObjects.at(i).GetLongitudinalPosition());
                }
            }
            else if ((kinkObject.GetLongitudinalPosition() - latestJumpPosition) > range)
            {
                latestJumpPosition = kinkObjects.at(i).GetLongitudinalPosition();
                latestJumpValue = kinkObjects.at(i).GetJumpValue();
                nAdded++;

                if (kinkObject.GetOpeningAngle() > 0.05)
                    splitPositions.push_back(kinkObjects.at(i).GetLongitudinalPosition());
            }
        }
    }
}
//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindPlateauSplit(HitChargeVector &hitChargeVector, std::vector<JumpObject> &jumpObjects)
{
    float trackLength(0.f);
    this->GetTrackLength(hitChargeVector, trackLength);

    float averageCharge(0.f);
    this->GetAverageQoverWTrackBody(hitChargeVector, averageCharge);

    std::vector<JumpObject> plateauObjects;


    /////////////CHANGE SETTINGS HERE/////////////

    float positionStepSize(0.1);
    float chargeStep(0.10);
    float trackScanRange(0.05 * trackLength);

    //////////////////////////////////////////////

    for (float currentPosition = 0; currentPosition < trackLength; currentPosition += positionStepSize)
    {
        int totalHitsLeft(0), totalHitsRight(0);

        for (HitCharge &hitCharge : hitChargeVector)
        {
            if (std::abs(currentPosition - hitCharge.GetLongitudinalPosition()) > trackScanRange)
                continue;

            if (hitCharge.GetLongitudinalPosition() <= currentPosition)
                totalHitsLeft++;
            else
                totalHitsRight++;
        }

        float bestHitFractionLeft(0), bestHitFractionRight(0);
        float bestChargeLeft(0), bestChargeRight(0);

        for (float currentCharge = chargeStep; currentCharge < 10.0; currentCharge += chargeStep)
        {
            float hitCountLeft(0.f), hitCountRight(0.f);

            for (HitCharge &hitCharge : hitChargeVector)
            {
                if (std::abs(currentPosition - hitCharge.GetLongitudinalPosition()) > trackScanRange)
                    continue;

                if (hitCharge.GetLongitudinalPosition() <= currentPosition && hitCharge.GetCharge() > (currentCharge - chargeStep) && hitCharge.GetCharge() < (currentCharge + chargeStep))
                    hitCountLeft += 1.0;

                if (hitCharge.GetLongitudinalPosition() > currentPosition && hitCharge.GetCharge() > (currentCharge - chargeStep) && hitCharge.GetCharge() < (currentCharge + chargeStep))
                    hitCountRight += 1.0;
            }

            float hitFractionLeft(hitCountLeft/totalHitsLeft), hitFractionRight(hitCountRight/totalHitsRight);

            if (hitFractionLeft > bestHitFractionLeft) // && hitCountLeft >= 5
            {
                bestHitFractionLeft = hitFractionLeft;
                bestChargeLeft = currentCharge;
            }

            if (hitFractionRight > bestHitFractionRight) //&& hitCountRight >= 5
            {
                bestHitFractionRight = hitFractionRight;
                bestChargeRight = currentCharge;
            }
        }

        float chargeRange(std::abs(bestChargeLeft - bestChargeRight));

        float currentScore(bestHitFractionLeft + bestHitFractionRight);
        currentScore *= chargeRange/averageCharge;
        JumpObject plateauObject(currentPosition, currentScore);

        plateauObjects.push_back(plateauObject);
    }

    std::sort(plateauObjects.begin(), plateauObjects.end(), SortJumpVector);
    int cutOff(3), nAdded(0);
    float latestJumpPosition(0.f), latestJumpValue(0.f), range(3.f);

    for (int i = 0; i < plateauObjects.size(); i++)
    {
        if (nAdded >= cutOff)
            break;

        JumpObject plateauObject(plateauObjects.at(i));

        if (plateauObject.GetJumpValue() < 0.5)
            continue;

        if (nAdded == 0)
        {
            jumpObjects.push_back(plateauObjects.at(i));
            latestJumpPosition = plateauObjects.at(i).GetLongitudinalPosition();
            latestJumpValue = plateauObjects.at(i).GetJumpValue();
            nAdded++;
        }
        else
        {
            if ((plateauObject.GetLongitudinalPosition() - latestJumpPosition) < range && plateauObject.GetJumpValue() > latestJumpValue)
            {
                jumpObjects.pop_back();
                jumpObjects.push_back(plateauObjects.at(i));
                latestJumpPosition = plateauObjects.at(i).GetLongitudinalPosition();
                latestJumpValue = plateauObjects.at(i).GetJumpValue();
            }
            else if ((plateauObject.GetLongitudinalPosition() - latestJumpPosition) > range)
            {
                jumpObjects.push_back(plateauObjects.at(i));
                latestJumpPosition = plateauObjects.at(i).GetLongitudinalPosition();
                latestJumpValue = plateauObjects.at(i).GetJumpValue();
                nAdded++;
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindJumpSplit(HitChargeVector &hitChargeVector, std::vector<JumpObject> &jumpObjects)
{
    //HitChargeVector binnedHitChargeVector;
    //BinHitChargeVector(hitChargeVector, binnedHitChargeVector, 1.0);

    std::vector<JumpObject> normalJumps, binnedJumps;

    std::vector<HitChargeVector> bothVectors;
    bothVectors.push_back(hitChargeVector);
    //bothVectors.push_back(binnedHitChargeVector);

    float minimumJumpOffset(3.f);

    for (HitChargeVector &vector : bothVectors)
    {
        for (int jumpRange = 1; jumpRange <= 5; jumpRange++)
        {
            for (int i = 0; i <= vector.size() - (jumpRange + 1); i++)
            {
                //float combinedUncertainty = vector.at(i).GetUncertainty() + vector.at(i + jumpRange).GetUncertainty();
                float jumpSize = (std::abs(vector.at(i).GetChargeOverWidth() - vector.at(i + jumpRange).GetChargeOverWidth()));
                //jumpSize /= combinedUncertainty;

                //if (jumpSize < 1.0)
                //    continue;

                //std::cout << normalJumps.size() << std::endl;

                float jumpPosition(vector.at(i + jumpSize).GetLongitudinalPosition());
                JumpObject jumpObject(jumpPosition, jumpSize);

                if (normalJumps.size() == 0)
                    normalJumps.push_back(jumpObject);
                else
                {
                    JumpObject closestJumpObject;
                    float smallestPositionOffset(1e6);

                    for (auto &jumpObjectIter : normalJumps)
                    {
                        if (std::abs(jumpObjectIter.GetLongitudinalPosition() - jumpPosition) < smallestPositionOffset)
                        {
                            closestJumpObject = jumpObjectIter;
                            smallestPositionOffset = std::abs(jumpObjectIter.GetLongitudinalPosition() - jumpPosition);
                        }
                    }

                    float closestPosition(closestJumpObject.GetLongitudinalPosition()), closestJumpValue(closestJumpObject.GetJumpValue());
                    
                    if (std::abs(jumpPosition - closestPosition) < minimumJumpOffset && jumpSize > closestJumpValue)
                    {
                        normalJumps.erase(std::remove(normalJumps.begin(), normalJumps.end(), closestJumpObject), normalJumps.end());
                        normalJumps.push_back(jumpObject);
                        closestPosition = jumpPosition;
                        closestJumpValue = jumpSize;
                    }
                    else if (std::abs(jumpPosition - closestPosition) >= minimumJumpOffset)
                    {
                        normalJumps.push_back(jumpObject);
                        closestPosition = jumpPosition;
                        closestJumpValue = jumpSize;
                    }
                }

                //if (vector.size() == hitChargeVector.size())
                //    normalJumps.push_back(jumpObject);

                //if (vector.size() == binnedHitChargeVector.size())
                //    binnedJumps.push_back(jumpObject);
            }
        }
    }

    jumpObjects = normalJumps;   

    /*
    std::sort(normalJumps.begin(), normalJumps.end(), SortJumpVector);
    std::sort(binnedJumps.begin(), binnedJumps.end(), SortJumpVector);

    int cutOff(10); 
    //int nAdded(0);
    //float latestJumpPosition(0.f), latestJumpValue(0.f), range(3.f);

    for (int i = 0; i < cutOff; i++)
    {
        JumpObject jumpObject(normalJumps.at(i));
        jumpObjects.push_back(normalJumps.at(i));
    }
     */

    /*
    for (int i = 0; i < normalJumps.size(); i++)
    {
        if (nAdded >= cutOff)
            break;

        std::cout << normalJumps.at(i).GetLongitudinalPosition() << std::endl;

        JumpObject jumpObject(normalJumps.at(i));

        if (jumpObject.GetJumpValue() < 1.0)
            continue;

        if (nAdded == 0)
        {
            jumpObjects.push_back(normalJumps.at(i));
            latestJumpPosition = normalJumps.at(i).GetLongitudinalPosition();
            latestJumpValue = normalJumps.at(i).GetJumpValue();
            nAdded++;
        }
        else
        {
            std::cout << std::abs(jumpObject.GetLongitudinalPosition() - latestJumpPosition) << std::endl;
            if (std::abs(jumpObject.GetLongitudinalPosition() - latestJumpPosition) < range && jumpObject.GetJumpValue() > latestJumpValue)
            {
                jumpObjects.pop_back();
                jumpObjects.push_back(normalJumps.at(i));
                latestJumpPosition = normalJumps.at(i).GetLongitudinalPosition();
                latestJumpValue = normalJumps.at(i).GetJumpValue();
            }
            else if ((jumpObject.GetLongitudinalPosition() - latestJumpPosition) > range)
            {
                jumpObjects.push_back(normalJumps.at(i));
                latestJumpPosition = normalJumps.at(i).GetLongitudinalPosition();
                latestJumpValue = normalJumps.at(i).GetJumpValue();
                nAdded++;
            }
        }
    }
    */

    /*
    nAdded = 0;

    for (int i = 0; i < binnedJumps.size(); i++)
    {
        if (nAdded >= cutOff)
            break;

        JumpObject jumpObject(binnedJumps.at(i));

        if (jumpObject.GetJumpValue() < 1.0)
            continue;

        if (nAdded == 0)
        {
            jumpObjects.push_back(normalJumps.at(i));
            latestJumpPosition = binnedJumps.at(i).GetLongitudinalPosition();
            latestJumpValue = binnedJumps.at(i).GetJumpValue();
            nAdded++;
        }
        else
        {
            if ((jumpObject.GetLongitudinalPosition() - latestJumpPosition) < range && jumpObject.GetJumpValue() > latestJumpValue)
            {
                jumpObjects.pop_back();
                jumpObjects.push_back(normalJumps.at(i));
                latestJumpPosition = binnedJumps.at(i).GetLongitudinalPosition();
                latestJumpValue = binnedJumps.at(i).GetJumpValue();
            }
            else if ((jumpObject.GetLongitudinalPosition() - latestJumpPosition) > range)
            {
                jumpObjects.push_back(normalJumps.at(i));
                latestJumpPosition = binnedJumps.at(i).GetLongitudinalPosition();
                latestJumpValue = binnedJumps.at(i).GetJumpValue();
                nAdded++;
            }
        }
    }
    */
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindSamplingSplit(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions)
{
    float trackLength(0.f);
    this->GetTrackLength(hitChargeVector, trackLength);

    for (float i = 0.05; i <= 0.95; i += 0.05)
        splitPositions.push_back(i * trackLength);
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FitHitChargeVector(HitChargeVector &hitChargeVector, TrackDirectionTool::DirectionFitObject &fitResult, int numberHitsToConsider)
{
    float particleForwardsChiSquared(0.f), particleBackwardsChiSquared(0.f);
    int numberHits(std::min(2 * numberHitsToConsider, (int)hitChargeVector.size())), particleForwardsFitStatus(-1), particleBackwardsFitStatus(-1);
    HitChargeVector forwardsFitPoints, backwardsFitPoints;
    FitParameters fitParameters;
    this->PerformFits(hitChargeVector, forwardsFitPoints, backwardsFitPoints, numberHitsToConsider, particleForwardsChiSquared, particleBackwardsChiSquared, particleForwardsFitStatus, particleBackwardsFitStatus, fitParameters);

    float mean_dEdx(0.f);
    HitChargeVector thisHitChargeVector = hitChargeVector;
    for (HitCharge &hitCharge : thisHitChargeVector)
        mean_dEdx += hitCharge.GetChargeOverWidth();
    mean_dEdx /= thisHitChargeVector.size();

    std::sort(thisHitChargeVector.begin(), thisHitChargeVector.end(), SortHitChargeVectorByRL);
    std::sort(forwardsFitPoints.begin(), forwardsFitPoints.end(), SortHitChargeVectorByRL);
    std::sort(backwardsFitPoints.begin(), backwardsFitPoints.end(), SortHitChargeVectorByRL);

    DirectionFitObject finalDirectionFitObject(thisHitChargeVector, forwardsFitPoints, backwardsFitPoints, numberHits, mean_dEdx, particleForwardsChiSquared, particleBackwardsChiSquared);
    finalDirectionFitObject.SetFitParameters(fitParameters);

    SplitObject tefObject(fitResult.GetTEFObject());
    finalDirectionFitObject.SetTEFObject(tefObject);

    fitResult = finalDirectionFitObject;
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FitHitChargeVector(HitChargeVector &hitChargeVector1, HitChargeVector &hitChargeVector2, TrackDirectionTool::DirectionFitObject &fitResult1, TrackDirectionTool::DirectionFitObject &fitResult2, int numberHitsToConsider)
{
    this->FitHitChargeVector(hitChargeVector1, fitResult1, numberHitsToConsider);
    this->FitHitChargeVector(hitChargeVector2, fitResult2, numberHitsToConsider);
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::ComputeProbability(DirectionFitObject &fitResult)
{
    float deltaChiSquaredPerHit(fitResult.GetDeltaChiSquaredPerHit()), forwardsChiSquaredPerHit(fitResult.GetForwardsChiSquaredPerHit()), backwardsChiSquaredPerHit(fitResult.GetBackwardsChiSquaredPerHit());

    std::string fileName(m_probabilityFileName.c_str());
    ifstream inputFile(fileName);

    if (deltaChiSquaredPerHit < -15.0 || deltaChiSquaredPerHit > 15.0)
    {
        float probability(deltaChiSquaredPerHit < -15.0 ? 0.75 : 0.25);
        fitResult.SetProbability(probability);
        return;
    }

    if (inputFile)
    {
        TFile *f = new TFile(m_probabilityFileName.c_str());

        TH1F* forwardsDeltaChiSquared = (TH1F*)f->Get("forwardsDeltaChiSquared"); 
        TH1F* backwardsDeltaChiSquared = (TH1F*)f->Get("backwardsDeltaChiSquared"); 

        float forwardsBinEntry = forwardsDeltaChiSquared->GetBinContent(forwardsDeltaChiSquared->GetBin(forwardsChiSquaredPerHit));
        float backwardsBinEntry = backwardsDeltaChiSquared->GetBinContent(backwardsDeltaChiSquared->GetBin(backwardsChiSquaredPerHit));
        float probability(forwardsBinEntry/(forwardsBinEntry + backwardsBinEntry));

        //TO DO: OUT OF RANGE PROBABILITIES

        fitResult.SetProbability(probability);
        f->Close();
    }
    else
    {
        //std::cout << "WARNING: using pre-defined probability values calibrated on CCQEL muons, because probability.root cannot be found. Define a probability m_probabilityFileName by including <FileName>m_probabilityFileName.root</FileName> in the Pandora XML settings file." << std::endl;

        std::map<float, int> deltaChiSquaredPerHitToBinMap = {
        {-15.0, 1}, {-14.625, 2}, {-14.25, 3}, {-13.875, 4}, {-13.5, 5}, {-13.125, 6}, {-12.75, 7}, {-12.375, 8}, {-12.0, 9}, {-11.625, 10},
        {-11.25, 11}, {-10.875, 12}, {-10.5, 13}, {-10.125, 14}, {-9.75, 15}, {-9.375, 16}, {-9.0, 17}, {-8.625, 18}, {-8.25, 19}, {-7.875, 20},
        {-7.5, 21}, {-7.125, 22}, {-6.75, 23}, {-6.375, 24}, {-6.0, 25}, {-5.625, 26}, {-5.25, 27}, {-4.875, 28}, {-4.5, 29}, {-4.125, 30},
        {-3.75, 31}, {-3.375, 33}, {-3.0, 33}, {-2.625, 34}, {-2.25, 35}, {-1.875, 36}, {-1.5, 37}, {-1.125, 38}, {-0.75, 39}, {-0.375, 40},
        {0.0, 41}, {0.375, 42}, {0.75, 43}, {1.125, 44}, {1.5, 45}, {1.875, 46}, {2.25, 47}, {2.625, 48}, {3.0, 49}, {3.375, 50},
        {3.75, 51}, {4.125, 52}, {4.5, 53}, {4.875, 55}, {5.25, 55}, {5.625, 56}, {6.0, 57}, {6.375, 58}, {6.75, 59}, {7.125, 60},
        {7.5, 61}, {7.875, 62}, {8.25, 63}, {8.625, 66}, {9.0, 66}, {9.375, 66}, {9.75, 67}, {10.125, 68}, {10.5, 69}, {10.875, 70},
        {11.25, 71}, {11.625, 72}, {12.0, 73}, {12.375, 77}, {12.75, 77}, {13.125, 77}, {13.5, 77}, {13.875, 78}, {14.25, 79}, {14.625, 80}
        };

        std::map<int, float> binToProbabilityMap = {
        {1, 0.396614}, {2, 0.396614}, {3, 0.567965}, {4, 0.677773}, {5, 0.630863}, {6, 0.567965}, {7, 0.66352}, {8, 0.612035}, {9, 0.66352}, {10, 0.773655},
        {11, 0.743075}, {12, 0.812674}, {13, 0.858101}, {14, 0.829472}, {15, 0.84969}, {16, 0.829472}, {17, 0.895234}, {18, 0.905632}, {19, 0.920437}, {20, 0.931227},
        {21, 0.940389}, {22, 0.945513}, {23, 0.958795}, {24, 0.961112}, {25, 0.965044}, {26, 0.969887}, {27, 0.975667}, {28, 0.981012}, {29, 0.982457}, {30, 0.983119},
        {31, 0.98561}, {32, 0.98807}, {33, 0.989574}, {34, 0.989973}, {35, 0.98897}, {36, 0.944622}, {37, 0.861042}, {38, 0.81822}, {39, 0.78381}, {40, 0.53081},
        {41, 0.31489}, {42, 0.175161}, {44, 0.157666}, {44, 0.081415}, {45, 0.0977991}, {46, 0.0102574}, {47, 0.0107648}, {48, 0.0078804}, {49, 0.00898676}, {50, 0.0112083},
        {51, 0.0108723}, {52, 0.0100676}, {53, 0.0100676}, {54, 0.0113249}, {55, 0.0124953}, {56, 0.0115656}, {57, 0.0124953}, {58, 0.0146878}, {59, 0.0153076}, {60, 0.0208913},
        {61, 0.0217255}, {62, 0.0293406}, {63, 0.0319228}, {64, 0.0271449}, {65, 0.0387419}, {66, 0.0492657}, {67, 0.0676391}, {68, 0.0471319}, {69, 0.041712}, {70, 0.0981396},
        {71, 0.107868}, {72, 0.0831429}, {73, 0.178738}, {74, 0.119737}, {75, 0.107868}, {76, 0.178738}, {77, 0.134541}, {78, 0.521117}, {79, 0.266179}, {80, 0.266179}
        };

        std::map<float, int>::iterator binIter = deltaChiSquaredPerHitToBinMap.lower_bound(deltaChiSquaredPerHit);
        if(binIter != deltaChiSquaredPerHitToBinMap.begin()) {--binIter;}
        int bin((*binIter).second);

        std::map<int, float>::iterator probabilityIter = binToProbabilityMap.lower_bound(bin);
        if(probabilityIter != binToProbabilityMap.begin()) {--probabilityIter;}
        float probability((*probabilityIter).second);

        fitResult.SetProbability(probability);
    }
}

//---------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SetGlobalMinuitPreliminaries(HitChargeVector &hitChargeVector)
{
    this->ClearGlobalVariables();

    for (HitCharge &hitCharge : hitChargeVector)
        pMinuitVector->push_back(hitCharge);

    for (HitCharge &hitCharge : *pMinuitVector)
    {
        globalTotalHitWidth += hitCharge.GetHitWidth();
        globalTotalCharge += hitCharge.GetCharge();
    }

    this->GetTrackLength(hitChargeVector, globalTrackLength);
}

//---------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::PerformFits(HitChargeVector &hitChargeVector, HitChargeVector &forwardsFitPoints, HitChargeVector &backwardsFitPoints, int numberHitsToConsider, float &forwardsChiSquared, float &backwardsChiSquared, int &fitStatus1, int &fitStatus2, FitParameters &fitParameters)
{
    this->SetGlobalMinuitPreliminaries(hitChargeVector);

    //---------------------------------------------------------------------------------------------------
    //Forwards Fit

    //for (const auto &hitCharge : hitChargeVector)
    //    std::cout << hitCharge.GetChargeOverWidth() << std::endl;

    //double maxScale(globalLookupTable.GetMaxRange()/globalTrackLength);
    double particleMass(globalLookupTable.GetMass());

    LookupTable lookupTable = globalLookupTable;

    const int nParameters = 3;
    const std::string parName[nParameters]   = {"ENDENERGY", "SCALE", "FUDGE"};
    const double vstart[nParameters] = {1.1, 1.01, 1.0};
    const double step[nParameters] = {0.01, 0.01, 0.01};
    const double lowphysbound[nParameters] = {1.0, 0.5, 0.8};
    const double highphysbound[nParameters] = {100.0, 2.0, 1.2};

    int ierflg(0);

    TMinuit *pMinuit = new TMinuit(nParameters);
    pMinuit->SetPrintLevel(-1);
    pMinuit->SetFCN(GetForwardsChiSquared);

    for (int j = 0 ; j < nParameters ; ++j)
        pMinuit->mnparm(j, parName[j].c_str(), vstart[j], step[j], lowphysbound[j], highphysbound[j], ierflg);

    double arglist[2];
    arglist[0] = 1000;
    arglist[1] = 1;
    pMinuit->mnexcm("MIGRAD", arglist, 1, fitStatus1);

    double outpar[nParameters], err[nParameters];

    for (int k = 0; k < nParameters; k++)
        pMinuit->GetParameter(k, outpar[k], err[k]);

    delete pMinuit;

    //---------------------------------------------------------------------------------
    //Backwards Fit

    const int nParameters2 = 3;
    const std::string parName2[nParameters2] = {"ENDENERGY", "SCALE", "FUDGE"};
    const double vstart2[nParameters2] = {1.1, 1.01, 1.0};
    const double step2[nParameters2] = {0.01, 0.01, 0.01};
    const double lowphysbound2[nParameters2] = {1.0, 0.5, 0.8};
    const double highphysbound2[nParameters2] = {100.0, 2.0, 1.2};

    int ierflg2(0);

    TMinuit *pMinuit2 = new TMinuit(nParameters2);
    pMinuit2->SetPrintLevel(-1);
    pMinuit2->SetFCN(GetBackwardsChiSquared);

    for (int j = 0 ; j < nParameters2 ; ++j)
        pMinuit2->mnparm(j, parName2[j].c_str(), vstart2[j], step2[j], lowphysbound2[j], highphysbound2[j], ierflg2);

    double arglist2[2];
    arglist2[0] = 1000;
    arglist2[1] = 1;
    pMinuit2->mnexcm("MIGRAD", arglist2, 1, fitStatus2);

    double outpar2[nParameters2], err2[nParameters2];

    for (int k = 0; k < nParameters2; k++)
        pMinuit2->GetParameter(k, outpar2[k], err2[k]);

    delete pMinuit2;

    //--------------------------------------------------------------------------

    double f_Ee(outpar[0]), f_L(outpar[1] * globalTrackLength);
    double f_Le(GetLengthfromEnergy(lookupTable, f_Ee));
    double f_Ls = f_Le - f_L;

    double f_Es = GetEnergyfromLength(lookupTable, f_Ls);
    double f_deltaE = f_Es - f_Ee;

    double f_alpha = f_deltaE/globalTotalCharge;
    double f_beta = f_L/globalTotalHitWidth;

    double b_Ee(outpar2[0]), b_L(outpar2[1] * globalTrackLength);
    double b_Le(GetLengthfromEnergy(lookupTable, b_Ee));
    double b_Ls = b_Le - b_L;

    double b_Es = GetEnergyfromLength(lookupTable, b_Ls);
    double b_deltaE = b_Es - b_Ee;

    double b_alpha = b_deltaE/globalTotalCharge;
    double b_beta = b_L/globalTotalHitWidth;

    //--------------------------------------------------------------------------
    
    //Here we set the individual HitCharge attributes
    int nHitsConsidered(0);

    for (HitCharge &hitCharge : hitChargeVector)
    {
        double f_L_i = f_Ls + (outpar[1] * hitCharge.GetLongitudinalPosition());
        double f_E_i = GetEnergyfromLength(lookupTable, f_L_i);
        double f_dEdx_2D = outpar[2] * (f_beta/f_alpha) * BetheBloch(f_E_i, particleMass);

        double b_L_i = b_Ls + (outpar2[1] * (globalTrackLength - hitCharge.GetLongitudinalPosition()));
        double b_E_i = GetEnergyfromLength(lookupTable, b_L_i);
        double b_dEdx_2D = outpar2[2] * (b_beta/b_alpha) * BetheBloch(b_E_i, particleMass);

        double Q_fit_f(f_dEdx_2D);
        double Q_fit_b(b_dEdx_2D);

        float forwardsDelta(hitCharge.GetChargeOverWidth() - f_dEdx_2D), backwardsDelta(hitCharge.GetChargeOverWidth() - b_dEdx_2D);

        float f_sigma(std::sqrt(0.85) * std::sqrt((0.00303236 * Q_fit_f * Q_fit_f) + (0.0038681 * Q_fit_f))); //70%
        float b_sigma(std::sqrt(0.75) * std::sqrt((0.00303236 * Q_fit_b * Q_fit_b) + (0.0038681 * Q_fit_b))); //70%

        //float f_sigma(std::sqrt((0.00164585 * f_dEdx_2D * f_dEdx_2D) + (0.0201838 * f_dEdx_2D))); //80%
        //float b_sigma(std::sqrt((0.00164585 * b_dEdx_2D * b_dEdx_2D) + (0.0201838 * b_dEdx_2D))); //80%

        /*
        std::cout << "Hit Q/w: " << hitCharge.GetChargeOverWidth() << std::endl; 
        std::cout << "Hit f_dEdx_2D: " << f_dEdx_2D << std::endl; 
        std::cout << "Hit f_sigma: " << f_sigma << std::endl; 
        std::cout << "Hit b_dEdx_2D: " << b_dEdx_2D << std::endl; 
        std::cout << "Hit b_sigma: " << b_sigma << std::endl; 
        std::cout << "Hit offset: " << backwardsDelta << std::endl; 
        std::cout << "------------------" << std::endl;
        */

        float lp(hitCharge.GetLongitudinalPosition()), hw(hitCharge.GetHitWidth());
        float f_Q_fit_f(Q_fit_f), f_Q_fit_b(Q_fit_b);
        float forwardsRecoChargeOverWidth(f_Q_fit_f * hitCharge.GetHitWidth()), backwardsRecoChargeOverWidth(f_Q_fit_b * hitCharge.GetHitWidth());

        HitCharge forwardsRecoHitCharge(hitCharge.GetCaloHit(), lp, hw, forwardsRecoChargeOverWidth, f_sigma);
        forwardsFitPoints.push_back(forwardsRecoHitCharge);
        HitCharge backwardsRecoHitCharge(hitCharge.GetCaloHit(), lp, hw, backwardsRecoChargeOverWidth, b_sigma);
        backwardsFitPoints.push_back(backwardsRecoHitCharge);

        float forwardsHitChisquared((forwardsDelta * forwardsDelta)/(f_sigma * f_sigma));
        float backwardsHitChisquared((backwardsDelta * backwardsDelta)/(b_sigma * b_sigma));

        float Q_fit_forwards(Q_fit_f * hitCharge.GetHitWidth()), Q_fit_backwards(Q_fit_b * hitCharge.GetHitWidth()); 

        hitCharge.SetForwardsFitCharge(Q_fit_forwards); 
        hitCharge.SetForwardsSigma(f_sigma);
        hitCharge.SetForwardsDelta(forwardsDelta);
        hitCharge.SetForwardsChiSquared(forwardsHitChisquared);

        hitCharge.SetBackwardsFitCharge(Q_fit_backwards); 
        hitCharge.SetBackwardsSigma(b_sigma);
        hitCharge.SetBackwardsDelta(backwardsDelta);
        hitCharge.SetBackwardsChiSquared(backwardsHitChisquared);

        if (!((pMinuitVector->size() >= 2 * numberHitsToConsider) && nHitsConsidered > numberHitsToConsider && nHitsConsidered < pMinuitVector->size() - numberHitsToConsider))
        {
            forwardsChiSquared += forwardsHitChisquared;
            backwardsChiSquared += backwardsHitChisquared;
        }

        nHitsConsidered++;
    }

    if (forwardsChiSquared <= backwardsChiSquared)
    {
        float parameterZero(outpar[0]), parameterOne(outpar[1]), parameterTwo(outpar[2]);
        FitParameters bestFitParameters(parameterZero, parameterOne, parameterTwo);
        fitParameters = bestFitParameters;
    }
    else
    {
        float parameterZero(outpar2[0]), parameterOne(outpar2[1]), parameterTwo(outpar2[2]);
        FitParameters bestFitParameters(parameterZero, parameterOne, parameterTwo);
        fitParameters = bestFitParameters;
    }

    std::sort(forwardsFitPoints.begin(), forwardsFitPoints.end(), SortHitChargeVectorByRL);
    std::sort(backwardsFitPoints.begin(), backwardsFitPoints.end(), SortHitChargeVectorByRL);

    this->ClearGlobalVariables();
}

//---------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::GetCalorimetricDirection(const Cluster* pTargetClusterW, DirectionFitObject &directionFitObject)
{
    HitChargeVector hitChargeVector;
    this->FillHitChargeVector(pTargetClusterW, hitChargeVector);

    if (hitChargeVector.size() < m_minClusterCaloHits)
    {
        //std::cout << "Direction fit error: invalid cluster" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    HitChargeVector filteredHitChargeVector;
    this->TrackInnerFilter(hitChargeVector, filteredHitChargeVector);

    if (filteredHitChargeVector.size() < m_minClusterCaloHits)
    {
        //std::cout << "Direction fit error: invalid cluster" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    DirectionFitObject beforeDirectionFitObject;
    this->FitHitChargeVector(filteredHitChargeVector, beforeDirectionFitObject);


    //if (beforeDirectionFitObject.GetBackwardsChiSquaredPerHit() >= 20.0)
    //    beforeDirectionFitObject.DrawEnhancedFit();

    if (m_enableBraggPeakFilter)
    {
        this->SimpleTrackEndFilter(filteredHitChargeVector);
        this->TrackEndFilter(filteredHitChargeVector, directionFitObject, beforeDirectionFitObject);
    }

    /*
    HitChargeVector regularisedHitChargeVector;
    this->Regularise(filteredHitChargeVector, regularisedHitChargeVector);

    if (regularisedHitChargeVector.size() < m_minClusterCaloHits)
    {
        std::cout << "Direction fit error: invalid cluster" << std::endl;
        throw STATUS_CODE_FAILURE;
    }
    */

    this->FitHitChargeVector(filteredHitChargeVector, directionFitObject);

    this->TestHypothesisOne(directionFitObject);

    if (m_enableSplitting)
        this->TestHypothesisTwo(pTargetClusterW, directionFitObject);

    //directionFitObject.DrawEnhancedFit();
    //this->TestHypothesisThree(directionFitObject);
    std::cout << "Min chi squared per hit: " << directionFitObject.GetMinChiSquaredPerHit() << std::endl;

}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::TestHypothesisOne(DirectionFitObject &directionFitObject)
{
    bool likelyForwards(directionFitObject.GetDirectionEstimate() == 1 && directionFitObject.GetHitChargeVector().size() >= 400 && directionFitObject.GetForwardsChiSquaredPerHit() <= 1.5);
    bool likelyBackwards(directionFitObject.GetDirectionEstimate() == 0 && directionFitObject.GetHitChargeVector().size() <= 200 && directionFitObject.GetBackwardsChiSquaredPerHit() <= 1.5);

    if (likelyForwards || likelyBackwards)
    {
        //std::cout << "Applied Hypothesis #1 (Single Clean Particle)" << std::endl;
        directionFitObject.SetHypothesis(1); 
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::TestHypothesisTwo(const Cluster *const pCluster, DirectionFitObject &directionFitObject)
{
    if (directionFitObject.GetHypothesis() == 1 || m_enableSplitting == false)
        return;

    DirectionFitObject backwardsSplitResult, forwardsSplitResult;
    HitChargeVector filteredHitChargeVector(directionFitObject.GetHitChargeVector());

    bool splitApplied(false);
    SplitObject splitObject;
    this->ParticleSplitting(pCluster, filteredHitChargeVector, backwardsSplitResult, forwardsSplitResult, splitApplied, splitObject);

    DirectionFitObject largestFitObject(forwardsSplitResult.GetNHits() > backwardsSplitResult.GetNHits() ? forwardsSplitResult : backwardsSplitResult);

    //To create a chi squared change scatter plot
    splitObject.SetAfterNHits(largestFitObject.GetNHits());

    //std::cout << "splitObject.GetSplitCorrect(): " << splitObject.GetSplitCorrect() << std::endl;

    directionFitObject.SetSplitObject(splitObject);

    if (splitApplied)
    {
        //std::cout << "Split applied" << std::endl;
        //std::cout << "Applied Hypothesis #2 (Split Particle)" << std::endl;
        directionFitObject.SetHypothesis(2); 

        //Forwards and backwards now refer to the best fits for the forwards and backwards particles
        directionFitObject.SetForwardsFitCharges(forwardsSplitResult.GetForwardsFitCharges());
        directionFitObject.SetBackwardsFitCharges(backwardsSplitResult.GetBackwardsFitCharges());

        //Delta chi squared should still make sense for distributions, so take the likely muon
        //DirectionFitObject largestFitObject(forwardsSplitResult.GetNHits() > backwardsSplitResult.GetNHits() ? forwardsSplitResult : backwardsSplitResult);
        directionFitObject.SetForwardsChiSquared(largestFitObject.GetForwardsChiSquared());
        directionFitObject.SetBackwardsChiSquared(largestFitObject.GetBackwardsChiSquared());
        directionFitObject.SetNHits(largestFitObject.GetNHits());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackDirectionTool::IsClusterTwoParticles(const Cluster *const pCluster, TrackDirectionTool::HitChargeVector forwardsFitCharges, TrackDirectionTool::HitChargeVector backwardsFitCharges)
{
    bool isTwoParticles(false);

    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);

    int nContributingParticles(0), contributionThreshold(10);
    int nSecondaryHits(0);
    float chargeContributionThreshold(0.5);
    std::map<int, int> primaryPDGToContributionMap;

    const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

    for (const CaloHit* pCaloHit : caloHitList)
    {
        MCParticleWeightMap mcParticleWeightMap(pCaloHit->GetMCParticleWeightMap());
        for (MCParticleWeightMap::const_iterator mapIter = mcParticleWeightMap.begin(), mapIterEnd = mcParticleWeightMap.end(); mapIter != mapIterEnd; ++mapIter)
        {
            int primaryPDG(mapIter->first->GetParticleId());
            if (primaryPDG == 11 && ((pCaloHit->GetPositionVector().GetZ() > pMCParticle->GetVertex().GetZ() && pCaloHit->GetPositionVector().GetZ() < pMCParticle->GetEndpoint().GetZ()) || (pCaloHit->GetPositionVector().GetZ() < pMCParticle->GetVertex().GetZ() && pCaloHit->GetPositionVector().GetZ() > pMCParticle->GetEndpoint().GetZ())))
                primaryPDG = 13;

            float contribution(mapIter->second);

            if (primaryPDGToContributionMap.find(primaryPDG) == primaryPDGToContributionMap.end())
            {
                if (contribution > chargeContributionThreshold)
                    primaryPDGToContributionMap[primaryPDG] = 1;
            }
            else
            {
                if (contribution > chargeContributionThreshold)
                    primaryPDGToContributionMap.at(primaryPDG)++;
            }
        }
    }

    for (auto &entry : primaryPDGToContributionMap)
    {
        if (entry.first != m_targetParticlePDG)
            nSecondaryHits += entry.second;

        if (entry.second >= contributionThreshold)
            nContributingParticles++;
    }  

    float chargeContributionThreshold2(0.25);
    int backwardsNonPrimaryHits(0), backwardsPrimaryHits(0);
    for (TrackDirectionTool::HitCharge &hitCharge : backwardsFitCharges)
    {        
        const pandora::CaloHit* pCaloHit(hitCharge.GetCaloHit());
        MCParticleWeightMap mcParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

        for (MCParticleWeightMap::const_iterator mapIter = mcParticleWeightMap.begin(), mapIterEnd = mcParticleWeightMap.end(); mapIter != mapIterEnd; ++mapIter)
        {        
            int primaryPDG(mapIter->first->GetParticleId());

            float contribution(mapIter->second);
            if (contribution >= chargeContributionThreshold2 && primaryPDG != m_targetParticlePDG)
                backwardsNonPrimaryHits++;
            else if (contribution >= chargeContributionThreshold2 && primaryPDG == m_targetParticlePDG)
                backwardsPrimaryHits++;
        }    
    }    

    int forwardsNonPrimaryHits(0), forwardsPrimaryHits(0);
    for (TrackDirectionTool::HitCharge &hitCharge : forwardsFitCharges)
    {        
        const pandora::CaloHit* pCaloHit(hitCharge.GetCaloHit());
        MCParticleWeightMap mcParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

        for (MCParticleWeightMap::const_iterator mapIter = mcParticleWeightMap.begin(), mapIterEnd = mcParticleWeightMap.end(); mapIter != mapIterEnd; ++mapIter)
        {        
            int primaryPDG(mapIter->first->GetParticleId());
            float contribution(mapIter->second);

            if (contribution >= chargeContributionThreshold2 && primaryPDG != m_targetParticlePDG)
                forwardsNonPrimaryHits++;
            else if (contribution >= chargeContributionThreshold2 && primaryPDG == m_targetParticlePDG)
                forwardsPrimaryHits++;
        }    
    }    

    float forwardsImpurityFraction((static_cast<float>(forwardsNonPrimaryHits))/forwardsFitCharges.size()), backwardsImpurityFraction(((float)backwardsNonPrimaryHits)/backwardsFitCharges.size());
    //float forwardsImpurityCompleteness(nSecondaryHits != 0 ? (static_cast<float>(forwardsNonPrimaryHits))/nSecondaryHits : 0), backwardsImpurityCompleteness(nSecondaryHits != 0 ? static_cast<float>(backwardsNonPrimaryHits)/nSecondaryHits : 0);

    /*
    std::cout << "nSecondaryHits: " << nSecondaryHits << std::endl;
    std::cout << "forwardsImpurityFraction: " << forwardsImpurityFraction << std::endl;
    std::cout << "backwardsImpurityFraction: " << backwardsImpurityFraction << std::endl;
    */

    if (nSecondaryHits >= 5 || forwardsImpurityFraction >= 0.5 || backwardsImpurityFraction >= 0.5)
        isTwoParticles= true;

    return isTwoParticles;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::TestHypothesisThree(DirectionFitObject &directionFitObject)
{
    if (directionFitObject.GetHypothesis() == 1 || directionFitObject.GetHypothesis() == 2 || m_enableFragmentRemoval == false)
        return;

    HitChargeVector filteredHitChargeVector(directionFitObject.GetHitChargeVector()), fragmentlessHitChargeVector;
    float splitPosition(0.f);
    this->FragmentRemoval(filteredHitChargeVector, fragmentlessHitChargeVector, splitPosition);

    DirectionFitObject fragmentRemovalDirectionFitObject;
    this->FitHitChargeVector(fragmentlessHitChargeVector, fragmentRemovalDirectionFitObject);

    bool likelyCorrectFragmentRemoval(directionFitObject.GetDirectionEstimate() != fragmentRemovalDirectionFitObject.GetDirectionEstimate() && directionFitObject.GetMinChiSquaredPerHit() - fragmentRemovalDirectionFitObject.GetMinChiSquaredPerHit() >= 2.0);

    SplitObject frObject(filteredHitChargeVector.size(), fragmentlessHitChargeVector.size(), directionFitObject.GetMinChiSquaredPerHit(), fragmentRemovalDirectionFitObject.GetMinChiSquaredPerHit(), directionFitObject.GetMinChiSquaredPerHit() - fragmentRemovalDirectionFitObject.GetMinChiSquaredPerHit(), 0.f);

    if (m_useMCInformation)
    {
        bool splitCorrect((directionFitObject.GetDirectionEstimate() != globalTrueDirection && fragmentRemovalDirectionFitObject.GetDirectionEstimate() == globalTrueDirection) ? true : false);
        frObject.SetSplitCorrect(splitCorrect);
        frObject.SetBeforeForwardsChiSquaredPerHit(directionFitObject.GetForwardsChiSquaredPerHit());
        frObject.SetBeforeBackwardsChiSquaredPerHit(directionFitObject.GetBackwardsChiSquaredPerHit());
    }

    directionFitObject.SetFRObject(frObject); 

    if (likelyCorrectFragmentRemoval)
    {
        //std::cout << "Applied Hypothesis #3: fragment removed." << std::endl;
        directionFitObject.SetHypothesis(3); 
        directionFitObject = fragmentRemovalDirectionFitObject;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::AddToSlidingFitCache(const Cluster *const pCluster)
{
    if (m_slidingFitResultMap.find(pCluster) != m_slidingFitResultMap.end())
        return;

    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFit(pCluster, m_slidingFitWindow, slidingFitPitch);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, slidingFit)).second)
    {
        std::cout << "Sliding fit failure" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const TwoDSlidingFitResult &TrackDirectionTool::GetCachedSlidingFit(const Cluster *const pCluster) const
{
    TwoDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
    {
        std::cout << "Sliding fit retrieval failure" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::TidyUp()
{
    m_slidingFitResultMap.clear();

    globalTrackLength = (0.f);
    globalTotalHitWidth = (0.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::ClearGlobalVariables()
{
    pMinuitVector->clear();
    globalTotalCharge = 0.f;
    globalTrackLength = 0.f;
    globalTotalHitWidth = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackDirectionTool::SortHitChargeVectorByRL(HitCharge &hitCharge1, HitCharge &hitCharge2)
{
    return hitCharge1.GetLongitudinalPosition() < hitCharge2.GetLongitudinalPosition();
}

//----------------------------------------------------------------------------------------------------------------------------------

bool TrackDirectionTool::SortHitChargeVectorByChargeOverWidth(HitCharge &hitCharge1, HitCharge &hitCharge2)
{
    return hitCharge1.GetChargeOverWidth() < hitCharge2.GetChargeOverWidth();
}

//----------------------------------------------------------------------------------------------------------------------------------

bool TrackDirectionTool::SortByDistanceToNN(HitCharge &hitCharge1, HitCharge &hitCharge2)
{
    return hitCharge1.GetDistanceToNN() < hitCharge2.GetDistanceToNN();
}

//----------------------------------------------------------------------------------------------------------------------------------

bool TrackDirectionTool::SortJumpVector(JumpObject &jumpObject1, JumpObject &jumpObject2)
{
    return jumpObject1.GetJumpValue() > jumpObject2.GetJumpValue();
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::DirectionFitObject::DrawEnhancedFit()
{
    float firstLengthPosition((*(m_hitchargevector.begin())).GetLongitudinalPosition());
    float lastLengthPosition((*(std::prev(m_hitchargevector.end(), 1))).GetLongitudinalPosition());
    float trackLength(lastLengthPosition - firstLengthPosition);

    HitChargeVector hitChargeVector(m_hitchargevector), forwardsRecoHits(m_forwardsrecohits), backwardsRecoHits(m_backwardsrecohits);

    TGraphErrors *pureHits = new TGraphErrors(hitChargeVector.size());
    TGraphErrors *impureHits = new TGraphErrors(hitChargeVector.size());
    TGraphErrors *fitHits = new TGraphErrors(forwardsRecoHits.size());
    int n(1), i(1);

    std::cout << "Drawing " << hitChargeVector.size() << " hits." << std::endl;

    float minCharge(1e6), maxCharge(0);
    float minL2D(1e6), maxL2D(0);
    bool drawFit(true);

    for (HitCharge hitCharge : hitChargeVector)
    {    
        const auto pCaloHit(hitCharge.GetCaloHit());
        const MCParticleWeightMap &hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

        const auto pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
        const float weight(hitMCParticleWeightMap.at(pMCParticle));

        if (hitCharge.GetChargeOverWidth() < minCharge)
            minCharge = hitCharge.GetChargeOverWidth();

        if (hitCharge.GetChargeOverWidth() > maxCharge)
            maxCharge = hitCharge.GetChargeOverWidth();

        if (hitCharge.GetLongitudinalPosition() < minL2D)
            minL2D = hitCharge.GetLongitudinalPosition();

        if (hitCharge.GetLongitudinalPosition() > maxL2D)
            maxL2D = hitCharge.GetLongitudinalPosition();

        if (n == 0)
        {
            std::cout << "weight: " << weight << std::endl;
            std::cout << "PDG: " << pMCParticle->GetParticleId() << std::endl;
        }

        if (weight >= 0.85)
            pureHits->SetPoint(n, hitCharge.GetLongitudinalPosition(), hitCharge.GetChargeOverWidth());
        else
            impureHits->SetPoint(n, hitCharge.GetLongitudinalPosition(), hitCharge.GetChargeOverWidth());
        n++; 
    }    

    if (drawFit && (m_forwardschisquared/m_nhits < m_backwardschisquared/m_nhits || m_hypothesis == 2))
    {    
        for (HitCharge hitCharge : forwardsRecoHits)
        {    
            fitHits->SetPoint(i, hitCharge.GetLongitudinalPosition(), hitCharge.GetChargeOverWidth());
            i++; 
        }    

    }    

    if (drawFit && (m_forwardschisquared/m_nhits > m_backwardschisquared/m_nhits || m_hypothesis == 2))
    {    
        for (HitCharge hitCharge : backwardsRecoHits)
        {    
            fitHits->SetPoint(i, hitCharge.GetLongitudinalPosition(), hitCharge.GetChargeOverWidth());
            i++; 
        }    
    }    

    if (hitChargeVector.size() != 0 && forwardsRecoHits.size() != 0 && backwardsRecoHits.size() != 0)
    {
        TCanvas *canvas = new TCanvas("canvas", "canvas", 1200, 900);
        canvas->cd();

        pureHits->GetXaxis()->SetLimits(minL2D -0.05 * trackLength, maxL2D + 0.05 * trackLength);
        pureHits->GetYaxis()->SetRangeUser(0.9 * minCharge, 0.1 * minCharge +  maxCharge);
        pureHits->SetMarkerStyle(20);
        pureHits->SetMarkerSize(0.5);
        pureHits->SetMarkerColor(kBlack);

        impureHits->SetMarkerStyle(20);
        impureHits->SetMarkerSize(0.5);
        impureHits->SetMarkerColor(kRed);

        fitHits->SetMarkerStyle(20);
        fitHits->SetMarkerSize(0.5);
        fitHits->SetMarkerColor(kMagenta);

        pureHits->SetTitle(";Longitudinal Coordinate L_{2D} (cm); #tilde{Q} (MeV/cm)"); //Bethe-Bloch Theory Fit
        pureHits->Draw("AP");
        impureHits->Draw("Psame");

        if (drawFit)
            fitHits->Draw("Psame");

        //PANDORA_MONITORING_API(Pause(this->GetPandora()));
        canvas->SaveAs("fit.pdf");
        delete canvas;
    }    

    delete pureHits;
    delete impureHits;
    delete fitHits;
}

StatusCode TrackDirectionTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NumberTrackEndHits", m_numberTrackEndHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EndpointProtectionFraction", m_endpointProtectionFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableFragmentRemoval", m_enableFragmentRemoval));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableBraggPeakFilter", m_enableBraggPeakFilter));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableSplitting", m_enableSplitting));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteTable", m_writeTable));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseMCInformation", m_useMCInformation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FileName", m_lookupTableFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputFile", m_fileName));

    return STATUS_CODE_SUCCESS;
}

//----------------------------------------------------------------------------------------------------------------------------------

} // namespac lar_content
