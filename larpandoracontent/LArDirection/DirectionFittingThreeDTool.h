/**
 *  @file   larpandoracontent/LArVertex/DirectionFittingThreeDTool.h
 *
 *  @brief  Header file for the candidate vertex creation AlgorithmTool class.
 *
 *  $Log: $
 */
#ifndef LAR_DIRECTION_FITTING_THREED_TOOL
#define LAR_DIRECTION_FITTING_THREED_TOOL 1

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArPfoObjects.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "Pandora/AlgorithmTool.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"

#include <utility>
#include <algorithm>
#include <deque>
#include <unordered_map>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <numeric>

namespace lar_content
{

class DirectionFittingThreeDTool : public pandora::AlgorithmTool
{
public:

    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

    DirectionFittingThreeDTool();

    ~DirectionFittingThreeDTool();

    //-----------------------------------------------------------------------------------------------

    class HitObject
    {
    public:

        HitObject(const pandora::CaloHit* caloHit, float &longitudinalPosition, float &hitEnergy, float &hitWidth, float &segmentLength, float &dQdx, float &dEdx);

        const pandora::CaloHit* GetCaloHit() const;
        void SetLongitudinalPosition(float longitudinalPosition);
        float GetLongitudinalPosition() const;
        float GetEnergy() const;
        float GetWidth() const;
        void SetSegmentLength(float segmentLength);
        float GetSegmentLength() const;
        float GetdQdx() const;
        float GetdEdx() const;

        void SetDistanceToNN(float &distance);
        float GetDistanceToNN() const;

        void SetForwardsFitdEdx(float forwardsFitdEdx);
        float GetForwardsFitdEdx() const;

        void SetBackwardsFitdEdx(float backwardsFitdEdx);
        float GetBackwardsFitdEdx() const;

    private:
        const pandora::CaloHit*                    m_calohit;
        float                                      m_longitudinalposition;
        float                                      m_hitenergy;
        float                                      m_hitwidth;
        float                                      m_segmentlength;
        float                                      m_dqdx;        
        float                                      m_dedx;        

        float                                      m_distancetonearestneighbour;
        float                                      m_forwardsfitenergy;
        float                                      m_backwardsfitenergy;
    };

    //-----------------------------------------------------------------------------------------------

    typedef std::vector<HitObject> HitObjectVector;

    class FitParameters
    {
        public:
        
            FitParameters();
            FitParameters(float parameterZero, float parameterOne, float parameterTwo);

            float GetParameterZero();
            float GetParameterOne();
            float GetParameterTwo();

        private:
            float   m_parameterzero;
            float   m_parameterone;
            float   m_parametertwo;
    };

    //-----------------------------------------------------------------------------------------------

    class LookupTable
    {
    public:

        LookupTable();
        LookupTable(double &initialEnergy, double &binWidth);

        std::map<int, double> GetMap();
        void SetMap(std::map<int, double> &map);

        double GetInitialEnergy();
        void SetInitialEnergy(double &initialEnergy);

        double GetMaxEnergy();
        void SetMaxEnergy(double &initialEnergy);

        int GetMaxBin();
        void SetMaxBin(int maxBin);

        double GetBinWidth();
        void SetBinWidth(double &binWidth);

        void SetMaxRange(double maxRange);
        double GetMaxRange();

        void SetMass(double &maxRange);
        double GetMass();

    private:
        std::map<int, double>                       m_map;
        std::map<double, int>                       m_reversemap;
        double                                      m_binwidth;
        double                                      m_initialenergy;
        double                                      m_maxenergy;
        int                                         m_maxbin;
        double                                      m_maxrange;
        double                                      m_mass;
    };

    //-----------------------------------------------------------------------------------------------

    class DirectionFitObject
    {
    public:
        DirectionFitObject();
        DirectionFitObject(HitObjectVector &hitObjectVector, int &numberHits, float &forwardsChiSquared, float &backwardsChiSquared, FitParameters &fitParameters, int &fitStatus);

        DirectionFittingThreeDTool::HitObjectVector GetHitObjectVector();

        float GetForwardsChiSquared();
        float GetBackwardsChiSquared();
        float GetForwardsChiSquaredPerHit();
        float GetBackwardsChiSquaredPerHit();

        int GetNHits();

        float GetMinChiSquared();
        float GetMinChiSquaredPerHit();
        float GetDeltaChiSquaredPerHit();

        void SetMCParent(int mcParentPdg);
        int GetMCParent();

        void SetContained(bool isContained);
        bool GetContained();

        void DrawFit();

        void SetFitParameters(FitParameters &fitParameters);
        FitParameters GetFitParameters();
        
        float GetFitMass();
        int GetFitStatus();

    private:
        HitObjectVector     m_hitobjectvector;

        int                 m_nhits;
        int                 m_mcparent;
        bool                m_contained;

        float               m_forwardschisquared;
        float               m_backwardschisquared;

        FitParameters       m_fitparameters;
        int                 m_fitstatus;
    };

    //-----------------------------------------------------------------------------------------------

    DirectionFittingThreeDTool::DirectionFitObject GetPfoDirection(const pandora::ParticleFlowObject *const pPfo);

    private:

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void FillHitObjectVector(const pandora::ParticleFlowObject *const pPfo, HitObjectVector &hitObjectVector);

    void GetTrackLength(HitObjectVector &hitEnergyVector, float &trackLength);

    float GetDriftCoordinateCorrection(const pandora::CartesianVector &positionVector);

    float GetYZCoordinateCorrection(const pandora::CartesianVector &positionVector);

    float CalculateModBoxdEdx(float &dQdx);

    void FilterHitObjectVector(HitObjectVector &hitObjectVector, HitObjectVector &filteredHitObjectVector);

    void DrawHitObjectVector(HitObjectVector &hitObjectVector);

    void SelectEndpoint(HitObjectVector &hitObjectVector, HitObjectVector &filteredHitObjectVector);

    void SetNearestNeighbourValues(HitObjectVector &innerHitObjectVector, int &nNeighboursToConsider);

    void FitHitObjectVector(HitObjectVector &hitObjectVector, DirectionFittingThreeDTool::DirectionFitObject &fitResult);

    void PerformFits(HitObjectVector &hitObjectVector, float &forwardsChiSquared, float &backwardsChiSquared, FitParameters &fitParameters, int &fitStatus);

    void SetFinalFitValues(HitObjectVector &hitObjectVector, float &forwardsChiSquared, float &backwardsChiSquared, double (&fitParameters)[2], bool forwards);

    void SetMCInformation(const pandora::ParticleFlowObject *const pPfo, DirectionFittingThreeDTool::DirectionFitObject &fitResult);

    bool IsParticleContained(const pandora::MCParticle* pMCParticle);

    static bool SortByLongitudinalPosition(HitObject &hitObject1, HitObject &hitObject2);
    static bool SortByDistanceToNN(HitObject &hitObject1, HitObject &hitObject2);
    static bool SortByEnergy(HitObject &hitObject1, HitObject &hitObject2);

    //-----------------------------------------------------------------------------------------------

    float                   m_endpointRange;

    unsigned int            m_slidingFitWindow;                
    unsigned int            m_minClusterCaloHits;             
    float                   m_minClusterLength;          

    double                  m_tableInitialEnergy;
    double                  m_tableStepSize;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *DirectionFittingThreeDTool::Factory::CreateAlgorithmTool() const
{
    return new DirectionFittingThreeDTool();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline DirectionFittingThreeDTool::HitObject::HitObject(const pandora::CaloHit* caloHit, float &longitudinalPosition, float &hitEnergy, float &hitWidth, float &segmentLength, float &dQdx, float &dEdx) : 
    m_calohit(caloHit),
    m_longitudinalposition(longitudinalPosition),
    m_hitenergy(hitEnergy),
    m_hitwidth(hitWidth),
    m_segmentlength(segmentLength),
    m_dqdx(dQdx),
    m_dedx(dEdx),
    m_distancetonearestneighbour(0.f),
    m_forwardsfitenergy(0.f),
    m_backwardsfitenergy(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHit* DirectionFittingThreeDTool::HitObject::GetCaloHit() const
{
    return m_calohit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::HitObject::SetLongitudinalPosition(float longitudinalPosition) 
{
    m_longitudinalposition = longitudinalPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::HitObject::GetLongitudinalPosition() const
{
    return m_longitudinalposition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::HitObject::GetEnergy() const
{
    return m_hitenergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::HitObject::GetWidth() const
{
    return m_hitwidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::HitObject::SetSegmentLength(float segmentLength) 
{
    m_segmentlength = segmentLength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::HitObject::GetSegmentLength() const
{
    return m_segmentlength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::HitObject::GetdQdx() const
{
    return m_dqdx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::HitObject::GetdEdx() const
{
    return m_dedx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::HitObject::SetDistanceToNN(float &distance)
{
    m_distancetonearestneighbour = distance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::HitObject::GetDistanceToNN() const
{
    return m_distancetonearestneighbour;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::HitObject::SetForwardsFitdEdx(float forwardsFitdEdx) 
{
    m_forwardsfitenergy = forwardsFitdEdx;
}
//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::HitObject::GetForwardsFitdEdx() const 
{
    return m_forwardsfitenergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::HitObject::SetBackwardsFitdEdx(float backwardsFitdEdx)
{
    m_backwardsfitenergy = backwardsFitdEdx;
}
//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::HitObject::GetBackwardsFitdEdx() const 
{
    return m_backwardsfitenergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline DirectionFittingThreeDTool::FitParameters::FitParameters() :
    m_parameterzero(0.f),
    m_parameterone(0.f),
    m_parametertwo(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline DirectionFittingThreeDTool::FitParameters::FitParameters(float parameterZero, float parameterOne, float parameterTwo) :
    m_parameterzero(parameterZero),
    m_parameterone(parameterOne),
    m_parametertwo(parameterTwo)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::FitParameters::GetParameterZero()
{
    return m_parameterzero;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::FitParameters::GetParameterOne()
{
    return m_parameterone;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::FitParameters::GetParameterTwo()
{
    return m_parametertwo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline DirectionFittingThreeDTool::LookupTable::LookupTable()
{
    std::map<int, double> emptyMap;

    m_map = emptyMap;
    m_binwidth = 0.f;
    m_initialenergy = 0.f;
    m_maxenergy = 0.f;
    m_maxbin = 0;
    m_maxrange = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline DirectionFittingThreeDTool::LookupTable::LookupTable(double &initialEnergy, double &binWidth)
{
    std::map<int, double> emptyMap;

    m_map = emptyMap;
    m_binwidth = binWidth;
    m_initialenergy = initialEnergy;
    m_maxenergy = 0.f;
    m_maxbin = 0;
    m_maxrange = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::map<int, double> DirectionFittingThreeDTool::LookupTable::GetMap()
{
    return m_map;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::LookupTable::SetMap(std::map<int, double> &map)
{
    m_map = map;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double DirectionFittingThreeDTool::LookupTable::GetInitialEnergy()
{
    return m_initialenergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::LookupTable::SetInitialEnergy(double &initialEnergy)
{
    m_initialenergy = initialEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double DirectionFittingThreeDTool::LookupTable::GetMaxEnergy()
{
    return m_maxenergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::LookupTable::SetMaxEnergy(double &maxEnergy)
{
    m_maxenergy = maxEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int DirectionFittingThreeDTool::LookupTable::GetMaxBin()
{
    return m_maxbin;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::LookupTable::SetMaxBin(int maxBin)
{
    m_maxbin = maxBin;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double DirectionFittingThreeDTool::LookupTable::GetBinWidth()
{
    return m_binwidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::LookupTable::SetBinWidth(double &binWidth)
{
    m_binwidth = binWidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::LookupTable::SetMaxRange(double maxRange)
{
    m_maxrange = maxRange;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double DirectionFittingThreeDTool::LookupTable::GetMaxRange()
{
    return m_maxrange;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::LookupTable::SetMass(double &mass)
{
    m_mass = mass;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double DirectionFittingThreeDTool::LookupTable::GetMass()
{
    return m_mass;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline DirectionFittingThreeDTool::DirectionFitObject::DirectionFitObject() 
{
    HitObjectVector emptyVector;
    FitParameters emptyFitParameters;

    m_hitobjectvector = emptyVector;
    m_nhits = 0;
    m_mcparent = 0;
    m_contained = false;
    m_forwardschisquared = 0.f;
    m_backwardschisquared = 0.f;
    m_fitparameters = emptyFitParameters;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline DirectionFittingThreeDTool::DirectionFitObject::DirectionFitObject(HitObjectVector &hitObjectVector, int &numberHits, float &forwardsChiSquared, float &backwardsChiSquared, FitParameters &fitParameters, int &fitStatus) :
    m_hitobjectvector(hitObjectVector),
    m_nhits(numberHits),
    m_forwardschisquared(forwardsChiSquared),
    m_backwardschisquared(backwardsChiSquared),
    m_fitparameters(fitParameters),
    m_fitstatus(fitStatus)
{
    m_mcparent = 0;
    m_contained = false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline DirectionFittingThreeDTool::HitObjectVector DirectionFittingThreeDTool::DirectionFitObject::GetHitObjectVector()
{
    return m_hitobjectvector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::DirectionFitObject::GetForwardsChiSquared()
{
    return m_forwardschisquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::DirectionFitObject::GetBackwardsChiSquared()
{
    return m_backwardschisquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::DirectionFitObject::GetForwardsChiSquaredPerHit()
{
    return m_forwardschisquared/m_nhits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::DirectionFitObject::GetBackwardsChiSquaredPerHit()
{
    return m_backwardschisquared/m_nhits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int DirectionFittingThreeDTool::DirectionFitObject::GetNHits()
{
    return m_hitobjectvector.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::DirectionFitObject::GetMinChiSquared()
{
    return std::min(m_forwardschisquared, m_backwardschisquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::DirectionFitObject::GetMinChiSquaredPerHit()
{
    return (m_nhits != 0 ? std::min(m_forwardschisquared, m_backwardschisquared)/m_nhits : std::numeric_limits<float>::max());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::DirectionFitObject::GetDeltaChiSquaredPerHit()
{
    return (m_nhits != 0 ? (m_forwardschisquared - m_backwardschisquared)/m_nhits : std::numeric_limits<float>::max());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::DirectionFitObject::SetMCParent(int mcParentPdg)
{
   m_mcparent = mcParentPdg; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int DirectionFittingThreeDTool::DirectionFitObject::GetMCParent()
{
    return m_mcparent; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::DirectionFitObject::SetContained(bool isContained)
{
   m_contained = isContained; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool DirectionFittingThreeDTool::DirectionFitObject::GetContained()
{
    return m_contained; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::DirectionFitObject::DrawFit()
{
    float firstLengthPosition(m_hitobjectvector.front().GetLongitudinalPosition());
    float lastLengthPosition(m_hitobjectvector.back().GetLongitudinalPosition());
    float trackLength(lastLengthPosition - firstLengthPosition);
    float minEnergy(1e6), maxEnergy(0.f);

    HitObjectVector hitObjectVector(m_hitobjectvector);

    TGraphErrors *Hits = new TGraphErrors(hitObjectVector.size());
    TGraphErrors *fitHits = new TGraphErrors(hitObjectVector.size());

    int n(0);
    
    for (HitObject hitObject : hitObjectVector)
    {    
        Hits->SetPoint(n, hitObject.GetLongitudinalPosition(), hitObject.GetdEdx());
    
        if (hitObject.GetdEdx() < minEnergy)
            minEnergy = hitObject.GetdEdx();

        if (hitObject.GetdEdx() > maxEnergy)
            maxEnergy = hitObject.GetdEdx();

        if (m_forwardschisquared/m_nhits < m_backwardschisquared/m_nhits)
        {
            fitHits->SetPoint(n, hitObject.GetLongitudinalPosition(), hitObject.GetForwardsFitdEdx());

            if (hitObject.GetForwardsFitdEdx() < minEnergy)
                minEnergy = hitObject.GetForwardsFitdEdx();

            if (hitObject.GetForwardsFitdEdx() > maxEnergy)
                maxEnergy = hitObject.GetForwardsFitdEdx();
        }
        else
        {
            fitHits->SetPoint(n, hitObject.GetLongitudinalPosition(), hitObject.GetBackwardsFitdEdx());

            if (hitObject.GetBackwardsFitdEdx() < minEnergy)
                minEnergy = hitObject.GetBackwardsFitdEdx();

            if (hitObject.GetBackwardsFitdEdx() > maxEnergy)
                maxEnergy = hitObject.GetBackwardsFitdEdx();
        }

        n++; 
    }    

    if (hitObjectVector.size() != 0) 
    {    
        TCanvas *canvas = new TCanvas("canvas", "canvas", 900, 600);
        canvas->cd();
         
        Hits->GetXaxis()->SetLimits(m_hitobjectvector.front().GetLongitudinalPosition() - 0.05 * trackLength, m_hitobjectvector.back().GetLongitudinalPosition() + 0.05 * trackLength);
        Hits->GetYaxis()->SetRangeUser(minEnergy - 0.25, maxEnergy + 0.25);
        Hits->SetMarkerStyle(20);
        Hits->SetMarkerSize(0.5);
        Hits->SetMarkerColor(kBlack); 
        
        fitHits->SetMarkerStyle(20);
        fitHits->SetMarkerSize(0.5);
        fitHits->SetMarkerColor(kMagenta);

        Hits->SetTitle(";Longitudinal Coordinate L_{3D} (cm); Hit Energy (MeV)");
        Hits->Draw("AP");
        fitHits->Draw("Psame");

        //PANDORA_MONITORING_API(Pause(this->GetPandora()));
        canvas->SaveAs("fit.png");
        delete canvas;
    } 

    delete Hits;
    delete fitHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DirectionFittingThreeDTool::DirectionFitObject::SetFitParameters(FitParameters &fitParameters)
{
    m_fitparameters = fitParameters;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline DirectionFittingThreeDTool::FitParameters DirectionFittingThreeDTool::DirectionFitObject::GetFitParameters()
{
    return m_fitparameters;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DirectionFittingThreeDTool::DirectionFitObject::GetFitMass()
{
    return m_fitparameters.GetParameterOne();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int DirectionFittingThreeDTool::DirectionFitObject::GetFitStatus()
{
    return m_fitparameters.GetParameterOne();
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // #ifndef LAR_DIRECTION_FITTING_THREED_TOOL
