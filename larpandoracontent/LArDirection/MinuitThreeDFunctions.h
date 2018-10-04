/**
 *  @file   larpandoracontent/LArVertex/MinuitFunctions.h
 *
 *  @brief  Header file for Minuit implementation.
 *
 *  $Log: $
 */

#ifndef TOOL_MINUIT_FUNCTIONS_H
#define TOOL_MINUIT_FUNCTIONS_H 1

//----------------------------------------------------------------------------------------------------------------------------------

//PHYSICAL CONSTANTS

const double K = 0.307075; // constant K in MeV cm mol^-1
const double z = 1; // charge in e
const double Z = 18; // Atomic number Z
const double A = 39.948; // Atomic mass in g mol-1
const double m_e = 0.511; // Mass of electron in MeV
const double rho = 1.396; // Density of material in g cm^-3 (here: argon density)
const double I = 0.000188; // Ionisation energy in MeV

const double C = 5.2146;
const double a = 0.19559;
const double m = 3.0;
const double X1 = 3.0;
const double X0 = 0.2000;
const double delta0 = 0.0;

//----------------------------------------------------------------------------------------------------------------------------------

void BinHitObjectVector(lar_content::DirectionFittingThreeDTool::HitObjectVector &hitObjectVector, lar_content::DirectionFittingThreeDTool::HitObjectVector &binnedHitObjectVector)
{   
    float binSize = (hitObjectVector.size() > 50 ? (0.5 + (hitObjectVector.size() - 50) * 2.5/300) : 0.5);
    float trackLength(hitObjectVector.back().GetLongitudinalPosition() - hitObjectVector.front().GetLongitudinalPosition());
    
    for (float i = binSize; i <= trackLength; i += binSize)
    {   
        int nHitsBin(0);
        float meanBinPosition(0.f), meanBinCharge(0.f), meanBinWidth(0.f), meanBinSegmentLength(0.f);;
        
        for (lar_content::DirectionFittingThreeDTool::HitObject &hitObject : hitObjectVector)
        {   
            if (!(hitObject.GetLongitudinalPosition() < i && hitObject.GetLongitudinalPosition() >= (i - binSize)))
                continue; 

            if (hitObject.GetLongitudinalPosition() > i)
                break;
            
            meanBinPosition += hitObject.GetLongitudinalPosition();
            meanBinCharge += hitObject.GetHitEnergy();
            meanBinWidth += hitObject.GetHitWidth();
            meanBinSegmentLength += hitObject.GetSegmentLength();

            ++nHitsBin;
        }
        
        if (nHitsBin == 0)
            continue;
        
        meanBinPosition /= nHitsBin;
        meanBinCharge /= nHitsBin;
        meanBinWidth /= nHitsBin;
        meanBinSegmentLength /= nHitsBin;
        
        lar_content::DirectionFittingThreeDTool::HitObject binnedHitObject(NULL, meanBinPosition, meanBinCharge, meanBinWidth, meanBinSegmentLength);
        binnedHitObjectVector.push_back(binnedHitObject);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------

double DensityCorrection3D(double &T, double &M)
{
    double p = std::sqrt((T*T) + 2*T*M);
    double gamma = std::sqrt(1 + ((p/M) * (p/M)));
    double beta = std::sqrt(1 - 1 / (gamma*gamma));
    double X = std::log10(beta*gamma);

    if (X < X0)
        return delta0;
    else if ((X > X0) && (X < X1))
        return 2 * X * std::log(10) - C + (a * (std::pow((X1 - X), m)));
    else
        return 2 * X * std::log(10) + C;
}

//----------------------------------------------------------------------------------------------------------------------------------

double BetheBloch3D(double &T, double &M)
{
    double p(std::sqrt((T*T) + 2*T*M));
    double gamma(std::sqrt(1 + ((p/M) * (p/M))));
    double beta(std::sqrt(1 - 1 / (gamma*gamma)));

    double T_max(2 * m_e * (beta*gamma*beta*gamma) / (1 + 2 * gamma * m_e / M + ((m_e/M) * (m_e/M))));
    //return rho * ((K * z * z * Z) / A) * (0.5*std::log(2 * m_e * T_max * (beta*gamma*beta*gamma) / (I*I) ) - (beta*beta) - (0.5*DensityCorrection3D(p, M))) / (beta*beta); //in MeV/cm

    double W_cut(0.9 * T_max);
    return rho * ((K * z * z * Z) / A) * (0.5*std::log(2 * m_e * W_cut * (beta*gamma*beta*gamma) / (I*I) ) - ((beta*beta*0.5) * (1 + (W_cut/T_max))) - (0.5*DensityCorrection3D(p, M))) / (beta*beta); 
}

//----------------------------------------------------------------------------------------------------------------------------------

void FillLookupTable3D(lar_content::DirectionFittingThreeDTool::LookupTable &lookupTable, double M, double L_offset)
{
    std::map<int, double> lookupMap, finalLookupMap;
    double endpointRange(globalEndpointRange);

    double currentEnergy(lookupTable.GetInitialEnergy()), binWidth(lookupTable.GetBinWidth());
    int maxBin(0);

    for (double n = 0; n < 100000; ++n)
    {
        double currentdEdx = BetheBloch3D(currentEnergy, M);

        if ((currentdEdx * binWidth) >= currentEnergy || currentdEdx < 0.0)
        {
            double maxRange = (n * binWidth) + (currentEnergy/currentdEdx);
            lookupTable.SetMaxRange(maxRange);
            maxBin = n;

            lookupMap.insert(std::pair<int, double>(n, 0.0));

            break;
        }
        else
        {
            lookupMap.insert(std::pair<int, double>(n, currentEnergy));
        }

        currentEnergy -= (currentdEdx * binWidth);
    }

    int i(0);
    for (const auto &entry : lookupMap)
    {
        double residualRange((maxBin - entry.first) * binWidth);

        if (residualRange <= (endpointRange + L_offset) && residualRange >= L_offset)
        {
            finalLookupMap.emplace(i, entry.second);
            ++i;
        }
    }

    lookupTable.SetMap(finalLookupMap);
    lookupTable.SetMass(M);
    lookupTable.SetMaxEnergy(finalLookupMap.at(1));
    lookupTable.SetMaxBin(i - 1);

    /*
    TGraphErrors *Hits = new TGraphErrors(lookupMap.size());

    int n(0);
    
    for (const auto &entry : finalLookupMap)
    {    
        Hits->SetPoint(n, entry.first * binWidth, entry.second);
        ++n; 
    }  

    TCanvas *canvas = new TCanvas("canvas", "canvas", 900, 600);
    canvas->cd();

    Hits->SetMarkerStyle(20);
    Hits->SetMarkerSize(0.5);
    Hits->SetMarkerColor(kBlack); 

    Hits->SetTitle("Lookup Table;Longitudinal Coordinate L_{3D} (cm); Energy (MeV)"); //Bethe-Bloch Theory Fit
    Hits->Draw("AP");

    canvas->SaveAs("lookuptable.png");
    delete canvas;
    */
}

//----------------------------------------------------------------------------------------------------------------------------------

double GetForwardsFitEnergy(lar_content::DirectionFittingThreeDTool::HitObject &hitObject)
{
    int binNumber(std::floor(hitObject.GetLongitudinalPosition()/globalLookupTable3D.GetBinWidth()));
    int binShift((0.5 * hitObject.GetSegmentLength())/globalLookupTable3D.GetBinWidth());

    int leftBinNumber(std::min(std::max(binNumber - binShift, 0), globalLookupTable3D.GetMaxBin())); 
    int rightBinNumber(std::max(std::min(globalLookupTable3D.GetMaxBin(), binNumber + binShift), 0));

    double leftEnergy(globalLookupTable3D.GetMap().at(leftBinNumber)), rightEnergy(globalLookupTable3D.GetMap().at(rightBinNumber));

    return (leftEnergy - rightEnergy); 
}

//----------------------------------------------------------------------------------------------------------------------------------

double GetBackwardsFitEnergy(lar_content::DirectionFittingThreeDTool::HitObject &hitObject)
{
    int binNumber(std::floor((globalTrackLength3D - hitObject.GetLongitudinalPosition())/globalLookupTable3D.GetBinWidth()));
    int binShift((0.5 * hitObject.GetSegmentLength())/globalLookupTable3D.GetBinWidth());

    int leftBinNumber(std::min(std::max(binNumber - binShift, 0), globalLookupTable3D.GetMaxBin())); 
    int rightBinNumber(std::max(std::min(globalLookupTable3D.GetMaxBin(), binNumber + binShift), 0));

    double leftEnergy(globalLookupTable3D.GetMap().at(leftBinNumber)), rightEnergy(globalLookupTable3D.GetMap().at(rightBinNumber));

    return (leftEnergy - rightEnergy); 
}

//----------------------------------------------------------------------------------------------------------------------------------

void GetForwardsChiSquared3D(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t )
{
    double L_offset(par[0]);
    double M = globalLookupTable3D.GetMass();

    if (!globalFixedMass)
    {
        M = par[1];
        FillLookupTable3D(globalLookupTable3D, M, L_offset);
    }

    double chiSquared(0.0);

    for (lar_content::DirectionFittingThreeDTool::HitObject &hitObject : *pMinuitVector3D)
    {
        double hitEnergy(hitObject.GetHitEnergy());
        double fitEnergy(GetForwardsFitEnergy(hitObject));
        double hitUncertainty(1.0);

        chiSquared += ((hitEnergy - fitEnergy) * (hitEnergy - fitEnergy) )/(hitUncertainty * hitUncertainty);
    }

    f = chiSquared;
}

//----------------------------------------------------------------------------------------------------------------------------------

void GetBackwardsChiSquared3D(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
    double L_offset(par[0]);
    double M = globalLookupTable3D.GetMass();

    if (!globalFixedMass)
    {
        M = par[1];
        FillLookupTable3D(globalLookupTable3D, M, L_offset);
    }

    double chiSquared(0.0);

    for (lar_content::DirectionFittingThreeDTool::HitObject &hitObject : *pMinuitVector3D)
    {
        double hitEnergy(hitObject.GetHitEnergy());
        double fitEnergy(GetBackwardsFitEnergy(hitObject));
        double hitUncertainty(1.0);

        chiSquared += ((hitEnergy - fitEnergy) * (hitEnergy - fitEnergy) )/(hitUncertainty * hitUncertainty);
    }

    f = chiSquared;
}

//----------------------------------------------------------------------------------------------------------------------------------

#endif
