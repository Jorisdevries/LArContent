/**
 *  @file   larpandoracontent/LArHelpers/LArSpaceChargeHelper.cc
 *
 *  @brief  Implementation of the cluster helper class.
 *
 *  $Log: $
 */

#include "LArSpaceChargeHelper.h"

using namespace pandora;

namespace lar_content
{

    std::array<Interpolator, 7U> g1_x;
    std::array<Interpolator, 7U> g2_x;
    std::array<Interpolator, 7U> g3_x;
    std::array<Interpolator, 7U> g4_x;
    std::array<Interpolator, 7U> g5_x;
      
    std::array<Interpolator, 7U> g1_y;
    std::array<Interpolator, 7U> g2_y;
    std::array<Interpolator, 7U> g3_y;
    std::array<Interpolator, 7U> g4_y;
    std::array<Interpolator, 7U> g5_y;
    std::array<Interpolator, 7U> g6_y;
      
    std::array<Interpolator, 7U> g1_z;
    std::array<Interpolator, 7U> g2_z;
    std::array<Interpolator, 7U> g3_z;
    std::array<Interpolator, 7U> g4_z;

//------------------------------------------------------------------------------------------------------------------------------------------

void LArSpaceChargeHelper::Configure(std::string fileName)
{
    TFile infile(fileName.c_str(), "READ");

    if(!infile.IsOpen())
    {
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        std::cout << "Could not find the space charge effect file '" << fileName << "'!\n";
    }


    for(int i = 0; i < 5; i++)
    {
        g1_x[i] = MakeInterpolator(infile, Form("deltaX/g1_%d",i));
        g2_x[i] = MakeInterpolator(infile, Form("deltaX/g2_%d",i));
        g3_x[i] = MakeInterpolator(infile, Form("deltaX/g3_%d",i));   
        g4_x[i] = MakeInterpolator(infile, Form("deltaX/g4_%d",i));
        g5_x[i] = MakeInterpolator(infile, Form("deltaX/g5_%d",i));
        g1_y[i] = MakeInterpolator(infile, Form("deltaY/g1_%d",i));
        g2_y[i] = MakeInterpolator(infile, Form("deltaY/g2_%d",i));
        g3_y[i] = MakeInterpolator(infile, Form("deltaY/g3_%d",i));   
        g4_y[i] = MakeInterpolator(infile, Form("deltaY/g4_%d",i));
        g5_y[i] = MakeInterpolator(infile, Form("deltaY/g5_%d",i));
        g6_y[i] = MakeInterpolator(infile, Form("deltaY/g6_%d",i));
        g1_z[i] = MakeInterpolator(infile, Form("deltaZ/g1_%d",i));
        g2_z[i] = MakeInterpolator(infile, Form("deltaZ/g2_%d",i));
        g3_z[i] = MakeInterpolator(infile, Form("deltaZ/g3_%d",i));   
        g4_z[i] = MakeInterpolator(infile, Form("deltaZ/g4_%d",i));
    }

    g1_x[5] = MakeInterpolator(infile, "deltaX/g1_5");
    g2_x[5] = MakeInterpolator(infile, "deltaX/g2_5");
    g3_x[5] = MakeInterpolator(infile, "deltaX/g3_5");   
    g4_x[5] = MakeInterpolator(infile, "deltaX/g4_5");
    g5_x[5] = MakeInterpolator(infile, "deltaX/g5_5");
    g1_y[5] = MakeInterpolator(infile, "deltaY/g1_5");
    g2_y[5] = MakeInterpolator(infile, "deltaY/g2_5");
    g3_y[5] = MakeInterpolator(infile, "deltaY/g3_5");   
    g4_y[5] = MakeInterpolator(infile, "deltaY/g4_5");
    g5_y[5] = MakeInterpolator(infile, "deltaY/g5_5");
    g6_y[5] = MakeInterpolator(infile, "deltaY/g6_5");

    g1_x[6] = MakeInterpolator(infile, "deltaX/g1_6");
    g2_x[6] = MakeInterpolator(infile, "deltaX/g2_6");
    g3_x[6] = MakeInterpolator(infile, "deltaX/g3_6");
    g4_x[6] = MakeInterpolator(infile, "deltaX/g4_6");
    g5_x[6] = MakeInterpolator(infile, "deltaX/g5_6");

    infile.Close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(CartesianVector &positionVector)
{
    return positionVector + GetPositionOffset(positionVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArSpaceChargeHelper::GetPositionOffset(CartesianVector &positionVector)
{
    CartesianVector transformedPosition(TransformPosition(positionVector));
    CartesianVector positionOffset(GetPositionOffsetX(transformedPosition), GetPositionOffsetY(transformedPosition), GetPositionOffsetZ(transformedPosition));
    return positionOffset;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArSpaceChargeHelper::GetPositionOffsetX(CartesianVector &positionVector)
{
  double parA[5][7];
  double const zValNew = positionVector.GetZ();

  for(size_t j = 0; j < 7; j++) 
  {
    parA[0][j] = g1_x[j].Eval(zValNew);
    parA[1][j] = g2_x[j].Eval(zValNew);
    parA[2][j] = g3_x[j].Eval(zValNew);
    parA[3][j] = g4_x[j].Eval(zValNew);
    parA[4][j] = g5_x[j].Eval(zValNew);
  }

  f1_x_poly_t f1_x;
  f2_x_poly_t f2_x;
  f3_x_poly_t f3_x;
  f4_x_poly_t f4_x;
  f5_x_poly_t f5_x;
  
  f1_x.SetParameters(parA[0]);
  f2_x.SetParameters(parA[1]);
  f3_x.SetParameters(parA[2]);
  f4_x.SetParameters(parA[3]);
  f5_x.SetParameters(parA[4]);

  double const aValNew = positionVector.GetY();

  double const parB[] = {
    f1_x.Eval(aValNew),
    f2_x.Eval(aValNew),
    f3_x.Eval(aValNew),
    f4_x.Eval(aValNew),
    f5_x.Eval(aValNew)
  };
  
  double const bValNew = positionVector.GetX();
  return 100.0*fFinal_x_poly_t::Eval(bValNew, parB);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArSpaceChargeHelper::GetPositionOffsetY(CartesianVector &positionVector)
{
  double parA[6][6];
  double const zValNew = positionVector.GetZ();

  for(size_t j = 0; j < 6; j++)
  {
    parA[0][j] = g1_y[j].Eval(zValNew);
    parA[1][j] = g2_y[j].Eval(zValNew);
    parA[2][j] = g3_y[j].Eval(zValNew);
    parA[3][j] = g4_y[j].Eval(zValNew);
    parA[4][j] = g5_y[j].Eval(zValNew);
    parA[5][j] = g6_y[j].Eval(zValNew);
  }
  
  f1_y_poly_t f1_y;
  f2_y_poly_t f2_y;
  f3_y_poly_t f3_y;
  f4_y_poly_t f4_y;
  f5_y_poly_t f5_y;
  f6_y_poly_t f6_y;
  
  f1_y.SetParameters(parA[0]);
  f2_y.SetParameters(parA[1]);
  f3_y.SetParameters(parA[2]);
  f4_y.SetParameters(parA[3]);
  f5_y.SetParameters(parA[4]);
  f6_y.SetParameters(parA[5]);
  
  double const aValNew = positionVector.GetX();

  double const parB[] = {
    f1_y.Eval(aValNew),
    f2_y.Eval(aValNew),
    f3_y.Eval(aValNew),
    f4_y.Eval(aValNew),
    f5_y.Eval(aValNew),
    f6_y.Eval(aValNew)
  };
  
  double const bValNew = positionVector.GetY();
  return 100.0*fFinal_y_poly_t::Eval(bValNew, parB);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArSpaceChargeHelper::GetPositionOffsetZ(CartesianVector &positionVector)
{
  double parA[4][5];
  double const zValNew = positionVector.GetZ();

  for(size_t j = 0; j < 5; j++)
  {
    parA[0][j] = g1_z[j].Eval(zValNew);
    parA[1][j] = g2_z[j].Eval(zValNew);
    parA[2][j] = g3_z[j].Eval(zValNew);
    parA[3][j] = g4_z[j].Eval(zValNew);
  }

  f1_z_poly_t f1_z;
  f2_z_poly_t f2_z;
  f3_z_poly_t f3_z;
  f4_z_poly_t f4_z;
  
  f1_z.SetParameters(parA[0]);
  f2_z.SetParameters(parA[1]);
  f3_z.SetParameters(parA[2]);
  f4_z.SetParameters(parA[3]);
  
  double const aValNew = positionVector.GetY();
  double const parB[] = {
    f1_z.Eval(aValNew),
    f2_z.Eval(aValNew),
    f3_z.Eval(aValNew),
    f4_z.Eval(aValNew)
  };
  
  double const bValNew = positionVector.GetX();
  return 100.0*fFinal_z_poly_t::Eval(bValNew, parB);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArSpaceChargeHelper::TransformPosition(CartesianVector &positionVector)
{
    CartesianVector transformedPosition(TransformX(positionVector.GetX()), TransformY(positionVector.GetY()), TransformZ(positionVector.GetZ()));
    return transformedPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArSpaceChargeHelper::TransformX(float xPosition)
{
    return 1.25 - (2.50/2.56)*(xPosition/100.0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArSpaceChargeHelper::TransformY(float yPosition)
{
    return (2.50/2.33)*((yPosition/100.0)+1.165) - 1.25;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArSpaceChargeHelper::TransformZ(float zPosition)
{
    return (10.0/10.37)*(zPosition/100.0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArSpaceChargeHelper::IsInsideDetectorBoundaries(CartesianVector &positionVector)
{
  return !((positionVector.GetX() <    0.0) || (positionVector.GetX() >  260.0)
        || (positionVector.GetY() < -120.0) || (positionVector.GetY() >  120.0)
        || (positionVector.GetZ() <    0.0) || (positionVector.GetZ() > 1050.0));
}

//------------------------------------------------------------------------------------------------------------------------------------------

Interpolator LArSpaceChargeHelper::MakeInterpolator(TFile& file, char const* graphName)
{
    SortedTGraphFromFile graph(file, graphName);
    return Interpolator(*graph);
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
