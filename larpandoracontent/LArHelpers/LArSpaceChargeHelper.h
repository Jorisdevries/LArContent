/**
 *  @file   larpandoracontent/LArHelpers/LArSpaceChargeHelper.h
 *
 *  @brief  Header file for the cluster helper class.
 *
 *  $Log: $
 */
#ifndef LAR_SPACE_CHARGE_HELPER
#define LAR_SPACE_CHARGE_HELPER 1

#include "gsl/gsl_interp.h"
#include "gsl/gsl_poly.h"

#include "TGraph.h"
#include "TFile.h"

#include <string>
#include <vector>
#include <array>
#include <algorithm> 

#include "Objects/CartesianVector.h"

namespace lar_content
{
//------------------------------------------------------------------------------------------------------------------------------------------
//Adapted from LArSoft module

    class SortedTGraphFromFile 
    {
        TGraph* graph = nullptr;
      
        public:

        SortedTGraphFromFile(TFile& file, char const* graphName)
        : graph(static_cast<TGraph*>(file.Get(graphName)))
        {
          graph->Sort();
        }
      
      ~SortedTGraphFromFile() { delete graph; }
      
      TGraph const& operator* () const { return *graph; }
    }; 

//------------------------------------------------------------------------------------------------------------------------------------------
//Adapted from LArSoft module

    class Interpolator {
    
    std::vector<double> xa, ya; 
    gsl_interp* fInterp = nullptr; 
    
    void swap(Interpolator& with)
      {
        std::swap(xa, with.xa);
        std::swap(ya, with.ya);
        std::swap(fInterp, with.fInterp);
      }
    
      public:
    
    Interpolator() = default;
    
    Interpolator(TGraph const& graph)
      : xa(graph.GetX(), graph.GetX() + graph.GetN())
      , ya(graph.GetY(), graph.GetY() + graph.GetN())
      , fInterp(gsl_interp_alloc(gsl_interp_linear, xa.size()))
      {
        gsl_interp_init(fInterp, xa.data(), ya.data(), xa.size());
      }
    
    ~Interpolator()
      { if (fInterp) { gsl_interp_free(fInterp); fInterp = nullptr; } }
    
    Interpolator(Interpolator const&) = delete;
    Interpolator(Interpolator&& from) { swap(from); }
    Interpolator& operator= (Interpolator const&) = delete;
    Interpolator& operator= (Interpolator&& from)
      { swap(from); return *this; }
    
    double Eval(double x) const
      {
        if (x <= xa.front()) return ya.front();
        if (x >= xa.back()) return ya.back();
        return gsl_interp_eval(fInterp, xa.data(), ya.data(), x, nullptr);
      }
  }; 
  
//------------------------------------------------------------------------------------------------------------------------------------------
  
  template <unsigned int N>
  struct PolynomialBase {
    static constexpr unsigned int Degree = N; 
    
    static constexpr unsigned int NParams = Degree + 1;
    
    using Params_t = std::array<double, NParams>; 
    
    static double Eval(double x, double const* params);
    
    static double Eval(double x, Params_t const& params)
      { return Eval(x, params.data()); }
    
  }; 

//------------------------------------------------------------------------------------------------------------------------------------------
  
  template <unsigned int N>
  class Polynomial: public PolynomialBase<N> {
    using PolyImpl_t = PolynomialBase<N>;
      public:
    static constexpr unsigned int Degree = PolyImpl_t::Degree;
    static constexpr unsigned int NParams = PolyImpl_t::NParams;
    
    using Params_t = typename PolyImpl_t::Params_t;
    
    Polynomial() = default;
    
    void SetParameters(double const* params)
      { std::copy_n(params, NParams, fParams.begin()); }
    
    using PolynomialBase<N>::Eval;
    
    double Eval(double x) const
      { return PolyImpl_t::Eval(x, fParams.data()); }
    
      protected:
    Params_t fParams; 
    
  }; 
    
//------------------------------------------------------------------------------------------------------------------------------------------

template <unsigned int N>
double PolynomialBase<N>::Eval
  (double x, double const* params)
{
  return gsl_poly_eval(params, NParams, x);
}

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArSpaceChargeHelper class
 */
class LArSpaceChargeHelper
{
public:
    /**
     *  @brief  Get the new, spacecharge-adjusted position 
     * 
     *  @param  positionVector the uncorrected position, by reference 
     * 
     *  @return the spacecharge-corrected position 
     */
    static void Configure(std::string fileName);

    static pandora::CartesianVector GetSpaceChargeCorrectedPosition(pandora::CartesianVector &positionVector);

    static pandora::CartesianVector GetPositionOffset(pandora::CartesianVector &positionVector);

    static float GetPositionOffsetX(pandora::CartesianVector &positionVector);

    static float GetPositionOffsetY(pandora::CartesianVector &positionVector);

    static float GetPositionOffsetZ(pandora::CartesianVector &positionVector);

    static pandora::CartesianVector TransformPosition(pandora::CartesianVector &positionVector);

    static float TransformX(float xPosition);

    static float TransformY(float yPosition);

    static float TransformZ(float zPosition);

    static float IsInsideDetectorBoundaries(pandora::CartesianVector &positionVector);

    static Interpolator MakeInterpolator(TFile& file, char const* graphName); 

    typedef Polynomial<6> f1_x_poly_t;
    typedef Polynomial<6> f2_x_poly_t;
    typedef Polynomial<6> f3_x_poly_t;
    typedef Polynomial<6> f4_x_poly_t;
    typedef Polynomial<6> f5_x_poly_t;
    typedef PolynomialBase<4> fFinal_x_poly_t;

    typedef Polynomial<5> f1_y_poly_t;
    typedef Polynomial<5> f2_y_poly_t;
    typedef Polynomial<5> f3_y_poly_t;
    typedef Polynomial<5> f4_y_poly_t;
    typedef Polynomial<5> f5_y_poly_t;
    typedef Polynomial<5> f6_y_poly_t;
    typedef PolynomialBase<5> fFinal_y_poly_t;

    typedef Polynomial<4> f1_z_poly_t;
    typedef Polynomial<4> f2_z_poly_t;
    typedef Polynomial<4> f3_z_poly_t;
    typedef Polynomial<4> f4_z_poly_t;
    typedef PolynomialBase<3> fFinal_z_poly_t;
};

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // #ifndef LAR_SPACE_CHARGE_HELPER
