/*****************************************************************//**
 * \file   fei1d_line_linear.h
 * \brief  
 * 
 * \author Leizs
 * \date   November 2023
 *********************************************************************/
#ifndef FEI1dLinear_h
#define FEI1dLinear_h

#include "feinterpol1d.h"

namespace fem 
{
/**
 * Class representing a 1d linear isoparametric interpolation.
 */
class FEM_EXPORT FEI1dLinear : public FEInterpolation1d
{
public:
    FEI1dLinear(): FEInterpolation1d(1) { }

    double giveLength() const override;

    static FloatArrayF<2> evalN(double ksi);
    std::pair<double, FloatMatrixF<1,2>> evaldNdx() const;

    void evalN(FloatArray &answer, const FloatArray &lcoords) const override;
    double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords) const override;
    void local2global(FloatArray &answer, const FloatArray &lcoords) const override;
    int  global2local(FloatArray &answer, const FloatArray &lcoords) const override;
    double giveTransformationJacobian(const FloatArray &lcoords) const override;
private:
    int cindx;
};
} // end namespace fem
#endif // FEI1dLinear_h
