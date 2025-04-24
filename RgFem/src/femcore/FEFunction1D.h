#pragma once
#include "femcore/FEObjectBase.h"
//#include "femcore/MathObject.h"

//-----------------------------------------------------------------------------
class FEModel;
class DumpStream;

//-----------------------------------------------------------------------------
// Class that represents a 1D function
// Derived classes need to implement:
//   - value(double): This evaluates the function at a particular point
//   - copy()       : Create a copy of the class
//   - derive(double) : Calculate the derivative. This is optional, but allows implementation of more efficient algorithm, since default implements forward difference
//
class FEM_EXPORT FEFunction1D : public FEObjectBase
{
    DECLARE_META_CLASS(FEFunction1D, FEObjectBase);

public:
	FEFunction1D(FEModel* pfem);

	// serialization
	void Serialize(DumpStream& ar);

	// this class requires a copy member
	virtual FEFunction1D* copy() = 0;

	// evaluate the function at x
	// must be defined by derived classes
	virtual double value(double x) const = 0;

	// value of first derivative of function at x
	// can be overridden by derived classes.
	// default implementation is a forward-difference
	virtual double derive(double x) const;

    virtual double deriv2(double x) const = 0;

	// value of first definite integral of funciton from a to b
	// can be overridden by derived classes.
	// default implementation is trapezoidal rule
	virtual double integrate(double a, double b) const;
    
	virtual void Clear() {}
    
    // invert function
    virtual bool invert(const double f0, double &x);
};

//-----------------------------------------------------------------------------
// A constant function
class FEM_EXPORT FEConstFunction : public FEFunction1D
{
public:
	FEConstFunction(FEModel* fem) : FEFunction1D(fem), m_value(0.0) {}
	FEFunction1D* copy() override { return new FEConstFunction(GetFEModel(), m_value); }

	double value(double t) const override { return m_value;	}
	double derive(double t) const override { return 0.0; }
	double deriv2(double t) const override { return 0.0; }

protected:
	FEConstFunction(FEModel* fem, double val) : FEFunction1D(fem), m_value(val) {}

private:
	double	m_value;

	DECLARE_PARAM_LIST();
};

//-----------------------------------------------------------------------------
// A linear function
class FEM_EXPORT FELinearFunction : public FEFunction1D
{
public:
	FELinearFunction(FEModel* fem) : FEFunction1D(fem), m_slope(0.0), m_intercept(0.0) {}
	FELinearFunction(FEModel* fem, double m, double y0) : FEFunction1D(fem), m_slope(m), m_intercept(y0) {}
	FEFunction1D* copy() override { return new FELinearFunction(GetFEModel(), m_slope, m_intercept); }

	double value(double t) const override
	{
		return m_slope*t + m_intercept;
	}

	double derive(double t) const override
	{
		return m_slope;
	}

    double deriv2(double t) const override
    {
        return 0;
    }

private:
	double	m_slope;
	double	m_intercept;

	DECLARE_PARAM_LIST();
};

//-----------------------------------------------------------------------------
// A step function
class FEM_EXPORT FEStepFunction : public FEFunction1D
{
public:
	FEStepFunction(FEModel* fem) : FEFunction1D(fem), m_x0(0.0), m_leftVal(0.0), m_rightVal(1.0) {}
	FEStepFunction(FEModel* fem, double x0, double lv, double rv) : FEFunction1D(fem), m_x0(x0), m_leftVal(lv), m_rightVal(rv) {}
	FEFunction1D* copy() override { return new FEStepFunction(GetFEModel(), m_x0, m_leftVal, m_rightVal); }

	double value(double t) const override
	{
		return (t < m_x0 ? m_leftVal : m_rightVal);
	}

	double derive(double t) const override
	{
		return 0.0;
	}

	double deriv2(double t) const override
	{
		return 0.0;
	}

    // invert function has no unique solution
    bool invert(const double f0, double &x) override { return false; }
    
private:
	double	m_x0;
	double	m_leftVal;
	double	m_rightVal;

	DECLARE_PARAM_LIST();
};

//
////-----------------------------------------------------------------------------
////! function defined via math expression
//class FEM_EXPORT FEMathFunction : public FEFunction1D
//{
//public:
//	FEMathFunction(FEModel* fem);
//
//	bool Init() override;
//
//	void Serialize(DumpStream& ar) override;
//
//	FEFunction1D* copy() override;
//
//	double value(double t) const override;
//
//	double derive(double t) const override;
//
//    double deriv2(double t) const override;
//
//	void SetMathString(const std::string& s);
//
//private:
//	void evalParams(std::vector<double>& val, double t) const;
//
//	bool BuildMathExpressions();
//
//private:
//	std::string			m_s;
//	int					m_ix;			// index of independent variable
//	std::vector<FEParamValue>	m_var;	// list of model parameters that are used as variables in expression.
//
//	MSimpleExpression	m_exp;
//	MSimpleExpression	m_dexp;
//    MSimpleExpression   m_d2exp;
//
//	DECLARE_PARAM_LIST();
//};
