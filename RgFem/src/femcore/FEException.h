#pragma once
#include "femcore/fem_export.h"
#include <string>

class FEElement;

class FEM_EXPORT FEException
{
public:
	enum {
		Error = 0,
		Warning = 1
	};

public:
	// level = 0, error
	// level = 1, warning
	FEException(const char* msg = nullptr, int level = 0);
	virtual ~FEException();

	const char* what();

	void what(const char* msg, ...);

	int level() const;

private:
	std::string	m_what;
	int	m_level;
};

class FEM_EXPORT NegativeJacobian : public FEException
{
public:
	NegativeJacobian(int iel, int ng, double vol, FEElement* pe = 0);

	int		m_iel;	// element where the jacobian was negative
	int		m_ng;	// integration point
	double	m_vol;	// volume
	FEElement*	m_pel;	// pointer to element

	static bool DoOutput();

	static void clearFlag();
	static bool IsThrown();

public:
	static bool m_boutput;	//!< set to false to suppress output of negative jacobians
	static bool m_bthrown;
};

class FEM_EXPORT ZeroDiagonal : public FEException
{
private:
	struct EQUATION
	{
		int	node;	// node
		int	dof;	// degree of node
	};

public:
	ZeroDiagonal(int node, int dof);
};

class FEM_EXPORT EnergyDiverging : public FEException {
public: EnergyDiverging() : FEException("Problem diverging uncontrollably.") {}
};

class FEM_EXPORT MaxStiffnessReformations : public FEException {
public: MaxStiffnessReformations() : FEException("Max nr of reformations reached.") {}
};

class FEM_EXPORT ZeroLinestepSize : public FEException {
public: ZeroLinestepSize() : FEException("Zero line step size.") {}
};

class FEM_EXPORT ForceConversion : public FEException {
public: ForceConversion() : FEException("User forced conversion.\nSolution might not be stable.", FEException::Warning) {}
};

class FEM_EXPORT IterationFailure : public FEException {
public: IterationFailure() : FEException("User forced iteration failure.", FEException::Warning) {}
};

class FEM_EXPORT MaxResidualError : public FEException {
public: MaxResidualError() : FEException("Maximum residual exceeded.", FEException::Warning) {}
};

struct FEM_EXPORT FENodalDofInfo;

class FEM_EXPORT NANInResidualDetected : public FEException {
public: 
	NANInResidualDetected() : FEException("NAN detected") {}
	NANInResidualDetected(const FENodalDofInfo& ndi);
};

class FEM_EXPORT NANInSolutionDetected : public FEException {
public:
	NANInSolutionDetected() : FEException("NAN detected") {}
	NANInSolutionDetected(const FENodalDofInfo& ndi);
};

class FEM_EXPORT FatalError : public FEException{
public: FatalError() : FEException("Fatal error") {}
};

class FEM_EXPORT DoRunningRestart : public FEException {
public: DoRunningRestart() : FEException("Running restart requested", FEException::Warning) {}
};

class FEM_EXPORT FEMultiScaleException : public FEException
{
public:
	FEMultiScaleException(int eid, int gpt);
};

class FEM_EXPORT LinearSolverFailed : public FEException {
public: LinearSolverFailed() : FEException("Linear solver failed to find solution. Aborting run.") {}
};

class FEM_EXPORT FactorizationError : public FEException{
public: FactorizationError() : FEException("Fatal error in factorization of stiffness Matrix. Aborting run.") {}
};

class FEM_EXPORT NegativeJacobianDetected : public FEException {
public: NegativeJacobianDetected() : FEException("Negative jacobian was detected.") {}
};

class FEM_EXPORT ConcentrationChangeDetected : public FEException {
public:
    ConcentrationChangeDetected() : FEException("Concentration change detected") {}
    ConcentrationChangeDetected(const FENodalDofInfo& ndi);
};

