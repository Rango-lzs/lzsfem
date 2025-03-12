/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "FEM_EXPORT.h"
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
public: FactorizationError() : FEException("Fatal error in factorization of stiffness matrix. Aborting run.") {}
};

class FEM_EXPORT NegativeJacobianDetected : public FEException {
public: NegativeJacobianDetected() : FEException("Negative jacobian was detected.") {}
};

class FEM_EXPORT ConcentrationChangeDetected : public FEException {
public:
    ConcentrationChangeDetected() : FEException("Concentration change detected") {}
    ConcentrationChangeDetected(const FENodalDofInfo& ndi);
};

