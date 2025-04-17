#pragma once

#include "femcore/FELinearSystem.h"

class FERigidSolver;

class FEM_EXPORT FESolidLinearSystem : public FELinearSystem
{
public:
	FESolidLinearSystem(FESolver* solver, FERigidSolver* rigidSolver, FEGlobalMatrix& K, std::vector<double>& F, std::vector<double>& u, bool bsymm, double alpha, int nreq);

	// Assembly routine
	// This assembles the element stiffness matrix ke into the global matrix.
	// The contributions of prescribed degrees of freedom will be stored in m_F
	void Assemble(const FEElementMatrix& ke) override;

	// scale factor for stiffness matrix
	void StiffnessAssemblyScaleFactor(double a);

private:
	FERigidSolver*	m_rigidSolver;
	double			m_alpha;
	int				m_nreq;

	double	m_stiffnessScale;
};
