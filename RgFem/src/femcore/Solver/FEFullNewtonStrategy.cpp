
#include "FEFullNewtonStrategy.h"
#include "femcore/NewtonSolver/SolverBase.h"
#include "femcore/FEException.h"
#include "femcore/NewtonSolver/NewtonSolver.h"

//-----------------------------------------------------------------------------
FEFullNewtonStrategy::FEFullNewtonStrategy() : FENewtonStrategy()
{
	m_plinsolve = nullptr;
	m_maxups = 0;
}

//-----------------------------------------------------------------------------
bool FEFullNewtonStrategy::Init()
{
	if (m_pns == nullptr) return false;
	m_plinsolve = m_pns->GetLinearSolver();
	return true;
}

//-----------------------------------------------------------------------------
bool FEFullNewtonStrategy::Update(double s, std::vector<double>& ui, std::vector<double>& R0, std::vector<double>& R1)
{
	// always return false to force a Matrix reformation
	return false;
}

//-----------------------------------------------------------------------------
void FEFullNewtonStrategy::SolveEquations(std::vector<double>& x, std::vector<double>& b)
{
	// perform a backsubstitution
	if (m_plinsolve->BackSolve(x, b) == false)
	{
		throw LinearSolverFailed();
	}
}
