#include "LinearSolver.h"

DEFINE_META_CLASS(LinearSolver, FEObjectBase, "");

//-----------------------------------------------------------------------------
LinearSolver::LinearSolver() : FEObjectBase()
{
	ResetStats();
}

//-----------------------------------------------------------------------------
LinearSolver::~LinearSolver()
{ 
	Destroy();
}

//-----------------------------------------------------------------------------
// returns whether this is an iterative solver or not
bool LinearSolver::IsIterative() const
{
	return false;
}

//-----------------------------------------------------------------------------
bool LinearSolver::PreProcess()
{ 
	return true; 
}

//-----------------------------------------------------------------------------
void LinearSolver::SetPartitions(const std::vector<int>& part)
{
	m_part = part;
}

//-----------------------------------------------------------------------------
void LinearSolver::SetPartitions(int npart0, int npart1)
{
	m_part.resize(2);
	m_part[0] = npart0;
	m_part[1] = npart1;
}

//-----------------------------------------------------------------------------
// nr of partitions
int LinearSolver::Partitions() const
{
	return (int)m_part.size();
}

//-----------------------------------------------------------------------------
// get the size of a partition
int LinearSolver::GetPartitionSize(int part) const
{
	return m_part[part];
}

//-----------------------------------------------------------------------------
const LinearSolverStats& LinearSolver::GetStats() const
{
	return	m_stats;
}

//-----------------------------------------------------------------------------
void LinearSolver::ResetStats()
{
	m_stats.backsolves = 0;
	m_stats.iterations = 0;
}

//-----------------------------------------------------------------------------
void LinearSolver::UpdateStats(int iterations)
{
	m_stats.backsolves++;
	m_stats.iterations += iterations;
}

//-----------------------------------------------------------------------------
void LinearSolver::Destroy()
{

}

//-----------------------------------------------------------------------------
//! helper function for when this solver is used as a preconditioner
bool LinearSolver::mult_vector(double* x, double* y)
{
	return BackSolve(y, x);
}

//-----------------------------------------------------------------------------
bool LinearSolver::SetSparseMatrix(SparseMatrix* pA)
{
	assert(false);
	return false;
}

//-----------------------------------------------------------------------------
//! convenience function for solving linear systems
bool LinearSolver::Solve(std::vector<double>& x, std::vector<double>& y)
{
	if (PreProcess() == false) return false;
	if (Factor() == false) return false;
	if (BackSolve(x, y) == false) return false;
	return true;
}

//-----------------------------------------------------------------------------
IterativeLinearSolver::IterativeLinearSolver() : LinearSolver() 
{

}

//! convenience function for solving linear systems with an iterative solver
bool IterativeLinearSolver::Solve(SparseMatrix& A, std::vector<double>& x, std::vector<double>& b, LinearSolver* pc)
{
	if (SetSparseMatrix(&A) == false) return false;
	SetLeftPreconditioner(pc);
	if (PreProcess() == false) return false;
	if (Factor() == false) return false;
	return BackSolve(x, b);
}

// returns whether this is an iterative solver or not
bool IterativeLinearSolver::IsIterative() const
{
	return true;
}

void IterativeLinearSolver::SetLeftPreconditioner(LinearSolver* pc) { assert(false); }
void IterativeLinearSolver::SetRightPreconditioner(LinearSolver* pc) { assert(false); }

// get the preconditioner
LinearSolver* IterativeLinearSolver::GetLeftPreconditioner() { return nullptr; }
LinearSolver* IterativeLinearSolver::GetRightPreconditioner() { return nullptr; }
