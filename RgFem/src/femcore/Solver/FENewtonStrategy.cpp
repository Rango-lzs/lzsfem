#include "FENewtonStrategy.h"
#include "femcore/NewtonSolver/NewtonSolver.h"
#include "LinearSolver.h"
#include "../fecore_enum.h"
#include "basicio/DumpStream.h"

DEFINE_META_CLASS(FENewtonStrategy, FEObjectBase, "");

FENewtonStrategy::FENewtonStrategy() : FEObjectBase()
{
	m_pns = nullptr;

	m_maxups = 10;
	m_cmax = 1e5;
	m_max_buf_size = 0; // when zero, it should default to m_maxups
	m_cycle_buffer = true;

	m_nups = 0;
}

FENewtonStrategy::~FENewtonStrategy()
{
}

void FENewtonStrategy::SetNewtonSolver(FENewtonSolver* solver)
{
	m_pns = solver;
}

//! initialize the linear system
SparseMatrix* FENewtonStrategy::CreateSparseMatrix(const MatrixType& mtype)
{
	if (m_pns == 0) return 0;

	LinearSolver* plinsolve = m_pns->GetLinearSolver();

    SparseMatrix* pS = plinsolve->CreateSparseMatrix(mtype);

	return pS;
}

bool FENewtonStrategy::ReformStiffness()
{
	return m_pns->ReformStiffness();
}

//! calculate the residual
bool FENewtonStrategy::Residual(std::vector<double>& R, bool binit)
{
	TRACK_TIME(TimerID::Timer_Residual);
	return m_pns->Residual(R);
}

void FENewtonStrategy::Serialize(DumpStream& ar)
{
	FEObjectBase::Serialize(ar);
	ar & m_nups;
	ar & m_pns;
}

//! reset data for new run
void FENewtonStrategy::Reset()
{

}
