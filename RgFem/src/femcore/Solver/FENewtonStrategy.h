#pragma once
#include <vector>
#include "femcore/FEObjectBase.h"

//-----------------------------------------------------------------------------
class FENewtonSolver;
class SparseMatrix;
class LinearSolver;
enum MatrixType;

//-----------------------------------------------------------------------------
//! A Base class for newton-type solution strategies
class FEM_EXPORT FENewtonStrategy : public FEObjectBase
{
    DECLARE_META_CLASS(FENewtonStrategy, FEObjectBase);

public:
	FENewtonStrategy();
	virtual ~FENewtonStrategy();

	void SetNewtonSolver(FENewtonSolver* solver);

	void Serialize(DumpStream& ar) override;

	//! reset data for new run
	virtual void Reset();

public:
	//! initialize the linear system
	virtual SparseMatrix* CreateSparseMatrix(const MatrixType& mtype);

	//! Presolve update
	virtual void PreSolveUpdate() {}

	//! perform a Newton udpate (returning false will force a Matrix reformations)
	virtual bool Update(double s, std::vector<double>& ui, std::vector<double>& R0, std::vector<double>& R1) = 0;

	//! solve the equations
	virtual void SolveEquations(std::vector<double>& x, std::vector<double>& b) = 0;

	//! reform the stiffness Matrix
	virtual bool ReformStiffness();

	//! calculate the residual
	virtual bool Residual(std::vector<double>& R, bool binit);

public:
	int		m_maxups;		//!< max nr of QN iters permitted between stiffness reformations
	int		m_max_buf_size;	//!< max buffer size for update vector storage
	bool	m_cycle_buffer;	//!< recycle the buffer when updates is larger than buffer size
	double	m_cmax;			//!< maximum value for the condition number
	int		m_nups;			//!< nr of stiffness updates

protected:
	FENewtonSolver*	m_pns;
};
