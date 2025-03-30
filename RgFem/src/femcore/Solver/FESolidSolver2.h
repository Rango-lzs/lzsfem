#pragma once

#include "femcore/Solver/FENewtonSolver.h"
#include "femcore/FETimeInfo.h"
#include "femcore/FEGlobalVector.h"
#include "femcore/Solver/FERigidSolver.h"
#include "femcore/FEDofList.h"

//-----------------------------------------------------------------------------
//! The FESolidSolver2 class solves large deformation solid mechanics problems
//! It can deal with quasi-static and dynamic problems
//! 
class FEM_EXPORT FESolidSolver2 : public FENewtonSolver
{
	enum ARC_LENGTH_METHOD {
		NONE,
		CRISFIELD,
	};

public:
	//! constructor
	FESolidSolver2(FEModel* pfem);

	//! destructor
	virtual ~FESolidSolver2();

	//! serialize data to/from dump file
	void Serialize(DumpStream& ar) override;

	//! Initializes data structures
	bool Init() override;

	//! initialize the step
	bool InitStep(double time) override;

	//! Initialize linear equation system
	bool InitEquations() override;

    //! Generate warnings if needed
    void SolverWarnings();

	//! Return the rigid solver
	FERigidSolver* GetRigidSolver();

public:
	//{ --- evaluation and update ---
		//! Perform an update
		void Update(std::vector<double>& ui) override;

		//! perform an updated where ui also contains displacement increments of prescribed displacements
		//! NOTE: This is a temporary hack that is only by the JFNKMatrix
		void Update2(const std::vector<double>& ui) override;

		//! update nodal positions, velocities, accelerations, etc.
		virtual void UpdateKinematics(std::vector<double>& ui);

		//! Update EAS
		void UpdateEAS(std::vector<double>& ui);
		void UpdateIncrementsEAS(std::vector<double>& ui, const bool binc);

		//! update DOF increments
		virtual void UpdateIncrements(std::vector<double>& Ui, std::vector<double>& ui, bool emap);
	//}

	//{ --- Solution functions ---

		//! prepares the data for the first QN iteration
		void PrepStep() override;

		//! Performs a Newton-Raphson iteration
		bool Quasin() override;

		//! Apply arc-length
		void DoArcLength();
	//}

	//{ --- Stiffness matrix routines ---

		//! calculates the global stiffness matrix
		virtual bool StiffnessMatrix() override;

		//! contact stiffness
		void ContactStiffness(FELinearSystem& LS);

		//! calculates stiffness contributon of nonlinear constraints
		void NonLinearConstraintStiffness(FELinearSystem& LS, const FETimeInfo& tp);
	//}

	//{ --- Residual routines ---

		//! Calculate the contact forces
		void ContactForces(FEGlobalVector& R);

		//! Calculates residual
		virtual bool Residual(std::vector<double>& R) override;

		//! Calculate nonlinear constraint forces
		void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);

		//! Internal forces
		void InternalForces(FEGlobalVector& R);

		//! external forces
		void ExternalForces(FEGlobalVector& R);
	//}

public:
	// convergence tolerances
	double	m_Dtol;			//!< displacement tolerance

	bool	m_logSolve;		//!< flag to use Aggarwal's log method

	// equation numbers
	int		m_nreq;			//!< start of rigid body equations

public:
	std::vector<double> m_Fr;	//!< nodal reaction forces
	std::vector<double> m_Fint;	//!< internal load std::vector
	std::vector<double> m_Fext;	//!< external load std::vector
	std::vector<double> m_Uip;	//!< previous converged displacement increment

    // generalized alpha method (for dynamic analyses)
    double  m_rhoi;         //!< spectral radius
    double  m_alphaf;       //!< alpha step for Y={v,e}
    double  m_alpham;       //!< alpha step for Ydot={∂v/∂t,∂e/∂t}
	double	m_alpha;		//!< Newmark parameter alpha (force integration)
	double	m_beta;			//!< Newmark parameter beta (displacement integration)
	double	m_gamma;		//!< Newmark parameter gamme (velocity integration)

	// arc-length parameters
	int		m_arcLength;	//!< arc-length method flag (0 = off, 1 = Crisfield)
	double	m_al_scale;		//!< arc-length scaling parameter (i.e. psi).
	double	m_al_lam;		//!< current arc-length lambda
	double	m_al_inc;		//!< arc-length lambda increment at current timestep
	double	m_al_ds;		//!< arc-length constraint
	double	m_al_gamma;		//!< acr-length increment at current iteration

protected:
	FEDofList	m_dofU, m_dofV;
	FEDofList	m_dofSQ;
	FEDofList	m_dofRQ;
	FEDofList	m_dofSU, m_dofSV, m_dofSA;
    
protected:
    FERigidSolverNew	m_rigidSolver;

	// declare the parameter list
	DECLARE_PARAM_LIST();
};
