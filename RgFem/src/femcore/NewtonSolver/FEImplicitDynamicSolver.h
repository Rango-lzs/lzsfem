/*********************************************************************
 * \file   FEImplicitDynamicSolver.h
 * \brief  Implicit dynamic analysis solver
 *
 * \author Refactored
 * \date   February 2025
 *********************************************************************/

#pragma once
#include "FENewtonSolver_refactored.h"
#include "femcore/FEGlobalVector.h"

//-----------------------------------------------------------------------------
//! Time integration scheme
enum TimeIntegrationScheme
{
    NEWMARK,            //!< Newmark-beta method
    GENERALIZED_ALPHA,  //!< Generalized-alpha method
    HHT_ALPHA          //!< Hilber-Hughes-Taylor alpha method
};

//-----------------------------------------------------------------------------
//! Solver for implicit dynamic analysis
//! 
//! This solver handles dynamic problems with inertia effects using
//! implicit time integration (Newmark, Generalized-alpha, HHT).
//! It solves: M * a + C * v + K * u = F
class FEM_EXPORT FEImplicitDynamicSolver : public FENewtonSolver
{
    DECLARE_META_CLASS(FEImplicitDynamicSolver, FENewtonSolver);

public:
    //! Constructor
    FEImplicitDynamicSolver();

    //! Destructor
    virtual ~FEImplicitDynamicSolver();

public:
    //! Serialization
    void Serialize(DumpStream& ar) override;

    //! Initialization
    bool Init() override;

    //! Initialize step
    bool InitStep(double time) override;

    //! Initialize equations
    bool InitEquations() override;

    //! Generate solver warnings
    void SolverWarnings();

public:
    //! Prepare for Newton iteration
    void PrepStep() override;

    //! Calculate effective stiffness matrix
    bool StiffnessMatrix() override;

    //! Calculate residual vector (dynamic equilibrium)
    bool Residual(std::vector<double>& R) override;

    //! Update kinematics (positions, velocities, accelerations)
    void Update(std::vector<double>& ui) override;

protected:
    //! Update kinematics
    virtual void UpdateKinematics(std::vector<double>& ui);

    //! Update DOF increments
    virtual void UpdateIncrements(std::vector<double>& Ui, std::vector<double>& ui);

protected:
    //! Calculate internal forces
    void InternalForces(FEGlobalVector& R);

    //! Calculate external forces
    void ExternalForces(FEGlobalVector& R);

    //! Calculate contact forces
    void ContactForces(FEGlobalVector& R);

    //! Calculate nonlinear constraint forces
    void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);

    //! Calculate inertial forces
    void InertialForces(FEGlobalVector& R);

    //! Calculate damping forces (if applicable)
    void DampingForces(FEGlobalVector& R);

protected:
    //! Calculate contact stiffness
    void ContactStiffness(FELinearSystem& LS);

    //! Calculate nonlinear constraint stiffness
    void NonLinearConstraintStiffness(FELinearSystem& LS, const FETimeInfo& tp);

    //! Calculate mass matrix contribution
    void MassMatrix(FELinearSystem& LS);

    //! Calculate damping matrix contribution
    void DampingMatrix(FELinearSystem& LS);

protected:
    //! Initialize time integration scheme
    void InitializeTimeIntegration();

    //! Predict values for new time step
    void PredictValues();

public:
    // Time integration parameters
    TimeIntegrationScheme m_timeScheme;  //!< time integration scheme

    // Newmark parameters
    double m_beta;   //!< Newmark beta (displacement integration)
    double m_gamma;  //!< Newmark gamma (velocity integration)

    // Generalized-alpha parameters
    double m_rhoi;    //!< spectral radius
    double m_alphaf;  //!< alpha for Y={v,e}
    double m_alpham;  //!< alpha for Ydot={∂v/∂t,∂e/∂t}

    // Convergence tolerance
    double m_Dtol;  //!< displacement tolerance

    // Damping parameters
    double m_dampingRatio;  //!< damping ratio (for Rayleigh damping)
    bool m_useDamping;      //!< enable damping

protected:
    FEDofList m_dofU;  //!< displacement DOFs
    FEDofList m_dofV;  //!< velocity DOFs
    FEDofList m_dofA;  //!< acceleration DOFs (if needed)

    // Force vectors
    std::vector<double> m_Fr;    //!< reaction forces
    std::vector<double> m_Fint;  //!< internal forces
    std::vector<double> m_Fext;  //!< external forces
    std::vector<double> m_Fin;   //!< inertial forces
    std::vector<double> m_Fd;    //!< damping forces

    // Previous converged values
    std::vector<double> m_Uip;  //!< previous displacement increment

    DECLARE_PARAM_LIST();
};
