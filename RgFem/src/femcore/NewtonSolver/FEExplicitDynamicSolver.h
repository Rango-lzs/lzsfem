/*********************************************************************
 * \file   FEExplicitDynamicSolver.h
 * \brief  Explicit dynamic analysis solver
 *
 * \author Refactored
 * \date   February 2025
 *********************************************************************/

#pragma once
#include "FESolver_refactored.h"
#include "femcore/FEGlobalVector.h"

//-----------------------------------------------------------------------------
//! Explicit time integration scheme
enum ExplicitScheme
{
    CENTRAL_DIFFERENCE,  //!< Central difference method
    FORWARD_EULER,       //!< Forward Euler method
    VELOCITY_VERLET      //!< Velocity Verlet method
};

//-----------------------------------------------------------------------------
//! Solver for explicit dynamic analysis
//! 
//! This solver handles dynamic problems using explicit time integration.
//! It is conditionally stable and requires small time steps.
//! Advantages: No matrix inversion required, suitable for wave propagation.
//! 
//! The central difference method solves:
//! M * a(n) = F(n) - Fint(n)
//! v(n+1/2) = v(n-1/2) + a(n) * dt
//! u(n+1) = u(n) + v(n+1/2) * dt
class FEM_EXPORT FEExplicitDynamicSolver : public FESolver
{
    DECLARE_META_CLASS(FEExplicitDynamicSolver, FESolver);

public:
    //! Constructor
    FEExplicitDynamicSolver();

    //! Destructor
    virtual ~FEExplicitDynamicSolver();

public:
    //! Serialization
    void Serialize(DumpStream& ar) override;

    //! Initialization
    bool Init() override;

    //! Initialize step
    bool InitStep(double time) override;

    //! Initialize equations
    bool InitEquations() override;

    //! Solve a time step
    bool SolveStep() override;

    //! Update model state
    void Update(std::vector<double>& u) override;

public:
    //! Get stiffness matrix (not used in explicit)
    FEGlobalMatrix* GetStiffnessMatrix() override;

    //! Get load vector
    std::vector<double> GetLoadVector() override;

    //! Get linear solver (not used in explicit)
    LinearSolver* GetLinearSolver() override;

protected:
    //! Calculate mass matrix (lumped or consistent)
    virtual bool CalculateMassMatrix();

    //! Calculate internal forces
    virtual void InternalForces(FEGlobalVector& R);

    //! Calculate external forces
    virtual void ExternalForces(FEGlobalVector& R);

    //! Calculate contact forces
    virtual void ContactForces(FEGlobalVector& R);

    //! Calculate accelerations
    virtual void CalculateAccelerations();

    //! Update velocities
    virtual void UpdateVelocities();

    //! Update displacements
    virtual void UpdateDisplacements();

    //! Check stability (CFL condition)
    virtual bool CheckStability();

    //! Calculate critical time step
    virtual double CalculateCriticalTimeStep();

protected:
    //! Initialize velocity for central difference
    void InitializeVelocity();

    //! Perform central difference integration
    void CentralDifferenceStep();

    //! Perform forward Euler integration
    void ForwardEulerStep();

    //! Perform velocity Verlet integration
    void VelocityVerletStep();

public:
    // Solver parameters
    ExplicitScheme m_scheme;     //!< time integration scheme
    bool m_lumpedMass;           //!< use lumped mass matrix
    bool m_autoTimeStep;         //!< automatic time step control
    double m_timeStepScale;      //!< time step scale factor (safety factor)
    double m_maxTimeStep;        //!< maximum allowed time step
    bool m_checkStability;       //!< check stability condition
    
    // Damping
    bool m_useDamping;           //!< enable numerical damping
    double m_dampingCoeff;       //!< damping coefficient

protected:
    FEDofList m_dofU;  //!< displacement DOFs
    FEDofList m_dofV;  //!< velocity DOFs
    FEDofList m_dofA;  //!< acceleration DOFs

    // State vectors
    std::vector<double> m_Un;    //!< displacement at time n
    std::vector<double> m_Vn;    //!< velocity at time n
    std::vector<double> m_An;    //!< acceleration at time n
    std::vector<double> m_Vnhalf; //!< velocity at time n-1/2 (for central difference)

    // Force vectors
    std::vector<double> m_Fint;  //!< internal forces
    std::vector<double> m_Fext;  //!< external forces
    std::vector<double> m_Fcon;  //!< contact forces

    // Mass matrix
    std::vector<double> m_Mi;    //!< lumped mass matrix (diagonal)
    FEGlobalMatrix* m_pM;        //!< consistent mass matrix (if used)

    // Time step
    double m_dtCritical;         //!< critical time step
    bool m_initialized;          //!< initialization flag

    DECLARE_PARAM_LIST();
};
