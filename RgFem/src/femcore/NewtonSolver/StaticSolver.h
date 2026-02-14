/*********************************************************************
 * \file   FEStaticSolver.h
 * \brief  Static analysis solver
 *
 * \author Refactored
 * \date   February 2025
 *********************************************************************/

#pragma once
#include "femcore/NewtonSolver/NewtonSolver.h"
#include "femcore/FEGlobalVector.h"

//-----------------------------------------------------------------------------
//! Solver for static (quasi-static) analysis
//! 
//! This solver handles equilibrium problems without inertia effects.
//! It solves: K * u = F
class FEM_EXPORT FEStaticSolver : public FENewtonSolver
{
    DECLARE_META_CLASS(FEStaticSolver, FENewtonSolver);

public:
    //! Constructor
    FEStaticSolver();

    //! Destructor
    virtual ~FEStaticSolver();

public:
    //! Serialization
    void Serialize(DumpStream& ar) override;

    //! Initialization
    bool Init() override;

    //! Initialize step
    bool InitStep(double time) override;

    //! Initialize equations
    bool InitEquations() override;

public:
    //! Prepare for Newton iteration
    void PrepStep() override;

    //! Calculate stiffness matrix
    bool StiffnessMatrix() override;

    //! Calculate residual vector
    bool Residual(std::vector<double>& R) override;

    //! Update model state
    void Update(std::vector<double>& ui) override;

protected:
    //! Calculate internal forces
    void InternalForces(FEGlobalVector& R);

    //! Calculate external forces
    void ExternalForces(FEGlobalVector& R);

    //! Calculate contact forces
    void ContactForces(FEGlobalVector& R);

    //! Calculate nonlinear constraint forces
    void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);

    //! Calculate contact stiffness
    void ContactStiffness(FELinearSystem& LS);

    //! Calculate nonlinear constraint stiffness
    void NonLinearConstraintStiffness(FELinearSystem& LS, const FETimeInfo& tp);

public:
    // Solver parameters
    double m_Dtol;  //!< displacement tolerance

protected:
    FEDofList m_dofU;  //!< displacement DOFs

    // Force vectors
    std::vector<double> m_Fr;    //!< reaction forces
    std::vector<double> m_Fint;  //!< internal forces
    std::vector<double> m_Fext;  //!< external forces

    DECLARE_PARAM_LIST();
};
