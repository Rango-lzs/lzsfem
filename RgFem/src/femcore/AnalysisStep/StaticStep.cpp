#include "StaticStep.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>

StaticStep::StaticStep(const std::string& name)
    : loadMultiplier_(1.0)
    , nonlinear_(false)
    , largeDisplacement_(false)
    , useLineSearch_(false)
    , lineSearchMaxIter_(10)
    , displacementNorm_(0.0)
    , forceNorm_(0.0)
    , energyNorm_(0.0)
    , initialForceNorm_(0.0)
{
    stepName_ = name;
    currentTime_ = 0.0;
    currentIncrement_ = 0;
    outputFrequency_ = 1;
}

void StaticStep::initialize(FEModel* model)
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "Initializing Static Step: " << stepName_ << std::endl;
    std::cout << "========================================" << std::endl;
    
    if (nonlinear_) {
        std::cout << "Analysis type: Nonlinear Static" << std::endl;
        if (largeDisplacement_) {
            std::cout << "Geometric nonlinearity: Enabled" << std::endl;
        }
    } else {
        std::cout << "Analysis type: Linear Static" << std::endl;
    }
    
    currentTime_ = 0.0;
    currentIncrement_ = 0;
    
    // Initialize solution vectors
    int nDOF = model->getNumberOfDOFs();
    displacement_.assign(nDOF, 0.0);
    internalForce_.assign(nDOF, 0.0);
    externalForce_.assign(nDOF, 0.0);
    residual_.assign(nDOF, 0.0);
    reactionForces_.assign(nDOF, 0.0);
    
    std::cout << "Number of DOFs: " << nDOF << std::endl;
    std::cout << "Initialization complete.\n" << std::endl;
}

bool StaticStep::execute(FEModel* model)
{
    std::cout << "Executing static analysis..." << std::endl;
    
    bool success = false;
    
    if (nonlinear_) {
        success = solveNonlinearSystem(model);
    } else {
        success = solveLinearSystem(model);
    }
    
    if (success) {
        calculateReactionForces(model);
        updateStressStrain(model);
    }
    
    return success;
}

void StaticStep::finalize(FEModel* model)
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "Static Step Completed: " << stepName_ << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Total increments: " << currentIncrement_ << std::endl;
    std::cout << "Final load factor: " << currentTime_ << std::endl;
    
    // Output summary statistics
    double maxDisp = 0.0;
    double maxReaction = 0.0;
    
    for (double d : displacement_) {
        maxDisp = std::max(maxDisp, std::abs(d));
    }
    for (double r : reactionForces_) {
        maxReaction = std::max(maxReaction, std::abs(r));
    }
    
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "Max displacement: " << maxDisp << std::endl;
    std::cout << "Max reaction force: " << maxReaction << std::endl;
    std::cout << std::endl;
}

bool StaticStep::solveLinearSystem(FEModel* model)
{
    std::cout << "\n--- Linear Static Analysis ---\n" << std::endl;
    
    // Step 1: Assemble global stiffness matrix
    std::cout << "Assembling global stiffness matrix..." << std::endl;
    assembleStiffnessMatrix(model, false);
    
    // Step 2: Assemble load vector
    std::cout << "Assembling load vector..." << std::endl;
    assembleLoadVector(model, loadMultiplier_);
    
    // Step 3: Apply boundary conditions
    std::cout << "Applying boundary conditions..." << std::endl;
    applyBoundaryConditions(model);
    
    // Step 4: Solve linear system K*u = F
    std::cout << "Solving linear system..." << std::endl;
    
    // TODO: Call actual linear solver
    // For now, placeholder - in real implementation would call:
    // model->getSolver()->solve(K, F, displacement_);
    
    std::cout << "Linear solution obtained." << std::endl;
    
    currentIncrement_ = 1;
    currentTime_ = loadMultiplier_;
    reportProgress(1);
    
    return true;
}

bool StaticStep::solveNonlinearSystem(FEModel* model)
{
    std::cout << "\n--- Nonlinear Static Analysis (Newton-Raphson) ---\n" << std::endl;
    
    double timeIncrement = control_.initialTimeIncrement;
    currentTime_ = 0.0;
    
    std::cout << std::setw(10) << "Increment" 
              << std::setw(12) << "Time/Load" 
              << std::setw(10) << "Iters" 
              << std::setw(15) << "Disp Norm"
              << std::setw(15) << "Force Norm" << std::endl;
    std::cout << std::string(72, '-') << std::endl;
    
    while (currentTime_ < control_.totalTime) {
        currentIncrement_++;
        double nextTime = std::min(currentTime_ + timeIncrement, control_.totalTime);
        double lambda = nextTime / control_.totalTime * loadMultiplier_;
        
        // Newton-Raphson iteration loop
        bool converged = false;
        int iteration = 0;
        
        for (iteration = 0; iteration < control_.maxIterations; ++iteration) {
            // Perform one Newton-Raphson iteration
            converged = performNewtonIteration(model, lambda);
            
            reportProgress(iteration + 1);
            
            if (converged) {
                currentTime_ = nextTime;
                
                std::cout << std::setw(10) << currentIncrement_ 
                          << std::setw(12) << std::fixed << std::setprecision(4) << currentTime_
                          << std::setw(10) << (iteration + 1)
                          << std::setw(15) << std::scientific << displacementNorm_
                          << std::setw(15) << forceNorm_ << std::endl;
                break;
            }
        }
        
        // Check if converged
        if (!converged) {
            std::cout << "\n*** WARNING: Failed to converge in " 
                      << control_.maxIterations << " iterations ***" << std::endl;
            
            // Try adaptive time stepping
            if (control_.useAdaptiveTimeStep && 
                timeIncrement > control_.minimumTimeIncrement) {
                
                timeIncrement *= control_.cutbackFactor;
                currentIncrement_--;
                
                std::cout << "Reducing time increment to " << timeIncrement 
                          << " and retrying..." << std::endl;
                
                // Restore previous solution
                // In real implementation: restore displacement_ from backup
                continue;
            } else {
                std::cerr << "\nERROR: Analysis failed - no convergence!" << std::endl;
                return false;
            }
        }
        
        // Adaptive time step increase for good convergence
        if (control_.useAdaptiveTimeStep && 
            iteration < control_.minIterationsForIncrease &&
            timeIncrement < control_.maximumTimeIncrement) {
            
            timeIncrement = std::min(timeIncrement * control_.increaseFactorGood, 
                                    control_.maximumTimeIncrement);
        }
        
        // Check if we've reached the end
        if (currentTime_ >= control_.totalTime) {
            break;
        }
    }
    
    std::cout << std::string(72, '-') << std::endl;
    std::cout << "Nonlinear analysis completed successfully." << std::endl;
    
    return true;
}

bool StaticStep::performNewtonIteration(FEModel* model, double lambda)
{
    // Step 1: Assemble tangent stiffness matrix
    assembleStiffnessMatrix(model, largeDisplacement_);
    
    // Step 2: Calculate residual: R = F_ext - F_int
    calculateResidual(model, lambda);
    
    // Step 3: Check convergence before solving
    if (currentIncrement_ > 1 && checkConvergence(0)) {
        return true;
    }
    
    // Step 4: Apply boundary conditions
    applyBoundaryConditions(model);
    
    // Step 5: Solve for displacement increment: K_t * delta_u = R
    std::vector<double> deltaU(displacement_.size(), 0.0);
    
    // TODO: Call actual solver
    // model->getSolver()->solve(K_tangent, residual_, deltaU);
    
    // Step 6: Line search (optional)
    double alpha = 1.0;
    if (useLineSearch_) {
        alpha = performLineSearch(model, deltaU);
    }
    
    // Step 7: Update solution
    for (size_t i = 0; i < displacement_.size(); ++i) {
        displacement_[i] += alpha * deltaU[i];
    }
    
    // Update displacement norm for convergence check
    displacementNorm_ = 0.0;
    for (double du : deltaU) {
        displacementNorm_ += du * du;
    }
    displacementNorm_ = std::sqrt(displacementNorm_);
    
    // Step 8: Check convergence
    return checkConvergence(1);
}

void StaticStep::assembleStiffnessMatrix(FEModel* model, bool geometric)
{
    // Assemble global stiffness matrix from element stiffness matrices
    
    // TODO: Actual implementation would:
    // 1. Loop over all elements
    // 2. Calculate element stiffness matrix (material + geometric if needed)
    // 3. Assemble into global matrix
    
    /*
    model->clearStiffnessMatrix();
    
    for (auto& element : model->getMesh()->getElements()) {
        MatrixXd Ke = element->getStiffnessMatrix();
        
        if (geometric) {
            MatrixXd Kg = element->getGeometricStiffness();
            Ke += Kg;
        }
        
        std::vector<int> dofs = element->getDOFs();
        model->assembleElementMatrix(Ke, dofs);
    }
    */
}

void StaticStep::assembleLoadVector(FEModel* model, double lambda)
{
    // Assemble global load vector
    
    // TODO: Actual implementation would:
    // 1. Apply nodal loads
    // 2. Apply distributed loads
    // 3. Apply body forces
    // 4. Multiply by load factor lambda
    
    /*
    externalForce_.assign(externalForce_.size(), 0.0);
    
    // Nodal loads
    for (auto& load : model->getNodalLoads()) {
        int dof = load.getDOF();
        externalForce_[dof] += lambda * load.getValue();
    }
    
    // Distributed loads on elements
    for (auto& element : model->getMesh()->getElements()) {
        VectorXd Fe = element->getEquivalentNodalForces();
        std::vector<int> dofs = element->getDOFs();
        
        for (size_t i = 0; i < dofs.size(); ++i) {
            externalForce_[dofs[i]] += lambda * Fe[i];
        }
    }
    */
}

void StaticStep::calculateResidual(FEModel* model, double lambda)
{
    // Calculate residual: R = F_ext - F_int
    
    assembleLoadVector(model, lambda);
    
    // Calculate internal forces from current displacement
    internalForce_.assign(internalForce_.size(), 0.0);
    
    // TODO: Actual implementation would:
    /*
    for (auto& element : model->getMesh()->getElements()) {
        VectorXd Fint = element->getInternalForces(displacement_);
        std::vector<int> dofs = element->getDOFs();
        
        for (size_t i = 0; i < dofs.size(); ++i) {
            internalForce_[dofs[i]] += Fint[i];
        }
    }
    */
    
    // Calculate residual
    forceNorm_ = 0.0;
    for (size_t i = 0; i < residual_.size(); ++i) {
        residual_[i] = externalForce_[i] - internalForce_[i];
        forceNorm_ += residual_[i] * residual_[i];
    }
    forceNorm_ = std::sqrt(forceNorm_);
    
    // Store initial force norm for relative convergence check
    if (initialForceNorm_ == 0.0) {
        initialForceNorm_ = forceNorm_;
    }
}

void StaticStep::applyBoundaryConditions(FEModel* model)
{
    // Apply displacement boundary conditions
    
    // TODO: Actual implementation would:
    /*
    for (auto& bc : model->getBoundaryConditions()) {
        int dof = bc.getDOF();
        double value = bc.getValue();
        
        // Modify stiffness matrix and load vector for prescribed DOF
        // Common methods:
        // 1. Penalty method
        // 2. Lagrange multiplier
        // 3. Direct elimination
    }
    */
}

bool StaticStep::checkConvergence(int iteration)
{
    bool converged = true;
    
    // Displacement norm convergence
    if (control_.checkDisplacementNorm) {
        double dispNormRef = 0.0;
        for (double d : displacement_) {
            dispNormRef += d * d;
        }
        dispNormRef = std::sqrt(dispNormRef);
        
        if (dispNormRef > 1e-10) {
            double relDispNorm = displacementNorm_ / dispNormRef;
            if (relDispNorm > control_.displacementTolerance) {
                converged = false;
            }
        } else {
            if (displacementNorm_ > control_.displacementTolerance) {
                converged = false;
            }
        }
    }
    
    // Force norm convergence
    if (control_.checkForceNorm) {
        double relForceNorm = initialForceNorm_ > 1e-10 ? 
                             forceNorm_ / initialForceNorm_ : forceNorm_;
        
        if (relForceNorm > control_.forceTolerance) {
            converged = false;
        }
    }
    
    // Energy norm convergence
    if (control_.checkEnergyNorm) {
        energyNorm_ = 0.0;
        for (size_t i = 0; i < residual_.size(); ++i) {
            energyNorm_ += residual_[i] * displacement_[i];
        }
        energyNorm_ = std::abs(energyNorm_);
        
        if (energyNorm_ > control_.energyTolerance) {
            converged = false;
        }
    }
    
    return converged;
}

void StaticStep::updateSolution(FEModel* model, const std::vector<double>& deltaU)
{
    for (size_t i = 0; i < displacement_.size(); ++i) {
        displacement_[i] += deltaU[i];
    }
}

double StaticStep::performLineSearch(FEModel* model, const std::vector<double>& direction)
{
    // Simple backtracking line search
    double alpha = 1.0;
    double rho = 0.5;      // Reduction factor
    double c = 0.5;        // Armijo constant
    
    // Save current displacement
    std::vector<double> u_old = displacement_;
    
    // Calculate initial residual norm
    double residual0 = forceNorm_;
    
    for (int i = 0; i < lineSearchMaxIter_; ++i) {
        // Try step with current alpha
        displacement_ = u_old;
        for (size_t j = 0; j < displacement_.size(); ++j) {
            displacement_[j] += alpha * direction[j];
        }
        
        // Recalculate residual
        calculateResidual(model, currentTime_ / control_.totalTime * loadMultiplier_);
        
        // Check Armijo condition
        if (forceNorm_ <= (1.0 - c * alpha) * residual0) {
            return alpha;
        }
        
        // Reduce alpha
        alpha *= rho;
    }
    
    // Restore original displacement if line search failed
    displacement_ = u_old;
    return 1.0;
}

void StaticStep::calculateReactionForces(FEModel* model)
{
    // Calculate reaction forces at constrained DOFs
    
    // TODO: Actual implementation would:
    /*
    reactionForces_.assign(reactionForces_.size(), 0.0);
    
    for (auto& bc : model->getBoundaryConditions()) {
        int dof = bc.getDOF();
        
        // Reaction = Internal force at constrained DOF
        reactionForces_[dof] = internalForce_[dof];
    }
    */
}

void StaticStep::updateStressStrain(FEModel* model)
{
    // Update stresses and strains in all elements
    
    // TODO: Actual implementation would:
    /*
    for (auto& element : model->getMesh()->getElements()) {
        element->updateStressStrain(displacement_);
    }
    */
}
