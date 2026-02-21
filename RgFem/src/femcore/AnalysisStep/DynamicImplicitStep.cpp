#include "DynamicImplicitStep.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <corecrt_math_defines.h>

DynamicImplicitStep::DynamicImplicitStep(const std::string& name)
    : beta_(0.25)        // Average acceleration method (unconditionally stable)
    , gamma_(0.5)
    , alphaM_(0.0)
    , betaK_(0.0)
    , currentTimeStep_(0.0)
{
    stepName_ = name;
    currentTime_ = 0.0;
    currentIncrement_ = 0;
    outputFrequency_ = 1;
}

void DynamicImplicitStep::initialize(FEModel* model)
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "Initializing Dynamic Step: " << stepName_ << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Time integration: Newmark method" << std::endl;
    std::cout << "  Beta = " << beta_ << ", Gamma = " << gamma_ << std::endl;
    
    if (beta_ == 0.25 && gamma_ == 0.5) {
        std::cout << "  (Average acceleration - unconditionally stable)" << std::endl;
    } else if (beta_ == 1.0/6.0 && gamma_ == 0.5) {
        std::cout << "  (Linear acceleration)" << std::endl;
    }
    
    if (alphaM_ != 0.0 || betaK_ != 0.0) {
        std::cout << "Rayleigh damping: C = " << alphaM_ << "*M + " 
                  << betaK_ << "*K" << std::endl;
    }
    
    currentTime_ = 0.0;
    currentIncrement_ = 0;
    currentTimeStep_ = control_.initialTimeIncrement;
    
    // Initialize solution vectors
    int nDOF = model->GetDofSchema().GetDofsPerNode();
    
    // Current state
    if (displacement_.empty()) {
        displacement_.assign(nDOF, 0.0);
    }
    if (velocity_.empty()) {
        velocity_.assign(nDOF, 0.0);
    }
    if (acceleration_.empty()) {
        acceleration_.assign(nDOF, 0.0);
    }
    
    // Previous state
    displacement_old_ = displacement_;
    velocity_old_ = velocity_;
    acceleration_old_ = acceleration_;
    
    // Force vectors
    externalForce_.assign(nDOF, 0.0);
    effectiveForce_.assign(nDOF, 0.0);
    
    // Calculate Newmark integration coefficients
    double dt = currentTimeStep_;
    a0_ = 1.0 / (beta_ * dt * dt);
    a1_ = gamma_ / (beta_ * dt);
    a2_ = 1.0 / (beta_ * dt);
    a3_ = 1.0 / (2.0 * beta_) - 1.0;
    a4_ = gamma_ / beta_ - 1.0;
    a5_ = dt / 2.0 * (gamma_ / beta_ - 2.0);
    a6_ = dt * (1.0 - gamma_);
    a7_ = gamma_ * dt;
    
    std::cout << "Number of DOFs: " << nDOF << std::endl;
    std::cout << "Time step: " << dt << " s" << std::endl;
    std::cout << "Total time: " << control_.totalTime << " s" << std::endl;
    std::cout << "Initialization complete.\n" << std::endl;
}

bool DynamicImplicitStep::execute(FEModel* model)
{
    std::cout << "Executing dynamic implicit analysis..." << std::endl;
    
    bool success = performTimeIntegration(model);
    
    return success;
}

void DynamicImplicitStep::finalize(FEModel* model)
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "Dynamic Step Completed: " << stepName_ << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Total time steps: " << currentIncrement_ << std::endl;
    std::cout << "Final time: " << currentTime_ << " s" << std::endl;
    
    // Output summary statistics
    double maxDisp = 0.0;
    double maxVel = 0.0;
    double maxAcc = 0.0;
    
    for (size_t i = 0; i < displacement_.size(); ++i) {
        maxDisp = std::max(maxDisp, std::abs(displacement_[i]));
        maxVel = std::max(maxVel, std::abs(velocity_[i]));
        maxAcc = std::max(maxAcc, std::abs(acceleration_[i]));
    }
    
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "Max displacement: " << maxDisp << " m" << std::endl;
    std::cout << "Max velocity: " << maxVel << " m/s" << std::endl;
    std::cout << "Max acceleration: " << maxAcc << " m/s^2" << std::endl;
    std::cout << std::endl;
}

bool DynamicImplicitStep::performTimeIntegration(FEModel* model)
{
    std::cout << "\n--- Time Integration (Newmark Method) ---\n" << std::endl;
    
    double dt = currentTimeStep_;
    int totalSteps = static_cast<int>(control_.totalTime / dt);
    
    std::cout << std::setw(8) << "Step" 
              << std::setw(12) << "Time (s)"
              << std::setw(15) << "Max Disp"
              << std::setw(15) << "Max Vel"
              << std::setw(15) << "Max Acc" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    while (currentTime_ < control_.totalTime) {
        currentIncrement_++;
        
        // Store previous solution
        displacement_old_ = displacement_;
        velocity_old_ = velocity_;
        acceleration_old_ = acceleration_;
        
        // Advance time
        currentTime_ += dt;
        if (currentTime_ > control_.totalTime) {
            currentTime_ = control_.totalTime;
        }
        
        // Step 1: Assemble effective stiffness matrix
        // K_eff = K + a0*M + a1*C
        assembleEffectiveStiffness(model);
        
        // Step 2: Assemble effective load vector
        // F_eff = F(t+dt) + M*(a0*u_n + a2*v_n + a3*a_n) + C*(a1*u_n + a4*v_n + a5*a_n)
        assembleEffectiveLoad(model, currentTime_);
        
        // Step 3: Apply boundary conditions
        applyBoundaryConditions(model);
        
        // Step 4: Solve for displacement at t + dt
        // K_eff * u_{n+1} = F_eff
        std::vector<double> deltaU(displacement_.size(), 0.0);
        
        // TODO: Call actual solver
        // model->getSolver()->solve(K_effective, effectiveForce_, displacement_);
        
        // For now, use a simple explicit update (for demonstration)
        // In real implementation, this should be an implicit solve
        for (size_t i = 0; i < displacement_.size(); ++i) {
            deltaU[i] = displacement_[i] - displacement_old_[i];
        }
        
        // Step 5: Update velocity and acceleration
        updateVelocityAcceleration(deltaU);
        
        // Report progress
        reportProgress();
        
        // Output results at specified frequency
        if (currentIncrement_ % outputFrequency_ == 0) {
            double maxDisp = 0.0, maxVel = 0.0, maxAcc = 0.0;
            
            for (size_t i = 0; i < displacement_.size(); ++i) {
                maxDisp = std::max(maxDisp, std::abs(displacement_[i]));
                maxVel = std::max(maxVel, std::abs(velocity_[i]));
                maxAcc = std::max(maxAcc, std::abs(acceleration_[i]));
            }
            
            std::cout << std::setw(8) << currentIncrement_
                      << std::setw(12) << std::fixed << std::setprecision(4) << currentTime_
                      << std::setw(15) << std::scientific << maxDisp
                      << std::setw(15) << maxVel
                      << std::setw(15) << maxAcc << std::endl;
        }
    }
    
    std::cout << std::string(70, '-') << std::endl;
    std::cout << "Time integration completed successfully." << std::endl;
    
    return true;
}

void DynamicImplicitStep::assembleEffectiveStiffness(FEModel* model)
{
    // Assemble: K_eff = K + a0*M + a1*C
    //         = K + a0*M + a1*(alphaM*M + betaK*K)
    //         = (1 + a1*betaK)*K + (a0 + a1*alphaM)*M
    
    // TODO: Actual implementation
    /*
    // Assemble stiffness matrix
    assembleStiffnessMatrix(model);
    
    // Assemble mass matrix
    assembleMassMatrix(model);
    
    // Compute effective stiffness
    K_effective = (1.0 + a1_ * betaK_) * K + (a0_ + a1_ * alphaM_) * M;
    */
}

void DynamicImplicitStep::assembleEffectiveLoad(FEModel* model, double time)
{
    // Assemble: F_eff = F(t+dt) + M*(a0*u_n + a2*v_n + a3*a_n) 
    //                           + C*(a1*u_n + a4*v_n + a5*a_n)
    
    // Get external load at time t+dt
    getExternalLoad(model, time, externalForce_);
    
    effectiveForce_ = externalForce_;
    
    // Add mass contributions
    for (size_t i = 0; i < effectiveForce_.size(); ++i) {
        // M * (a0*u_n + a2*v_n + a3*a_n)
        double mass_contrib = a0_ * displacement_old_[i] 
                            + a2_ * velocity_old_[i]
                            + a3_ * acceleration_old_[i];
        
        // TODO: Multiply by actual mass matrix
        // effectiveForce_[i] += M[i] * mass_contrib;
    }
    
    // Add damping contributions
    for (size_t i = 0; i < effectiveForce_.size(); ++i) {
        // C * (a1*u_n + a4*v_n + a5*a_n)
        double damp_contrib = a1_ * displacement_old_[i]
                            + a4_ * velocity_old_[i]
                            + a5_ * acceleration_old_[i];
        
        // TODO: Multiply by actual damping matrix
        // effectiveForce_[i] += C[i] * damp_contrib;
    }
}

void DynamicImplicitStep::updateVelocityAcceleration(const std::vector<double>& deltaU)
{
    // Newmark formulas:
    // a_{n+1} = a0*(u_{n+1} - u_n) - a2*v_n - a3*a_n
    // v_{n+1} = v_n + a6*a_n + a7*a_{n+1}
    
    for (size_t i = 0; i < displacement_.size(); ++i) {
        // Update acceleration
        acceleration_[i] = a0_ * deltaU[i]
                         - a2_ * velocity_old_[i]
                         - a3_ * acceleration_old_[i];
        
        // Update velocity
        velocity_[i] = velocity_old_[i]
                     + a6_ * acceleration_old_[i]
                     + a7_ * acceleration_[i];
    }
}

void DynamicImplicitStep::assembleMassMatrix(FEModel* model)
{
    // Assemble global mass matrix from element mass matrices
    
    // TODO: Actual implementation
    /*
    model->clearMassMatrix();
    
    for (auto& element : model->getMesh()->getElements()) {
        MatrixXd Me = element->getMassMatrix();
        std::vector<int> dofs = element->getDOFs();
        model->assembleElementMatrix(Me, dofs, model->getMassMatrix());
    }
    */
}

void DynamicImplicitStep::assembleDampingMatrix(FEModel* model)
{
    // Rayleigh damping: C = alphaM * M + betaK * K
    
    // TODO: Actual implementation
    /*
    assembleMassMatrix(model);
    assembleStiffnessMatrix(model);
    
    C = alphaM_ * M + betaK_ * K;
    */
}

void DynamicImplicitStep::assembleStiffnessMatrix(FEModel* model)
{
    // Assemble global stiffness matrix
    
    // TODO: Actual implementation
    /*
    model->clearStiffnessMatrix();
    
    for (auto& element : model->getMesh()->getElements()) {
        MatrixXd Ke = element->getStiffnessMatrix();
        std::vector<int> dofs = element->getDOFs();
        model->assembleElementMatrix(Ke, dofs, model->getStiffnessMatrix());
    }
    */
}

void DynamicImplicitStep::getExternalLoad(FEModel* model, double time, 
                                          std::vector<double>& load)
{
    // Get external load at given time
    load.assign(load.size(), 0.0);
    
    // TODO: Actual implementation
    /*
    // Apply time-dependent nodal loads
    for (auto& load : model->getNodalLoads()) {
        int dof = load.getDOF();
        double value = load.getValueAtTime(time);
        load[dof] += value;
    }
    
    // Apply time-dependent distributed loads
    for (auto& element : model->getMesh()->getElements()) {
        VectorXd Fe = element->getEquivalentNodalForcesAtTime(time);
        std::vector<int> dofs = element->getDOFs();
        
        for (size_t i = 0; i < dofs.size(); ++i) {
            load[dofs[i]] += Fe[i];
        }
    }
    */
}

void DynamicImplicitStep::applyBoundaryConditions(FEModel* model)
{
    // Apply displacement boundary conditions
    
    // TODO: Actual implementation
    /*
    for (auto& bc : model->getBoundaryConditions()) {
        int dof = bc.getDOF();
        double value = bc.getValue();
        
        // Modify effective stiffness and load for prescribed DOF
    }
    */
}

void DynamicImplicitStep::setNewmarkParameters(double beta, double gamma)
{
    beta_ = beta;
    gamma_ = gamma;
    
    // Recalculate integration coefficients if already initialized
    if (currentTimeStep_ > 0.0) {
        double dt = currentTimeStep_;
        a0_ = 1.0 / (beta_ * dt * dt);
        a1_ = gamma_ / (beta_ * dt);
        a2_ = 1.0 / (beta_ * dt);
        a3_ = 1.0 / (2.0 * beta_) - 1.0;
        a4_ = gamma_ / beta_ - 1.0;
        a5_ = dt / 2.0 * (gamma_ / beta_ - 2.0);
        a6_ = dt * (1.0 - gamma_);
        a7_ = gamma_ * dt;
    }
}

void DynamicImplicitStep::setRayleighDamping(double alphaM, double betaK)
{
    alphaM_ = alphaM;
    betaK_ = betaK;
}

void DynamicImplicitStep::setDampingFromFrequencies(double freq1, double freq2, 
                                                     double zeta1, double zeta2)
{
    // Convert frequencies to angular frequencies
    double omega1 = 2.0 * M_PI * freq1;
    double omega2 = 2.0 * M_PI * freq2;
    
    // Solve for Rayleigh damping coefficients
    // zeta_i = (alphaM / (2*omega_i)) + (betaK * omega_i / 2)
    
    double denom = omega1 * omega1 - omega2 * omega2;
    if (std::abs(denom) < 1e-10) {
        std::cerr << "Warning: Frequencies too close, using default damping" << std::endl;
        return;
    }
    
    alphaM_ = 2.0 * (zeta1 * omega1 * omega2 * omega2 - zeta2 * omega1 * omega1 * omega2) / denom;
    betaK_ = 2.0 * (zeta2 * omega2 - zeta1 * omega1) / denom;
    
    std::cout << "Rayleigh damping computed from frequencies:" << std::endl;
    std::cout << "  alphaM = " << alphaM_ << std::endl;
    std::cout << "  betaK = " << betaK_ << std::endl;
}

void DynamicImplicitStep::setInitialDisplacement(const std::vector<double>& u0)
{
    displacement_ = u0;
}

void DynamicImplicitStep::setInitialVelocity(const std::vector<double>& v0)
{
    velocity_ = v0;
}

void DynamicImplicitStep::setInitialAcceleration(const std::vector<double>& a0)
{
    acceleration_ = a0;
}
