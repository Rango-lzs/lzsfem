#ifndef DYNAMICIMPLICITSTEP_H
#define DYNAMICIMPLICITSTEP_H

#include "AnalysisStep.h"
#include <vector>

/**
 * @brief Dynamic implicit analysis step using Newmark time integration
 * 
 * Solves the equation of motion: M*a + C*v + K*u = F(t)
 * using the Newmark-beta method for time integration
 */
class DynamicImplicitStep : public AnalysisStep {
public:
    /**
     * @brief Constructor
     * @param name Step name
     */
    explicit DynamicImplicitStep(const std::string& name = "Dynamic");
    
    /**
     * @brief Destructor
     */
    ~DynamicImplicitStep() override = default;
    
    // AnalysisStep interface implementation
    bool execute(FEModel* model) override;
    void initialize(FEModel* model) override;
    void finalize(FEModel* model) override;
    AnalysisType getAnalysisType() const override { return AnalysisType::DYNAMIC_IMPLICIT; }
    
    // Newmark parameters
    void setNewmarkParameters(double beta, double gamma);
    void getNewmarkParameters(double& beta, double& gamma) const {
        beta = beta_; gamma = gamma_;
    }
    
    // Damping parameters (Rayleigh damping: C = alphaM * M + betaK * K)
    void setRayleighDamping(double alphaM, double betaK);
    void getRayleighDamping(double& alphaM, double& betaK) const {
        alphaM = alphaM_; betaK = betaK_;
    }
    
    // Set damping from frequencies
    void setDampingFromFrequencies(double freq1, double freq2, double zeta1, double zeta2);
    
    // Initial conditions
    void setInitialDisplacement(const std::vector<double>& u0);
    void setInitialVelocity(const std::vector<double>& v0);
    void setInitialAcceleration(const std::vector<double>& a0);
    
    // Get results
    const std::vector<double>& getDisplacement() const { return displacement_; }
    const std::vector<double>& getVelocity() const { return velocity_; }
    const std::vector<double>& getAcceleration() const { return acceleration_; }
    
private:
    /**
     * @brief Perform time integration using Newmark method
     * @param model FE model
     * @return true if successful
     */
    bool performTimeIntegration(FEModel* model);
    
    /**
     * @brief Assemble effective stiffness matrix
     * K_eff = K + a0*M + a1*C
     * @param model FE model
     */
    void assembleEffectiveStiffness(FEModel* model);
    
    /**
     * @brief Assemble effective load vector
     * F_eff = F(t+dt) + M*(a0*u + a2*v + a3*a) + C*(a1*u + a4*v + a5*a)
     * @param model FE model
     * @param time Current time
     */
    void assembleEffectiveLoad(FEModel* model, double time);
    
    /**
     * @brief Update velocity and acceleration using Newmark formulas
     * @param deltaU Displacement increment
     */
    void updateVelocityAcceleration(const std::vector<double>& deltaU);
    
    /**
     * @brief Assemble mass matrix
     * @param model FE model
     */
    void assembleMassMatrix(FEModel* model);
    
    /**
     * @brief Assemble damping matrix (Rayleigh damping)
     * @param model FE model
     */
    void assembleDampingMatrix(FEModel* model);
    
    /**
     * @brief Assemble stiffness matrix
     * @param model FE model
     */
    void assembleStiffnessMatrix(FEModel* model);
    
    /**
     * @brief Get external load at given time
     * @param model FE model
     * @param time Time value
     * @param load Output load vector
     */
    void getExternalLoad(FEModel* model, double time, std::vector<double>& load);
    
    /**
     * @brief Apply boundary conditions
     * @param model FE model
     */
    void applyBoundaryConditions(FEModel* model);
    
    // Newmark parameters
    double beta_;     ///< Newmark beta parameter (default: 0.25 for average acceleration)
    double gamma_;    ///< Newmark gamma parameter (default: 0.5)
    
    // Damping parameters
    double alphaM_;   ///< Mass proportional damping coefficient
    double betaK_;    ///< Stiffness proportional damping coefficient
    
    // Newmark integration coefficients
    double a0_, a1_, a2_, a3_, a4_, a5_, a6_, a7_;
    
    // Solution vectors
    std::vector<double> displacement_;      ///< Displacement at current time
    std::vector<double> velocity_;          ///< Velocity at current time
    std::vector<double> acceleration_;      ///< Acceleration at current time
    
    std::vector<double> displacement_old_;  ///< Displacement at previous time
    std::vector<double> velocity_old_;      ///< Velocity at previous time
    std::vector<double> acceleration_old_;  ///< Acceleration at previous time
    
    // System matrices and vectors (placeholders for actual sparse matrices)
    std::vector<double> externalForce_;     ///< External force vector
    std::vector<double> effectiveForce_;    ///< Effective force vector
    
    double currentTimeStep_;                ///< Current time step size
};

#endif // DYNAMICIMPLICITSTEP_H
