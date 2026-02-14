#ifndef STATICSTEP_H
#define STATICSTEP_H

#include "AnalysisStep.h"
#include <vector>

/**
 * @brief Static analysis step
 * 
 * Implements linear and nonlinear static analysis using:
 * - Direct solution for linear problems
 * - Newton-Raphson iteration for nonlinear problems
 * - Arc-length method for snap-through/snap-back problems
 */
class StaticStep : public AnalysisStep {
public:
    /**
     * @brief Constructor
     * @param name Step name
     */
    explicit StaticStep(const std::string& name = "Static");
    
    /**
     * @brief Destructor
     */
    ~StaticStep() override = default;
    
    // AnalysisStep interface implementation
    bool execute(FEModel* model) override;
    void initialize(FEModel* model) override;
    void finalize(FEModel* model) override;
    AnalysisType getAnalysisType() const override { return AnalysisType::STATIC; }
    
    // Configuration methods
    void setLoadMultiplier(double mult) { loadMultiplier_ = mult; }
    void enableNonlinear(bool enable) { nonlinear_ = enable; }
    void enableLargeDisplacement(bool enable) { largeDisplacement_ = enable; }
    void enableLineSearch(bool enable) { useLineSearch_ = enable; }
    void setLineSearchMaxIter(int maxIter) { lineSearchMaxIter_ = maxIter; }
    
    // Get results
    const std::vector<double>& getDisplacement() const { return displacement_; }
    const std::vector<double>& getReactionForces() const { return reactionForces_; }
    
private:
    /**
     * @brief Solve linear static problem: K*u = F
     * @param model FE model
     * @return true if successful
     */
    bool solveLinearSystem(FEModel* model);
    
    /**
     * @brief Solve nonlinear static problem using Newton-Raphson
     * @param model FE model
     * @return true if successful
     */
    bool solveNonlinearSystem(FEModel* model);
    
    /**
     * @brief Perform one Newton-Raphson iteration
     * @param model FE model
     * @param lambda Current load factor
     * @return true if converged
     */
    bool performNewtonIteration(FEModel* model, double lambda);
    
    /**
     * @brief Assemble global stiffness matrix
     * @param model FE model
     * @param geometric Include geometric stiffness if true
     */
    void assembleStiffnessMatrix(FEModel* model, bool geometric = false);
    
    /**
     * @brief Assemble global load vector
     * @param model FE model
     * @param lambda Load multiplier
     */
    void assembleLoadVector(FEModel* model, double lambda);
    
    /**
     * @brief Calculate residual force vector
     * @param model FE model
     * @param lambda Current load factor
     */
    void calculateResidual(FEModel* model, double lambda);
    
    /**
     * @brief Apply boundary conditions to system
     * @param model FE model
     */
    void applyBoundaryConditions(FEModel* model);
    
    /**
     * @brief Check convergence criteria
     * @param iteration Current iteration number
     * @return true if converged
     */
    bool checkConvergence(int iteration);
    
    /**
     * @brief Update solution with displacement increment
     * @param model FE model
     * @param deltaU Displacement increment
     */
    void updateSolution(FEModel* model, const std::vector<double>& deltaU);
    
    /**
     * @brief Perform line search to find optimal step length
     * @param model FE model
     * @param direction Search direction
     * @return Optimal step length factor (0.0 to 1.0)
     */
    double performLineSearch(FEModel* model, const std::vector<double>& direction);
    
    /**
     * @brief Calculate reaction forces at constrained DOFs
     * @param model FE model
     */
    void calculateReactionForces(FEModel* model);
    
    /**
     * @brief Update stresses and strains in all elements
     * @param model FE model
     */
    void updateStressStrain(FEModel* model);
    
    // Member variables
    double loadMultiplier_;              ///< Load multiplier factor
    bool nonlinear_;                     ///< Enable nonlinear analysis
    bool largeDisplacement_;             ///< Enable geometric nonlinearity
    bool useLineSearch_;                 ///< Use line search in Newton-Raphson
    int lineSearchMaxIter_;              ///< Max iterations for line search
    
    // Solution vectors
    std::vector<double> displacement_;   ///< Global displacement vector
    std::vector<double> internalForce_;  ///< Internal force vector
    std::vector<double> externalForce_;  ///< External force vector
    std::vector<double> residual_;       ///< Residual force vector
    std::vector<double> reactionForces_; ///< Reaction forces
    
    // Convergence tracking
    double displacementNorm_;            ///< Displacement increment norm
    double forceNorm_;                   ///< Force residual norm
    double energyNorm_;                  ///< Energy norm
    double initialForceNorm_;            ///< Initial force norm for relative check
};

#endif // STATICSTEP_H
