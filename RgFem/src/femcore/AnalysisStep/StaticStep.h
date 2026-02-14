/*********************************************************************
 * \file   StaticStep.h
 * \brief  Static analysis step (integrated with FEStaticSolver)
 *
 * \author Integration
 * \date   February 2025
 *********************************************************************/

#pragma once
#include "AnalysisStep.h"
#include <vector>

 // Forward declarations
class FEStaticSolver;
class FENewtonSolver;

/**
 * @brief Static analysis step
 *
 * Responsible for:
 * - Load step advancement and control
 * - Adaptive load stepping
 * - Calling FEStaticSolver to perform actual Newton-Raphson solution
 *
 * Does NOT contain Newton-Raphson algorithm (that's in FEStaticSolver).
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
    double getLoadMultiplier() const { return loadMultiplier_; }

    void enableNonlinear(bool enable) { nonlinear_ = enable; }
    bool isNonlinear() const { return nonlinear_; }

    // Get FEStaticSolver (casted from base solver)
    FEStaticSolver* getStaticSolver() const;
    FENewtonSolver* getNewtonSolver() const;

    // Get results from solver
    std::vector<double> getDisplacement() const;

private:
    /**
     * @brief Execute linear static analysis (single solve)
     * @param model FE model
     * @return true if successful
     */
    bool executeLinear(FEModel* model);

    /**
     * @brief Execute nonlinear static analysis with load stepping
     * @param model FE model
     * @return true if successful
     */
    bool executeNonlinear(FEModel* model);

    /**
     * @brief Perform one load increment
     * @param model FE model
     * @param time Current pseudo-time (load factor)
     * @return true if successful
     */
    bool performIncrement(FEModel* model, double time);

    /**
     * @brief Handle non-convergence (adaptive cutback)
     * @return true to retry, false to give up
     */
    bool handleNonConvergence();

    // Configuration
    double loadMultiplier_;              ///< Total load multiplier
    bool nonlinear_;                     ///< Enable nonlinear analysis

    // Load stepping for nonlinear analysis
    double currentTimeIncrement_;        ///< Current time/load increment size
};
