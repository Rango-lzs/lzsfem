/*********************************************************************
 * \file   AnalysisStep.h
 * \brief  Analysis step base class (integrated with FESolver)
 *
 * \author Integration
 * \date   February 2025
 *********************************************************************/

#pragma once
#include "StepControl.h"
#include <string>
#include <functional>

 // Forward declarations
class FEModel;
class FESolver;

/**
 * @brief Base class for all analysis steps
 *
 * AnalysisStep is responsible for:
 * 1. Managing step data (time, load, boundary conditions)
 * 2. Controlling load step advancement (incrementation, adaptive stepping)
 * 3. Delegating actual solving to FESolver
 *
 * It does NOT contain solving algorithms - those are in FESolver classes.
 */
class AnalysisStep {
public:
    virtual ~AnalysisStep() = default;

    /**
     * @brief Execute the analysis step (load step advancement loop)
     * @param model Pointer to the finite element model
     * @return true if successful, false otherwise
     */
    virtual bool execute(FEModel* model) = 0;

    /**
     * @brief Initialize the step before execution
     * @param model Pointer to the finite element model
     */
    virtual void initialize(FEModel* model) = 0;

    /**
     * @brief Finalize the step after execution
     * @param model Pointer to the finite element model
     */
    virtual void finalize(FEModel* model) = 0;

    /**
     * @brief Get the analysis type
     * @return Analysis type enumeration
     */
    virtual AnalysisType getAnalysisType() const = 0;

    // Solver management
    /**
     * @brief Set the solver for this step
     * @param solver Pointer to FESolver (not owned by step)
     */
    void setSolver(FESolver* solver) { solver_ = solver; }

    /**
     * @brief Get the solver
     * @return Pointer to FESolver
     */
    FESolver* getSolver() const { return solver_; }

    // Step control accessors
    void setStepControl(const StepControl& ctrl) { control_ = ctrl; }
    const StepControl& getStepControl() const { return control_; }

    // Step name accessors
    void setStepName(const std::string& name) { stepName_ = name; }
    const std::string& getStepName() const { return stepName_; }

    // Time/increment accessors
    double getCurrentTime() const { return currentTime_; }
    double getTotalTime() const { return control_.totalTime; }
    int getCurrentIncrement() const { return currentIncrement_; }

    // Output control
    void setOutputFrequency(int freq) { outputFrequency_ = freq; }
    int getOutputFrequency() const { return outputFrequency_; }

    // Progress callback
    using ProgressCallback = std::function<void(int increment, double time, int iteration)>;
    void setProgressCallback(ProgressCallback callback) { progressCallback_ = callback; }

protected:
    FESolver* solver_;                 ///< Solver for this step (not owned)
    StepControl control_;              ///< Step control parameters
    std::string stepName_;             ///< Name of this step
    double currentTime_;               ///< Current analysis time/load factor
    int currentIncrement_;             ///< Current increment number
    int outputFrequency_;              ///< Output frequency (every N increments)
    ProgressCallback progressCallback_; ///< Progress monitoring callback

    /**
     * @brief Report progress to callback if set
     */
    void reportProgress(int iteration = 0) {
        if (progressCallback_) {
            progressCallback_(currentIncrement_, currentTime_, iteration);
        }
    }
};
