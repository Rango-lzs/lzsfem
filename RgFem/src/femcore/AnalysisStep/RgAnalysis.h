/*********************************************************************
 * \file   RgAnalysis.h
 * \brief  Analysis manager (integrated with FESolver)
 *
 * \author Integration
 * \date   February 2025
 *********************************************************************/

#pragma once
#include "AnalysisStep.h"
#include "StepControl.h"
#include <memory>
#include <vector>
#include <string>
#include <functional>
#include <map>

 // Forward declarations
class FEModel;
class FESolver;
class FEStaticSolver;
class FENewtonSolver;
class StaticStep;

/**
 * @brief Main analysis manager class
 *
 * Responsibilities:
 * 1. Manage FEModel
 * 2. Manage analysis steps
 * 3. Manage FESolver instances and assign them to steps
 * 4. Control overall analysis workflow
 * 5. Handle output and checkpointing
 */
class RgAnalysis {
public:
    /**
     * @brief Default constructor
     */
    RgAnalysis();

    /**
     * @brief Constructor with model
     */
    explicit RgAnalysis(FEModel* model);

    /**
     * @brief Destructor
     */
    ~RgAnalysis();

    // Model management
    void setModel(FEModel* model);
    FEModel* getModel() const { return model_; }

    // Solver management
    /**
     * @brief Create and register a static solver
     * @return Pointer to created FEStaticSolver
     */
    FEStaticSolver* createStaticSolver();

    /**
     * @brief Get solver by name
     * @param name Solver name
     * @return Pointer to solver or nullptr
     */
    FESolver* getSolver(const std::string& name);

    // Step management
    void addStep(std::shared_ptr<AnalysisStep> step);
    void insertStep(size_t index, std::shared_ptr<AnalysisStep> step);
    void removeStep(size_t index);
    void clearSteps();

    size_t getNumberOfSteps() const { return steps_.size(); }
    std::shared_ptr<AnalysisStep> getStep(size_t index) const;

    // Convenience methods for creating common analysis types
    /**
     * @brief Create linear static analysis with automatic solver creation
     * @param stepName Step name
     * @param loadFactor Load multiplier
     * @return Pointer to created step
     */
    StaticStep* createStaticAnalysis(const std::string& stepName = "Static",
        double loadFactor = 1.0);

    /**
     * @brief Create nonlinear static analysis with automatic solver creation
     * @param stepName Step name
     * @param totalTime Total pseudo-time/load factor
     * @param timeIncrement Initial time increment
     * @return Pointer to created step
     */
    StaticStep* createNonlinearStaticAnalysis(const std::string& stepName = "Nonlinear Static",
        double totalTime = 1.0,
        double timeIncrement = 0.1);

    // Analysis execution
    bool run();
    bool runStep(size_t stepIndex);
    bool runSteps(size_t startIndex, size_t endIndex);

    void pause();
    void resume();
    void stop();
    void reset();

    // Output control
    void setOutputFrequency(int freq) { outputFrequency_ = freq; }
    int getOutputFrequency() const { return outputFrequency_; }

    void setOutputDirectory(const std::string& dir) { outputDir_ = dir; }
    const std::string& getOutputDirectory() const { return outputDir_; }

    void enableResultOutput(bool enable) { outputResults_ = enable; }
    bool isResultOutputEnabled() const { return outputResults_; }

    // Callbacks
    using ProgressCallback = std::function<void(int step, int increment, double time, const std::string& message)>;
    void setProgressCallback(ProgressCallback callback) { progressCallback_ = callback; }

    using ErrorCallback = std::function<void(const std::string& errorMessage)>;
    void setErrorCallback(ErrorCallback callback) { errorCallback_ = callback; }

    // Analysis state
    bool isRunning() const { return running_; }
    bool isPaused() const { return paused_; }
    int getCurrentStep() const { return currentStep_; }
    double getCurrentTime() const { return currentTime_; }
    const std::string& getLastError() const { return lastError_; }

    // Configuration
    void enableAutoRestart(bool enable) { autoRestart_ = enable; }
    void setMaxRestartAttempts(int maxAttempts) { maxRestartAttempts_ = maxAttempts; }

private:
    void initializeAnalysis();
    void finalizeAnalysis();
    void reportProgress(const std::string& message = "");
    void reportError(const std::string& message);

    // Model
    FEModel* model_;

    // Solvers (owned by analysis)
    std::vector<std::unique_ptr<FESolver>> solvers_;
    std::map<std::string, FESolver*> solverMap_;

    // Analysis steps
    std::vector<std::shared_ptr<AnalysisStep>> steps_;

    // Analysis state
    int currentStep_;
    double currentTime_;
    bool running_;
    bool paused_;
    bool stopRequested_;
    std::string lastError_;

    // Output settings
    int outputFrequency_;
    std::string outputDir_;
    bool outputResults_;

    // Restart and recovery
    bool autoRestart_;
    int maxRestartAttempts_;
    int currentRestartAttempt_;

    // Callbacks
    ProgressCallback progressCallback_;
    ErrorCallback errorCallback_;
};
