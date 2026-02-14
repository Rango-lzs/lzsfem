#ifndef RGANALYSIS_H
#define RGANALYSIS_H

#include "AnalysisStep.h"
#include "StepControl.h"
#include <memory>
#include <vector>
#include <string>
#include <functional>

// Forward declarations
class FEModel;
class StaticStep;
class DynamicImplicitStep;
class ModalStep;

/**
 * @brief Main analysis manager class for finite element simulations
 * 
 * This class manages the execution of multi-step finite element analyses.
 * It provides convenient methods for creating common analysis types and
 * controls the overall analysis workflow.
 */
class RgAnalysis {
public:
    /**
     * @brief Default constructor
     */
    RgAnalysis();
    
    /**
     * @brief Constructor with model
     * @param model Pointer to finite element model
     */
    explicit RgAnalysis(FEModel* model);
    
    /**
     * @brief Destructor
     */
    ~RgAnalysis();
    
    // Model management
    /**
     * @brief Set the finite element model
     * @param model Pointer to FE model
     */
    void setModel(FEModel* model);
    
    /**
     * @brief Get the finite element model
     * @return Pointer to FE model
     */
    FEModel* getModel() const { return model_; }
    
    // Step management
    /**
     * @brief Add an analysis step
     * @param step Shared pointer to analysis step
     */
    void addStep(std::shared_ptr<AnalysisStep> step);
    
    /**
     * @brief Insert step at specific position
     * @param index Position to insert
     * @param step Step to insert
     */
    void insertStep(size_t index, std::shared_ptr<AnalysisStep> step);
    
    /**
     * @brief Remove step at index
     * @param index Step index to remove
     */
    void removeStep(size_t index);
    
    /**
     * @brief Clear all steps
     */
    void clearSteps();
    
    /**
     * @brief Get number of steps
     * @return Number of steps
     */
    size_t getNumberOfSteps() const { return steps_.size(); }
    
    /**
     * @brief Get step at index
     * @param index Step index
     * @return Shared pointer to step
     */
    std::shared_ptr<AnalysisStep> getStep(size_t index) const;
    
    // Convenience methods for creating common analysis types
    /**
     * @brief Create a linear static analysis step
     * @param stepName Name of the step
     * @param loadFactor Load multiplier
     * @return Pointer to created step
     */
    StaticStep* createStaticAnalysis(const std::string& stepName = "Static",
                                     double loadFactor = 1.0);
    
    /**
     * @brief Create a nonlinear static analysis step
     * @param stepName Name of the step
     * @param totalTime Total pseudo-time/load factor
     * @param timeIncrement Initial time increment
     * @return Pointer to created step
     */
    StaticStep* createNonlinearStaticAnalysis(const std::string& stepName = "Nonlinear Static",
                                              double totalTime = 1.0,
                                              double timeIncrement = 0.1);
    
    /**
     * @brief Create a dynamic implicit analysis step
     * @param stepName Name of the step
     * @param totalTime Total simulation time
     * @param timeIncrement Time step size
     * @param beta Newmark beta parameter
     * @param gamma Newmark gamma parameter
     * @return Pointer to created step
     */
    DynamicImplicitStep* createDynamicAnalysis(const std::string& stepName = "Dynamic",
                                               double totalTime = 1.0,
                                               double timeIncrement = 0.01,
                                               double beta = 0.25,
                                               double gamma = 0.5);
    
    /**
     * @brief Create a modal analysis step
     * @param stepName Name of the step
     * @param numModes Number of modes to extract
     * @return Pointer to created step
     */
    ModalStep* createModalAnalysis(const std::string& stepName = "Modal",
                                   int numModes = 10);
    
    // Analysis execution
    /**
     * @brief Run all analysis steps
     * @return true if all steps completed successfully
     */
    bool run();
    
    /**
     * @brief Run specific step by index
     * @param stepIndex Index of step to run
     * @return true if successful
     */
    bool runStep(size_t stepIndex);
    
    /**
     * @brief Run steps in a range
     * @param startIndex First step to run
     * @param endIndex Last step to run (inclusive)
     * @return true if all steps successful
     */
    bool runSteps(size_t startIndex, size_t endIndex);
    
    /**
     * @brief Pause the analysis
     */
    void pause();
    
    /**
     * @brief Resume the analysis
     */
    void resume();
    
    /**
     * @brief Stop the analysis
     */
    void stop();
    
    /**
     * @brief Reset analysis to initial state
     */
    void reset();
    
    // Output control
    /**
     * @brief Set output frequency (write every N increments)
     * @param freq Output frequency
     */
    void setOutputFrequency(int freq) { outputFrequency_ = freq; }
    
    /**
     * @brief Get output frequency
     * @return Output frequency
     */
    int getOutputFrequency() const { return outputFrequency_; }
    
    /**
     * @brief Set output directory for results
     * @param dir Directory path
     */
    void setOutputDirectory(const std::string& dir) { outputDir_ = dir; }
    
    /**
     * @brief Get output directory
     * @return Directory path
     */
    const std::string& getOutputDirectory() const { return outputDir_; }
    
    /**
     * @brief Enable/disable result output
     * @param enable true to enable output
     */
    void enableResultOutput(bool enable) { outputResults_ = enable; }
    
    /**
     * @brief Check if output is enabled
     * @return true if enabled
     */
    bool isResultOutputEnabled() const { return outputResults_; }
    
    /**
     * @brief Set result file name prefix
     * @param prefix File name prefix
     */
    void setResultFilePrefix(const std::string& prefix) { resultFilePrefix_ = prefix; }
    
    // Callback for progress monitoring
    using ProgressCallback = std::function<void(int step, int increment, double time, const std::string& message)>;
    
    /**
     * @brief Set progress callback function
     * @param callback Callback function
     */
    void setProgressCallback(ProgressCallback callback) { progressCallback_ = callback; }
    
    // Error handling
    using ErrorCallback = std::function<void(const std::string& errorMessage)>;
    
    /**
     * @brief Set error callback function
     * @param callback Callback function
     */
    void setErrorCallback(ErrorCallback callback) { errorCallback_ = callback; }
    
    // Analysis state
    /**
     * @brief Check if analysis is currently running
     * @return true if running
     */
    bool isRunning() const { return running_; }
    
    /**
     * @brief Check if analysis is paused
     * @return true if paused
     */
    bool isPaused() const { return paused_; }
    
    /**
     * @brief Get current step index
     * @return Current step index
     */
    int getCurrentStep() const { return currentStep_; }
    
    /**
     * @brief Get current analysis time
     * @return Current time
     */
    double getCurrentTime() const { return currentTime_; }
    
    /**
     * @brief Get last error message
     * @return Error message string
     */
    const std::string& getLastError() const { return lastError_; }
    
    // Configuration
    /**
     * @brief Enable/disable automatic restart on failure
     * @param enable true to enable
     */
    void enableAutoRestart(bool enable) { autoRestart_ = enable; }
    
    /**
     * @brief Set maximum number of restart attempts
     * @param maxAttempts Maximum attempts
     */
    void setMaxRestartAttempts(int maxAttempts) { maxRestartAttempts_ = maxAttempts; }
    
    /**
     * @brief Enable/disable checkpointing
     * @param enable true to enable
     */
    void enableCheckpointing(bool enable) { useCheckpointing_ = enable; }
    
    /**
     * @brief Set checkpoint frequency
     * @param freq Checkpoint every N increments
     */
    void setCheckpointFrequency(int freq) { checkpointFrequency_ = freq; }
    
    /**
     * @brief Save checkpoint
     * @param filename Checkpoint file name
     * @return true if successful
     */
    bool saveCheckpoint(const std::string& filename);
    
    /**
     * @brief Load checkpoint
     * @param filename Checkpoint file name
     * @return true if successful
     */
    bool loadCheckpoint(const std::string& filename);
    
private:
    /**
     * @brief Initialize analysis before execution
     */
    void initializeAnalysis();
    
    /**
     * @brief Finalize analysis after execution
     */
    void finalizeAnalysis();
    
    /**
     * @brief Write results to file
     * @param stepIndex Step index
     * @param increment Increment number
     */
    void writeResults(int stepIndex, int increment);
    
    /**
     * @brief Create output directory if it doesn't exist
     */
    void createOutputDirectory();
    
    /**
     * @brief Report progress through callback
     * @param message Progress message
     */
    void reportProgress(const std::string& message = "");
    
    /**
     * @brief Report error through callback
     * @param message Error message
     */
    void reportError(const std::string& message);
    
    // Model
    FEModel* model_;                                  ///< Finite element model
    
    // Analysis steps
    std::vector<std::shared_ptr<AnalysisStep>> steps_; ///< List of analysis steps
    
    // Analysis state
    int currentStep_;                                  ///< Current step index
    double currentTime_;                               ///< Current analysis time
    bool running_;                                     ///< Analysis is running
    bool paused_;                                      ///< Analysis is paused
    bool stopRequested_;                               ///< Stop has been requested
    std::string lastError_;                            ///< Last error message
    
    // Output settings
    int outputFrequency_;                              ///< Output frequency
    std::string outputDir_;                            ///< Output directory
    bool outputResults_;                               ///< Enable result output
    std::string resultFilePrefix_;                     ///< Result file prefix
    
    // Restart and recovery
    bool autoRestart_;                                 ///< Auto restart on failure
    int maxRestartAttempts_;                           ///< Max restart attempts
    int currentRestartAttempt_;                        ///< Current restart attempt
    bool useCheckpointing_;                            ///< Use checkpointing
    int checkpointFrequency_;                          ///< Checkpoint frequency
    
    // Callbacks
    ProgressCallback progressCallback_;                ///< Progress callback
    ErrorCallback errorCallback_;                      ///< Error callback
};

#endif // RGANALYSIS_H
