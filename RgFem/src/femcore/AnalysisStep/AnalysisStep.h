#ifndef ANALYSISSTEP_H
#define ANALYSISSTEP_H

#include "StepControl.h"
#include <string>
#include <vector>
#include <functional>

class FEModel;

/**
 * @brief Base class for all analysis steps
 * 
 * This abstract class defines the interface for all analysis step types.
 * Each derived class implements specific analysis algorithms.
 */
class AnalysisStep {
public:
    virtual ~AnalysisStep() = default;
    
    /**
     * @brief Execute the analysis step
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
    
    // Step control accessors
    void setStepControl(const StepControl& ctrl) { control_ = ctrl; }
    const StepControl& getStepControl() const { return control_; }
    
    // Step name accessors
    void setStepName(const std::string& name) { stepName_ = name; }
    const std::string& getStepName() const { return stepName_; }
    
    // Time/increment accessors
    double getCurrentTime() const { return currentTime_; }
    int getCurrentIncrement() const { return currentIncrement_; }
    
    // Output control
    void setOutputFrequency(int freq) { outputFrequency_ = freq; }
    int getOutputFrequency() const { return outputFrequency_; }
    
    // Progress callback
    using ProgressCallback = std::function<void(int increment, double time, int iteration)>;
    void setProgressCallback(ProgressCallback callback) { progressCallback_ = callback; }
    
protected:
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

#endif // ANALYSISSTEP_H
