#include "RgAnalysis.h"
#include "StaticStep.h"
#include "DynamicImplicitStep.h"
#include "ModalStep.h"
#include "femcore/FEModel.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>

// Platform-specific includes for directory creation
#ifdef _WIN32
#include <direct.h>
#define CREATE_DIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define CREATE_DIR(path) mkdir(path, 0755)
#endif

RgAnalysis::RgAnalysis()
    : model_(nullptr)
    , currentStep_(0)
    , currentTime_(0.0)
    , running_(false)
    , paused_(false)
    , stopRequested_(false)
    , outputFrequency_(1)
    , outputDir_("./results")
    , outputResults_(true)
    , resultFilePrefix_("result")
    , autoRestart_(false)
    , maxRestartAttempts_(3)
    , currentRestartAttempt_(0)
    , useCheckpointing_(false)
    , checkpointFrequency_(10)
{
}

RgAnalysis::RgAnalysis(FEModel* model)
    : RgAnalysis()
{
    model_ = model;
}

RgAnalysis::~RgAnalysis()
{
    clearSteps();
}

void RgAnalysis::setModel(FEModel* model)
{
    model_ = model;
}

void RgAnalysis::addStep(std::shared_ptr<AnalysisStep> step)
{
    if (step) {
        steps_.push_back(step);
        step->setOutputFrequency(outputFrequency_);
    }
}

void RgAnalysis::insertStep(size_t index, std::shared_ptr<AnalysisStep> step)
{
    if (step && index <= steps_.size()) {
        steps_.insert(steps_.begin() + index, step);
        step->setOutputFrequency(outputFrequency_);
    }
}

void RgAnalysis::removeStep(size_t index)
{
    if (index < steps_.size()) {
        steps_.erase(steps_.begin() + index);
    }
}

void RgAnalysis::clearSteps()
{
    steps_.clear();
}

std::shared_ptr<AnalysisStep> RgAnalysis::getStep(size_t index) const
{
    if (index < steps_.size()) {
        return steps_[index];
    }
    return nullptr;
}

StaticStep* RgAnalysis::createStaticAnalysis(const std::string& stepName, 
                                              double loadFactor)
{
    auto step = std::make_shared<StaticStep>(stepName);
    step->setLoadMultiplier(loadFactor);
    
    addStep(step);
    return step.get();
}

StaticStep* RgAnalysis::createNonlinearStaticAnalysis(const std::string& stepName,
                                                       double totalTime,
                                                       double timeIncrement)
{
    auto step = std::make_shared<StaticStep>(stepName);
    step->enableNonlinear(true);
    step->enableLargeDisplacement(true);
    step->enableLineSearch(true);
    
    StepControl ctrl;
    ctrl.totalTime = totalTime;
    ctrl.initialTimeIncrement = timeIncrement;
    ctrl.useAdaptiveTimeStep = true;
    step->setStepControl(ctrl);
    
    addStep(step);
    return step.get();
}

DynamicImplicitStep* RgAnalysis::createDynamicAnalysis(const std::string& stepName,
                                                        double totalTime,
                                                        double timeIncrement,
                                                        double beta,
                                                        double gamma)
{
    auto step = std::make_shared<DynamicImplicitStep>(stepName);
    step->setNewmarkParameters(beta, gamma);
    
    StepControl ctrl;
    ctrl.totalTime = totalTime;
    ctrl.initialTimeIncrement = timeIncrement;
    step->setStepControl(ctrl);
    
    addStep(step);
    return step.get();
}

ModalStep* RgAnalysis::createModalAnalysis(const std::string& stepName,
                                           int numModes)
{
    auto step = std::make_shared<ModalStep>(stepName);
    step->setNumberOfModes(numModes);
    
    addStep(step);
    return step.get();
}

bool RgAnalysis::run()
{
    if (!model_) {
        reportError("No model assigned to analysis");
        return false;
    }
    
    if (steps_.empty()) {
        reportError("No analysis steps defined");
        return false;
    }
    
    running_ = true;
    paused_ = false;
    stopRequested_ = false;
    currentTime_ = 0.0;
    currentRestartAttempt_ = 0;
    
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "Starting Analysis" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
    //std::cout << "Model: " << model_->getName() << std::endl;
    std::cout << "Number of steps: " << steps_.size() << std::endl;
    std::cout << "Output directory: " << outputDir_ << std::endl;
    
    // Print timestamp
    std::time_t now = std::time(nullptr);
    std::cout << "Start time: " << std::ctime(&now);
    std::cout << std::string(80, '=') << std::endl;
    
    initializeAnalysis();
    
    bool success = true;
    for (size_t i = 0; i < steps_.size(); ++i) {
        if (stopRequested_) {
            std::cout << "\nAnalysis stopped by user request." << std::endl;
            success = false;
            break;
        }
        
        currentStep_ = static_cast<int>(i);
        
        reportProgress("Starting step " + std::to_string(i + 1));
        
        bool stepSuccess = runStep(i);
        
        if (!stepSuccess) {
            if (autoRestart_ && currentRestartAttempt_ < maxRestartAttempts_) {
                currentRestartAttempt_++;
                std::cout << "\nAttempting restart (attempt " 
                          << currentRestartAttempt_ << " of " 
                          << maxRestartAttempts_ << ")..." << std::endl;
                i--; // Retry the same step
                continue;
            } else {
                reportError("Step " + std::to_string(i + 1) + " failed");
                success = false;
                break;
            }
        }
        
        currentRestartAttempt_ = 0; // Reset on success
    }
    
    finalizeAnalysis();
    
    running_ = false;
    
    std::cout << "\n" << std::string(80, '=') << std::endl;
    if (success) {
        std::cout << "Analysis Completed Successfully" << std::endl;
    } else {
        std::cout << "Analysis Failed" << std::endl;
    }
    
    now = std::time(nullptr);
    std::cout << "End time: " << std::ctime(&now);
    std::cout << std::string(80, '=') << std::endl;
    
    return success;
}

bool RgAnalysis::runStep(size_t stepIndex)
{
    if (stepIndex >= steps_.size()) {
        reportError("Invalid step index: " + std::to_string(stepIndex));
        return false;
    }
    
    auto& step = steps_[stepIndex];
    
    // Set progress callback for the step
    step->setProgressCallback([this, stepIndex](int inc, double time, int iter) {
        currentTime_ = time;
        
        std::ostringstream msg;
        msg << "Step " << (stepIndex + 1) 
            << ", Increment " << inc 
            << ", Time " << std::fixed << std::setprecision(4) << time;
        
        if (iter > 0) {
            msg << ", Iteration " << iter;
        }
        
        reportProgress(msg.str());
        
        // Check for pause
        while (paused_ && !stopRequested_) {
            // Wait while paused
            // In a real implementation, use proper threading/synchronization
        }
        
        // Write checkpoint if enabled
        if (useCheckpointing_ && inc % checkpointFrequency_ == 0) {
            std::ostringstream filename;
            filename << outputDir_ << "/checkpoint_step" << (stepIndex + 1) 
                     << "_inc" << inc << ".chk";
            saveCheckpoint(filename.str());
        }
    });
    
    // Initialize step
    step->initialize(model_);
    
    // Execute step
    bool success = step->execute(model_);
    
    // Finalize step
    step->finalize(model_);
    
    return success;
}

bool RgAnalysis::runSteps(size_t startIndex, size_t endIndex)
{
    if (startIndex > endIndex || endIndex >= steps_.size()) {
        reportError("Invalid step range");
        return false;
    }
    
    bool success = true;
    for (size_t i = startIndex; i <= endIndex; ++i) {
        if (!runStep(i)) {
            success = false;
            break;
        }
    }
    
    return success;
}

void RgAnalysis::pause()
{
    if (running_) {
        paused_ = true;
        std::cout << "\nAnalysis paused." << std::endl;
    }
}

void RgAnalysis::resume()
{
    if (paused_) {
        paused_ = false;
        std::cout << "Analysis resumed." << std::endl;
    }
}

void RgAnalysis::stop()
{
    stopRequested_ = true;
    std::cout << "\nStop requested. Analysis will terminate after current increment." << std::endl;
}

void RgAnalysis::reset()
{
    currentStep_ = 0;
    currentTime_ = 0.0;
    running_ = false;
    paused_ = false;
    stopRequested_ = false;
    currentRestartAttempt_ = 0;
    lastError_.clear();
}

void RgAnalysis::initializeAnalysis()
{
    createOutputDirectory();
    
    if (outputResults_) {
        std::cout << "Output directory created: " << outputDir_ << std::endl;
        std::cout << "Result files will be written with prefix: " 
                  << resultFilePrefix_ << std::endl;
    }
    
    // TODO: Initialize result database, logging, etc.
}

void RgAnalysis::finalizeAnalysis()
{
    // TODO: Close result files, write summary, etc.
    
    if (outputResults_) {
        std::string summaryFile = outputDir_ + "/analysis_summary.txt";
        std::ofstream out(summaryFile);
        
        if (out.is_open()) {
            out << "Analysis Summary\n";
            out << "================\n\n";
            //out << "Model: " << model_->getName() << "\n";
            out << "Number of steps: " << steps_.size() << "\n";
            out << "Completed steps: " << currentStep_ + 1 << "\n";
            out << "Final time: " << currentTime_ << "\n\n";
            
            for (size_t i = 0; i < steps_.size() && i <= static_cast<size_t>(currentStep_); ++i) {
                out << "Step " << (i + 1) << ": " 
                    << steps_[i]->getStepName() << "\n";
            }
            
            out.close();
            std::cout << "Analysis summary written to: " << summaryFile << std::endl;
        }
    }
}

void RgAnalysis::writeResults(int stepIndex, int increment)
{
    if (!outputResults_) return;
    
    // TODO: Write results to file
    // This would typically write:
    // - Nodal displacements
    // - Element stresses/strains
    // - Reaction forces
    // - Energy values
    // etc.
    
    std::ostringstream filename;
    filename << outputDir_ << "/" << resultFilePrefix_ 
             << "_step" << stepIndex 
             << "_inc" << increment << ".res";
    
    // Actual file writing would go here
}

void RgAnalysis::createOutputDirectory()
{
    // Create output directory if it doesn't exist
    CREATE_DIR(outputDir_.c_str());
    
    // Note: In a production implementation, you'd want to check
    // the return value and handle errors appropriately
}

void RgAnalysis::reportProgress(const std::string& message)
{
    if (progressCallback_) {
        progressCallback_(currentStep_, 
                         steps_[currentStep_]->getCurrentIncrement(),
                         currentTime_,
                         message);
    }
}

void RgAnalysis::reportError(const std::string& message)
{
    lastError_ = message;
    std::cerr << "ERROR: " << message << std::endl;
    
    if (errorCallback_) {
        errorCallback_(message);
    }
}

bool RgAnalysis::saveCheckpoint(const std::string& filename)
{
    std::cout << "Saving checkpoint to: " << filename << std::endl;
    
    // TODO: Implement checkpoint saving
    /*
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        reportError("Failed to open checkpoint file for writing: " + filename);
        return false;
    }
    
    // Write checkpoint data:
    // - Current step and increment
    // - Solution state (displacements, velocities, etc.)
    // - Internal variables (stresses, strains, history variables)
    // - Time and load factor
    
    out.close();
    */
    
    return true;
}

bool RgAnalysis::loadCheckpoint(const std::string& filename)
{
    std::cout << "Loading checkpoint from: " << filename << std::endl;
    
    // TODO: Implement checkpoint loading
    /*
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        reportError("Failed to open checkpoint file for reading: " + filename);
        return false;
    }
    
    // Read checkpoint data and restore analysis state
    
    in.close();
    */
    
    return true;
}
