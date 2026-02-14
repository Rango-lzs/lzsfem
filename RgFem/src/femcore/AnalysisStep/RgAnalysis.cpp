/*********************************************************************
 * \file   RgAnalysis.cpp
 * \brief  Analysis manager implementation
 *
 * \author Integration
 * \date   February 2025
 *********************************************************************/

#include "RgAnalysis.h"
#include "StaticStep.h"
#include "femcore/NewtonSolver/StaticSolver.h"
#include "femcore/NewtonSolver/NewtonSolver.h"
#include "femcore/FEModel.h"
#include "logger/log.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>

 // Platform-specific includes
#ifdef _WIN32
#include <direct.h>
#define CREATE_DIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define CREATE_DIR(path) mkdir(path, 0755)
#endif
#include <iosfwd>
#include <iomanip>

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
    , autoRestart_(false)
    , maxRestartAttempts_(3)
    , currentRestartAttempt_(0)
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
    solvers_.clear();
    solverMap_.clear();
}

void RgAnalysis::setModel(FEModel* model)
{
    model_ = model;
}

// Solver creation and management
FEStaticSolver* RgAnalysis::createStaticSolver()
{
    auto solver = std::make_unique<FEStaticSolver>();
    FEStaticSolver* ptr = solver.get();

    // Set model for solver
    if (model_) {
        ptr->SetFEModel(model_);
    }

    std::string name = "StaticSolver_" + std::to_string(solvers_.size());
    solverMap_[name] = ptr;
    solvers_.push_back(std::move(solver));

    RgLog("Created %s\n", name.c_str());
    return ptr;
}

FESolver* RgAnalysis::getSolver(const std::string& name)
{
    auto it = solverMap_.find(name);
    if (it != solverMap_.end()) {
        return it->second;
    }
    return nullptr;
}

// Step management
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

// Convenience methods for creating analysis types
StaticStep* RgAnalysis::createStaticAnalysis(const std::string& stepName,
    double loadFactor)
{
    // Create static solver
    auto* solver = createStaticSolver();

    // Configure solver for linear analysis
    FENewtonSolver* newtonSolver = dynamic_cast<FENewtonSolver*>(solver);
    if (newtonSolver) {
        newtonSolver->m_maxref = 1;  // Only one reformation for linear
        newtonSolver->m_Rtol = 1e-6;
        newtonSolver->m_Etol = 1e-6;
    }

    // Create step
    auto step = std::make_shared<StaticStep>(stepName);
    step->setLoadMultiplier(loadFactor);
    step->setSolver(solver);
    step->enableNonlinear(false);

    addStep(step);
    return step.get();
}

StaticStep* RgAnalysis::createNonlinearStaticAnalysis(const std::string& stepName,
    double totalTime,
    double timeIncrement)
{
    // Create static solver
    auto* solver = createStaticSolver();

    // Configure solver for nonlinear analysis
    FENewtonSolver* newtonSolver = dynamic_cast<FENewtonSolver*>(solver);
    if (newtonSolver) {
        newtonSolver->m_maxref = 15;
        newtonSolver->m_Rtol = 1e-5;
        newtonSolver->m_Etol = 1e-4;
        newtonSolver->m_breformtimestep = true;
        newtonSolver->m_bdivreform = true;

        // Enable line search
        if (newtonSolver->m_lineSearch) {
            newtonSolver->m_lineSearch->m_LSiter = 5;
        }
    }

    // Create step
    auto step = std::make_shared<StaticStep>(stepName);
    step->enableNonlinear(true);
    step->setSolver(solver);

    StepControl ctrl;
    ctrl.totalTime = totalTime;
    ctrl.initialTimeIncrement = timeIncrement;
    ctrl.minimumTimeIncrement = 1e-6;
    ctrl.maximumTimeIncrement = 0.5;
    ctrl.useAdaptiveTimeStep = true;
    ctrl.cutbackFactor = 0.5;
    ctrl.increaseFactorGood = 1.25;
    ctrl.minIterationsForIncrease = 3;
    ctrl.maxIterationsForCutback = 8;
    step->setStepControl(ctrl);

    addStep(step);
    return step.get();
}

// Analysis execution
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
    std::cout << "Model: " << model_->GetName() << std::endl;
    std::cout << "Number of steps: " << steps_.size() << std::endl;
    std::cout << "Output directory: " << outputDir_ << std::endl;

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
                RgLog("\nAttempting restart (attempt %d of %d)...\n",
                    currentRestartAttempt_, maxRestartAttempts_);
                i--;
                continue;
            }
            else {
                reportError("Step " + std::to_string(i + 1) + " failed");
                success = false;
                break;
            }
        }

        currentRestartAttempt_ = 0;
    }

    finalizeAnalysis();

    running_ = false;

    std::cout << "\n" << std::string(80, '=') << std::endl;
    if (success) {
        std::cout << "Analysis Completed Successfully" << std::endl;
    }
    else {
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
        }
        });

    try {
        // Initialize step
        step->initialize(model_);

        // Execute step
        bool success = step->execute(model_);

        // Finalize step
        step->finalize(model_);

        return success;
    }
    catch (const std::exception& e) {
        reportError(std::string("Exception in step: ") + e.what());
        return false;
    }
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
        RgLog("\nAnalysis paused.\n");
    }
}

void RgAnalysis::resume()
{
    if (paused_) {
        paused_ = false;
        RgLog("Analysis resumed.\n");
    }
}

void RgAnalysis::stop()
{
    stopRequested_ = true;
    RgLog("\nStop requested.\n");
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
    CREATE_DIR(outputDir_.c_str());

    if (outputResults_) {
        RgLog("Output directory created: %s\n", outputDir_.c_str());
    }
}

void RgAnalysis::finalizeAnalysis()
{
    if (outputResults_) {
        std::string summaryFile = outputDir_ + "/analysis_summary.txt";
        std::ofstream out(summaryFile);

        if (out.is_open()) {
            out << "Analysis Summary\n";
            out << "================\n\n";
            out << "Model: " << model_->GetName() << "\n";
            out << "Number of steps: " << steps_.size() << "\n";
            out << "Completed steps: " << currentStep_ + 1 << "\n\n";

            for (size_t i = 0; i < steps_.size() && i <= static_cast<size_t>(currentStep_); ++i) {
                out << "Step " << (i + 1) << ": "
                    << steps_[i]->getStepName() << "\n";
            }

            out.close();
            RgLog("Analysis summary written to: %s\n", summaryFile.c_str());
        }
    }
}

void RgAnalysis::reportProgress(const std::string& message)
{
    if (progressCallback_) {
        progressCallback_(currentStep_,
            steps_.size() > 0 ? steps_[currentStep_]->getCurrentIncrement() : 0,
            currentTime_,
            message);
    }
}

void RgAnalysis::reportError(const std::string& message)
{
    lastError_ = message;
    RgLogError("%s\n", message.c_str());

    if (errorCallback_) {
        errorCallback_(message);
    }
}