/*********************************************************************
 * \file   StaticStep.cpp
 * \brief  Static analysis step implementation
 *
 * \author Integration
 * \date   February 2025
 *********************************************************************/

#include "StaticStep.h"
#include "femcore/NewtonSolver/StaticSolver.h"
#include "femcore/NewtonSolver/NewtonSolver.h"
#include "femcore/FEModel.h"
#include "logger/log.h"
#include <iostream>
#include <cmath>
#include <iomanip>

StaticStep::StaticStep(const std::string& name)
    : loadMultiplier_(1.0)
    , nonlinear_(false)
    , currentTimeIncrement_(0.0)
{
    stepName_ = name;
    currentTime_ = 0.0;
    currentIncrement_ = 0;
    outputFrequency_ = 1;
    solver_ = nullptr;
}

void StaticStep::initialize(FEModel* model)
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "Initializing Static Step: " << stepName_ << std::endl;
    std::cout << "========================================" << std::endl;

    if (nonlinear_) {
        std::cout << "Analysis type: Nonlinear Static" << std::endl;
    }
    else {
        std::cout << "Analysis type: Linear Static" << std::endl;
    }

    // Check if solver is assigned
    if (!solver_) {
        RgLogError("No solver assigned to step!");
        throw std::runtime_error("No solver assigned");
    }

    // Verify solver type
    FEStaticSolver* staticSolver = getStaticSolver();
    if (!staticSolver) {
        RgLogError("Solver is not FEStaticSolver!");
        throw std::runtime_error("Wrong solver type");
    }

    // Initialize solver
    if (!solver_->Init()) {
        RgLogError("Solver initialization failed!");
        throw std::runtime_error("Solver initialization failed");
    }

    currentTime_ = 0.0;
    currentIncrement_ = 0;
    currentTimeIncrement_ = control_.initialTimeIncrement;

    RgLog("Load multiplier: %.4f\n", loadMultiplier_);
    RgLog("Initialization complete.\n");
}

bool StaticStep::execute(FEModel* model)
{
    RgLog("Executing static analysis...\n");

    bool success = false;

    if (nonlinear_) {
        success = executeNonlinear(model);
    }
    else {
        success = executeLinear(model);
    }

    return success;
}

void StaticStep::finalize(FEModel* model)
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "Static Step Completed: " << stepName_ << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Total increments: " << currentIncrement_ << std::endl;
    std::cout << "Final load factor: " << currentTime_ << std::endl;

    // Get statistics from solver
    FENewtonSolver* newtonSolver = getNewtonSolver();
    if (newtonSolver) {
        std::cout << "Total iterations: " << newtonSolver->m_niter << std::endl;
        std::cout << "Total reformations: " << newtonSolver->m_ntotref << std::endl;
    }

    std::cout << std::endl;
}

bool StaticStep::executeLinear(FEModel* model)
{
    RgLog("\n--- Linear Static Analysis ---\n");

    FEStaticSolver* staticSolver = getStaticSolver();

    // Initialize the step in solver
    if (!staticSolver->InitStep(loadMultiplier_)) {
        RgLogError("InitStep failed!");
        return false;
    }

    // Solve single step (solver does Newton iteration internally)
    bool success = staticSolver->SolveStep();

    if (success) {
        RgLog("Linear solution obtained.\n");
        currentIncrement_ = 1;
        currentTime_ = loadMultiplier_;
        reportProgress(staticSolver->m_niter);
    }
    else {
        RgLogError("Linear solve failed!\n");
    }

    return success;
}

bool StaticStep::executeNonlinear(FEModel* model)
{
    RgLog("\n--- Nonlinear Static Analysis ---\n");
    RgLog("Load stepping with adaptive control\n\n");

    FEStaticSolver* staticSolver = getStaticSolver();
    FENewtonSolver* newtonSolver = getNewtonSolver();

    // Print header
    std::cout << std::setw(10) << "Increment"
        << std::setw(12) << "Time"
        << std::setw(10) << "Iters"
        << std::setw(15) << "Res Norm"
        << std::setw(15) << "Energy Norm"
        << std::setw(12) << "dt" << std::endl;
    std::cout << std::string(74, '-') << std::endl;

    currentTime_ = 0.0;
    currentTimeIncrement_ = control_.initialTimeIncrement;

    // Load stepping loop
    while (currentTime_ < control_.totalTime) {
        currentIncrement_++;

        // Calculate time for this increment
        double nextTime = std::min(currentTime_ + currentTimeIncrement_,
            control_.totalTime);

        // Perform increment
        bool success = performIncrement(model, nextTime);

        if (!success) {
            // Handle non-convergence
            if (!handleNonConvergence()) {
                RgLogError("\nAnalysis failed!\n");
                return false;
            }
            // Retry will happen in next loop iteration
            continue;
        }

        // Print results
        std::cout << std::setw(10) << currentIncrement_
            << std::setw(12) << std::fixed << std::setprecision(4) << currentTime_
            << std::setw(10) << newtonSolver->m_niter
            << std::setw(15) << std::scientific << newtonSolver->m_residuNorm.normi
            << std::setw(15) << newtonSolver->m_energyNorm.normi
            << std::setw(12) << std::fixed << currentTimeIncrement_
            << std::endl;

        // Update time
        currentTime_ = nextTime;

        // Adaptive time step increase
        if (control_.useAdaptiveTimeStep &&
            newtonSolver->m_niter <= control_.minIterationsForIncrease &&
            currentTimeIncrement_ < control_.maximumTimeIncrement) {

            currentTimeIncrement_ = std::min(
                currentTimeIncrement_ * control_.increaseFactorGood,
                control_.maximumTimeIncrement
            );
        }

        // Check if done
        if (currentTime_ >= control_.totalTime) {
            break;
        }
    }

    std::cout << std::string(74, '-') << std::endl;
    RgLog("Nonlinear analysis completed successfully.\n");

    return true;
}

bool StaticStep::performIncrement(FEModel* model, double time)
{
    FEStaticSolver* staticSolver = getStaticSolver();

    // Scale load by time factor
    double scaleFactor = time / control_.totalTime * loadMultiplier_;

    // Initialize the step in solver
    if (!staticSolver->InitStep(scaleFactor)) {
        RgLogError("InitStep failed at time %.4f\n", time);
        return false;
    }

    // Solve step (Newton iteration happens inside)
    bool success = staticSolver->SolveStep();

    if (success) {
        reportProgress(staticSolver->m_niter);
    }

    return success;
}

bool StaticStep::handleNonConvergence()
{
    // Try adaptive cutback
    if (control_.useAdaptiveTimeStep &&
        currentTimeIncrement_ > control_.minimumTimeIncrement) {

        currentTimeIncrement_ *= control_.cutbackFactor;
        currentIncrement_--;  // Decrement because we'll retry

        RgLog("*** Cutback: reducing increment to %.6e\n", currentTimeIncrement_);
        return true;  // Retry
    }

    return false;  // Give up
}

FEStaticSolver* StaticStep::getStaticSolver() const
{
    return dynamic_cast<FEStaticSolver*>(solver_);
}

FENewtonSolver* StaticStep::getNewtonSolver() const
{
    return dynamic_cast<FENewtonSolver*>(solver_);
}

std::vector<double> StaticStep::getDisplacement() const
{
    FENewtonSolver* newtonSolver = getNewtonSolver();
    if (newtonSolver) {
        return newtonSolver->GetSolutionVector();
    }
    return std::vector<double>();
}