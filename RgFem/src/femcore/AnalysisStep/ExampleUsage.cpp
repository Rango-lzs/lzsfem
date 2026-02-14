#include "RgAnalysis.h"
#include "StaticStep.h"
#include "DynamicImplicitStep.h"
#include "ModalStep.h"
#include "FEModel.h"
#include <iostream>

/**
 * @brief Example 1: Simple linear static analysis
 */
void example1_LinearStatic()
{
    std::cout << "\n=== Example 1: Linear Static Analysis ===" << std::endl;
    
    // Create model (placeholder - in real use, load from file or create programmatically)
    FEModel* model = new FEModel();
    model->setName("Cantilever Beam");
    
    // Create analysis
    RgAnalysis analysis(model);
    
    // Create a simple static step
    analysis.createStaticAnalysis("Static Load", 1.0);
    
    // Configure output
    analysis.setOutputDirectory("./results/example1");
    analysis.setOutputFrequency(1);
    
    // Run analysis
    bool success = analysis.run();
    
    if (success) {
        std::cout << "Analysis completed successfully!" << std::endl;
    }
    
    delete model;
}

/**
 * @brief Example 2: Nonlinear static analysis with multiple load steps
 */
void example2_NonlinearStatic()
{
    std::cout << "\n=== Example 2: Nonlinear Static Analysis ===" << std::endl;
    
    FEModel* model = new FEModel();
    model->setName("Nonlinear Frame");
    
    RgAnalysis analysis(model);
    
    // Create nonlinear static step with load ramping
    auto* step = analysis.createNonlinearStaticAnalysis(
        "Nonlinear Loading",
        1.0,    // Total load factor
        0.1     // Initial increment
    );
    
    // Customize convergence criteria
    StepControl ctrl = step->getStepControl();
    ctrl.maxIterations = 50;
    ctrl.convergenceTolerance = 1e-6;
    ctrl.useAdaptiveTimeStep = true;
    ctrl.cutbackFactor = 0.5;
    ctrl.increaseFactorGood = 1.25;
    step->setStepControl(ctrl);
    
    // Enable line search for better convergence
    step->enableLineSearch(true);
    
    // Set up progress monitoring
    analysis.setProgressCallback([](int step, int inc, double time, const std::string& msg) {
        std::cout << "Progress: " << msg << std::endl;
    });
    
    // Run analysis
    analysis.run();
    
    delete model;
}

/**
 * @brief Example 3: Dynamic analysis with Newmark integration
 */
void example3_DynamicAnalysis()
{
    std::cout << "\n=== Example 3: Dynamic Analysis ===" << std::endl;
    
    FEModel* model = new FEModel();
    model->setName("Dynamic Frame");
    
    RgAnalysis analysis(model);
    
    // Create dynamic step
    auto* dynStep = analysis.createDynamicAnalysis(
        "Transient Dynamic",
        10.0,    // 10 seconds total
        0.01,    // 0.01 second time step
        0.25,    // Newmark beta (average acceleration)
        0.5      // Newmark gamma
    );
    
    // Set Rayleigh damping
    dynStep->setRayleighDamping(0.1, 0.001);
    
    // Or set damping based on frequencies and damping ratios
    // dynStep->setDampingFromFrequencies(
    //     1.0,   // First mode frequency (Hz)
    //     10.0,  // Second mode frequency (Hz)
    //     0.05,  // 5% damping at first mode
    //     0.05   // 5% damping at second mode
    // );
    
    // Set initial conditions
    int nDOF = model->getNumberOfDOFs();
    std::vector<double> initialVelocity(nDOF, 0.0);
    // Set some initial velocity...
    dynStep->setInitialVelocity(initialVelocity);
    
    // Configure output (output every 10 time steps)
    analysis.setOutputFrequency(10);
    analysis.setOutputDirectory("./results/example3");
    
    analysis.run();
    
    delete model;
}

/**
 * @brief Example 4: Modal analysis
 */
void example4_ModalAnalysis()
{
    std::cout << "\n=== Example 4: Modal Analysis ===" << std::endl;
    
    FEModel* model = new FEModel();
    model->setName("Bridge Structure");
    
    RgAnalysis analysis(model);
    
    // Create modal analysis step
    auto* modalStep = analysis.createModalAnalysis("Modal", 20);
    
    // Set frequency range of interest
    modalStep->setFrequencyRange(0.0, 100.0);  // 0-100 Hz
    
    // Set normalization type
    modalStep->setNormalizationType(ModalStep::NormalizationType::MASS);
    
    // Set solver type
    modalStep->setSolverType(ModalStep::SolverType::SUBSPACE_ITERATION);
    
    // Run analysis
    bool success = analysis.run();
    
    if (success) {
        // Get modal results
        const auto& frequencies = modalStep->getFrequencies();
        const auto& periods = modalStep->getPeriods();
        
        std::cout << "\nModal Analysis Results:" << std::endl;
        std::cout << "Mode\tFrequency (Hz)\tPeriod (s)" << std::endl;
        
        for (size_t i = 0; i < frequencies.size(); ++i) {
            std::cout << (i+1) << "\t" 
                      << frequencies[i] << "\t\t" 
                      << periods[i] << std::endl;
        }
    }
    
    delete model;
}

/**
 * @brief Example 5: Multi-step analysis (Modal + Dynamic)
 */
void example5_MultiStepAnalysis()
{
    std::cout << "\n=== Example 5: Multi-Step Analysis ===" << std::endl;
    
    FEModel* model = new FEModel();
    model->setName("Building Structure");
    
    RgAnalysis analysis(model);
    
    // Step 1: Modal analysis to find natural frequencies
    auto* modalStep = analysis.createModalAnalysis("Modal Analysis", 10);
    modalStep->setFrequencyRange(0.0, 50.0);
    
    // Step 2: Dynamic analysis using modal results for damping
    auto* dynStep = analysis.createDynamicAnalysis(
        "Seismic Response",
        20.0,   // 20 seconds
        0.01    // 10 ms time step
    );
    
    // Configure analysis
    analysis.setOutputDirectory("./results/example5");
    analysis.enableCheckpointing(true);
    analysis.setCheckpointFrequency(100);
    
    // Run all steps
    bool success = analysis.run();
    
    if (success) {
        // Can access results from both steps
        const auto& frequencies = modalStep->getFrequencies();
        const auto& displacement = dynStep->getDisplacement();
        
        std::cout << "\nFirst natural frequency: " << frequencies[0] << " Hz" << std::endl;
        std::cout << "Final displacement vector size: " << displacement.size() << std::endl;
    }
    
    delete model;
}

/**
 * @brief Example 6: Custom step control and error handling
 */
void example6_AdvancedControl()
{
    std::cout << "\n=== Example 6: Advanced Control ===" << std::endl;
    
    FEModel* model = new FEModel();
    model->setName("Complex Structure");
    
    RgAnalysis analysis(model);
    
    // Create nonlinear step with custom control
    auto* step = analysis.createNonlinearStaticAnalysis("Pushover", 2.0, 0.05);
    
    // Detailed convergence control
    StepControl ctrl;
    ctrl.totalTime = 2.0;
    ctrl.initialTimeIncrement = 0.05;
    ctrl.minimumTimeIncrement = 1e-6;
    ctrl.maximumTimeIncrement = 0.2;
    
    ctrl.maxIterations = 100;
    ctrl.useAdaptiveTimeStep = true;
    ctrl.cutbackFactor = 0.5;
    ctrl.increaseFactorGood = 1.5;
    ctrl.minIterationsForIncrease = 3;
    ctrl.maxIterationsForCutback = 10;
    
    // Multiple convergence criteria
    ctrl.checkDisplacementNorm = true;
    ctrl.checkForceNorm = true;
    ctrl.checkEnergyNorm = true;
    ctrl.displacementTolerance = 1e-6;
    ctrl.forceTolerance = 1e-6;
    ctrl.energyTolerance = 1e-8;
    
    step->setStepControl(ctrl);
    
    // Enable auto-restart on failure
    analysis.enableAutoRestart(true);
    analysis.setMaxRestartAttempts(3);
    
    // Set error callback
    analysis.setErrorCallback([](const std::string& error) {
        std::cerr << "Analysis Error: " << error << std::endl;
        // Could send email, write to log file, etc.
    });
    
    // Set progress callback
    analysis.setProgressCallback([](int step, int inc, double time, const std::string& msg) {
        static int lastPercent = -1;
        int percent = static_cast<int>(time / 2.0 * 100.0);
        
        if (percent != lastPercent && percent % 10 == 0) {
            std::cout << "Progress: " << percent << "% complete" << std::endl;
            lastPercent = percent;
        }
    });
    
    analysis.run();
    
    delete model;
}

/**
 * @brief Main function to run all examples
 */
int main()
{
    std::cout << "RgAnalysis Framework Examples\n";
    std::cout << "==============================\n" << std::endl;
    
    try {
        // Uncomment the examples you want to run:
        
        example1_LinearStatic();
        // example2_NonlinearStatic();
        // example3_DynamicAnalysis();
        // example4_ModalAnalysis();
        // example5_MultiStepAnalysis();
        // example6_AdvancedControl();
        
    } catch (const std::exception& e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "\n=== All examples completed ===" << std::endl;
    
    return 0;
}
