# RgAnalysis - Finite Element Analysis Framework

A modern C++ framework for finite element analysis with support for multiple analysis types and multi-step simulations.

## Overview

RgAnalysis is a refactored and modular finite element analysis framework that provides:

- **Multiple analysis types**: Static, Nonlinear Static, Dynamic Implicit, Modal
- **Flexible multi-step workflows**: Chain multiple analysis steps together
- **Advanced convergence control**: Adaptive time stepping, multiple convergence criteria
- **Progress monitoring**: Callbacks for GUI integration
- **Checkpointing and restart**: Save and resume analyses
- **Clean architecture**: Polymorphic step design for easy extension

## File Structure

```
RgAnalysis Framework
├── Core Classes
│   ├── StepControl.h                 # Step control parameters and enums
│   ├── AnalysisStep.h                # Base class for all analysis steps
│   ├── RgAnalysis.h/cpp              # Main analysis manager
│   
├── Analysis Steps
│   ├── StaticStep.h/cpp              # Linear/nonlinear static analysis
│   ├── DynamicImplicitStep.h/cpp     # Implicit dynamic (Newmark)
│   └── ModalStep.h/cpp               # Modal/frequency analysis
│
└── Examples
    └── ExampleUsage.cpp              # Usage examples
```

## Quick Start

### 1. Linear Static Analysis

```cpp
#include "RgAnalysis.h"

FEModel* model = new FEModel();
// ... set up model ...

RgAnalysis analysis(model);
analysis.createStaticAnalysis("Static Load", 1.0);
analysis.run();
```

### 2. Nonlinear Static Analysis

```cpp
RgAnalysis analysis(model);

auto* step = analysis.createNonlinearStaticAnalysis(
    "Nonlinear",
    1.0,    // Total load factor
    0.1     // Initial increment
);

// Customize convergence
StepControl ctrl = step->getStepControl();
ctrl.maxIterations = 50;
ctrl.convergenceTolerance = 1e-6;
step->setStepControl(ctrl);

analysis.run();
```

### 3. Dynamic Analysis

```cpp
RgAnalysis analysis(model);

auto* step = analysis.createDynamicAnalysis(
    "Dynamic",
    10.0,   // 10 seconds
    0.01,   // 0.01 s time step
    0.25,   // Newmark beta
    0.5     // Newmark gamma
);

// Set damping
step->setRayleighDamping(0.1, 0.001);

// Or from frequencies
step->setDampingFromFrequencies(1.0, 10.0, 0.05, 0.05);

analysis.run();
```

### 4. Modal Analysis

```cpp
RgAnalysis analysis(model);

auto* step = analysis.createModalAnalysis("Modal", 20);
step->setFrequencyRange(0.0, 100.0);
step->setNormalizationType(ModalStep::NormalizationType::MASS);

analysis.run();

// Get results
const auto& frequencies = step->getFrequencies();
```

### 5. Multi-Step Analysis

```cpp
RgAnalysis analysis(model);

// Step 1: Modal analysis
analysis.createModalAnalysis("Modal", 10);

// Step 2: Dynamic analysis
analysis.createDynamicAnalysis("Dynamic", 20.0, 0.01);

// Run all steps
analysis.run();
```

## Analysis Types

### Static Analysis

**Linear Static**
- Direct solution of K*u = F
- Single load step

**Nonlinear Static**
- Newton-Raphson iteration
- Adaptive time stepping
- Line search for improved convergence
- Geometric nonlinearity support

**Features:**
- Load multiplier control
- Multiple convergence criteria
- Automatic cutback on non-convergence

### Dynamic Implicit Analysis

**Newmark Time Integration**
- Unconditionally stable (average acceleration)
- Rayleigh damping
- Mass and stiffness proportional damping

**Features:**
- Initial conditions (displacement, velocity, acceleration)
- Time-dependent loading
- Configurable beta and gamma parameters

### Modal Analysis

**Eigenvalue Solvers**
- Subspace Iteration
- Lanczos (planned)
- Jacobi-Davidson (planned)

**Features:**
- Mass normalization
- Modal participation factors
- Effective modal masses
- Frequency range filtering

## Step Control Parameters

```cpp
StepControl ctrl;

// Time stepping
ctrl.totalTime = 1.0;
ctrl.initialTimeIncrement = 0.01;
ctrl.minimumTimeIncrement = 1e-8;
ctrl.maximumTimeIncrement = 1.0;

// Convergence
ctrl.maxIterations = 100;
ctrl.convergenceTolerance = 1e-6;

// Adaptive control
ctrl.useAdaptiveTimeStep = true;
ctrl.cutbackFactor = 0.5;
ctrl.increaseFactorGood = 1.25;

// Multiple criteria
ctrl.checkDisplacementNorm = true;
ctrl.checkForceNorm = true;
ctrl.checkEnergyNorm = true;
ctrl.displacementTolerance = 1e-6;
ctrl.forceTolerance = 1e-6;
ctrl.energyTolerance = 1e-8;
```

## Advanced Features

### Progress Monitoring

```cpp
analysis.setProgressCallback([](int step, int inc, double time, const std::string& msg) {
    std::cout << "Step " << step << ", Increment " << inc 
              << ", Time " << time << std::endl;
});
```

### Error Handling

```cpp
analysis.setErrorCallback([](const std::string& error) {
    std::cerr << "Error: " << error << std::endl;
});
```

### Checkpointing

```cpp
analysis.enableCheckpointing(true);
analysis.setCheckpointFrequency(100);  // Every 100 increments

// Manual checkpoint
analysis.saveCheckpoint("checkpoint.chk");

// Restart from checkpoint
analysis.loadCheckpoint("checkpoint.chk");
```

### Auto-Restart

```cpp
analysis.enableAutoRestart(true);
analysis.setMaxRestartAttempts(3);
```

### Output Control

```cpp
analysis.setOutputDirectory("./results");
analysis.setOutputFrequency(10);  // Output every 10 increments
analysis.setResultFilePrefix("result");
analysis.enableResultOutput(true);
```

## Architecture

### Class Hierarchy

```
AnalysisStep (abstract base)
    ├── StaticStep
    ├── DynamicImplicitStep
    ├── ModalStep
    └── [Future: DynamicExplicitStep, ThermalStep, etc.]

RgAnalysis (manager)
    └── manages vector of AnalysisStep*
```

### Adding New Analysis Types

To add a new analysis type:

1. Inherit from `AnalysisStep`
2. Implement required virtual methods:
   - `initialize()`
   - `execute()`
   - `finalize()`
   - `getAnalysisType()`
3. Add factory method to `RgAnalysis` (optional)

Example:

```cpp
class BucklingStep : public AnalysisStep {
public:
    bool execute(FEModel* model) override {
        // Solve (K - lambda*Kg)*phi = 0
        return true;
    }
    
    void initialize(FEModel* model) override {
        // Setup
    }
    
    void finalize(FEModel* model) override {
        // Cleanup
    }
    
    AnalysisType getAnalysisType() const override {
        return AnalysisType::BUCKLING;
    }
};
```

## Integration Points

The framework currently uses placeholder implementations for:

- **Matrix Assembly**: `assembleStiffnessMatrix()`, `assembleMassMatrix()`
- **Linear Solver**: Solver interface needed
- **Element Routines**: Element stiffness, mass, internal forces
- **Boundary Conditions**: Application to system

These should be connected to your existing FEModel implementation.

## TODO: Implementation Hooks

To fully integrate with your FE solver, implement:

```cpp
// In StaticStep::assembleStiffnessMatrix()
model->clearStiffnessMatrix();
for (auto& element : model->getMesh()->getElements()) {
    MatrixXd Ke = element->getStiffnessMatrix();
    std::vector<int> dofs = element->getDOFs();
    model->assembleElementMatrix(Ke, dofs);
}

// In StaticStep::solveLinearSystem()
model->getSolver()->solve(K, F, displacement_);

// In StaticStep::applyBoundaryConditions()
for (auto& bc : model->getBoundaryConditions()) {
    // Apply BC using penalty, Lagrange multiplier, or elimination
}
```

## Compilation

```bash
# Compile individual files
g++ -c -std=c++11 RgAnalysis.cpp
g++ -c -std=c++11 StaticStep.cpp
g++ -c -std=c++11 DynamicImplicitStep.cpp
g++ -c -std=c++11 ModalStep.cpp

# Link example
g++ -std=c++11 -o example ExampleUsage.cpp RgAnalysis.o StaticStep.o \
    DynamicImplicitStep.o ModalStep.o -lYourFELibrary
```

Or use your existing build system (CMake, Make, etc.)

## License

[Your license here]

## Contact

[Your contact information]
