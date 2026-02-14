#ifndef STEPCONTROL_H
#define STEPCONTROL_H

/**
 * @brief Step control parameters for analysis steps
 */
struct StepControl {
    // Time stepping parameters
    double initialTimeIncrement;     ///< Initial time/load increment
    double minimumTimeIncrement;     ///< Minimum allowed time increment
    double maximumTimeIncrement;     ///< Maximum allowed time increment
    double totalTime;                ///< Total analysis time/load factor
    
    // Iteration control
    int maxIterations;               ///< Maximum iterations per increment
    double convergenceTolerance;     ///< Convergence tolerance
    
    // Adaptive control
    bool useAdaptiveTimeStep;        ///< Enable adaptive time stepping
    double cutbackFactor;            ///< Factor for reducing time step on non-convergence
    double increaseFactorGood;       ///< Factor for increasing time step on good convergence
    int minIterationsForIncrease;    ///< Min iterations to allow time step increase
    int maxIterationsForCutback;     ///< Max iterations before cutback
    
    // Convergence criteria
    bool checkDisplacementNorm;      ///< Check displacement norm convergence
    bool checkForceNorm;             ///< Check force residual norm convergence
    bool checkEnergyNorm;            ///< Check energy norm convergence
    double displacementTolerance;    ///< Displacement convergence tolerance
    double forceTolerance;           ///< Force convergence tolerance
    double energyTolerance;          ///< Energy convergence tolerance
    
    /**
     * @brief Default constructor with reasonable defaults
     */
    StepControl()
        : initialTimeIncrement(0.01)
        , minimumTimeIncrement(1e-8)
        , maximumTimeIncrement(1.0)
        , totalTime(1.0)
        , maxIterations(100)
        , convergenceTolerance(1e-6)
        , useAdaptiveTimeStep(true)
        , cutbackFactor(0.5)
        , increaseFactorGood(1.25)
        , minIterationsForIncrease(3)
        , maxIterationsForCutback(8)
        , checkDisplacementNorm(true)
        , checkForceNorm(true)
        , checkEnergyNorm(false)
        , displacementTolerance(1e-6)
        , forceTolerance(1e-6)
        , energyTolerance(1e-8)
    {}
};

/**
 * @brief Analysis type enumeration
 */
enum class AnalysisType {
    STATIC,                      ///< Linear/nonlinear static analysis
    DYNAMIC_IMPLICIT,            ///< Implicit dynamic analysis
    DYNAMIC_EXPLICIT,            ///< Explicit dynamic analysis
    MODAL,                       ///< Modal/frequency analysis
    BUCKLING,                    ///< Buckling analysis
    THERMAL,                     ///< Thermal analysis
    COUPLED_THERMO_MECHANICAL,   ///< Coupled thermo-mechanical
    NONLINEAR_STATIC            ///< Nonlinear static with geometric nonlinearity
};

#endif // STEPCONTROL_H
