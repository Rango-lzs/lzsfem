#ifndef MODALSTEP_H
#define MODALSTEP_H

#include "AnalysisStep.h"
#include <vector>

/**
 * @brief Modal/frequency analysis step
 * 
 * Solves the generalized eigenvalue problem:
 * (K - omega^2 * M) * phi = 0
 * 
 * to find natural frequencies and mode shapes
 */
class ModalStep : public AnalysisStep {
public:
    /**
     * @brief Constructor
     * @param name Step name
     */
    explicit ModalStep(const std::string& name = "Modal");
    
    /**
     * @brief Destructor
     */
    ~ModalStep() override = default;
    
    // AnalysisStep interface implementation
    bool execute(FEModel* model) override;
    void initialize(FEModel* model) override;
    void finalize(FEModel* model) override;
    AnalysisType getAnalysisType() const override { return AnalysisType::MODAL; }
    
    // Configuration
    void setNumberOfModes(int n) { numModes_ = n; }
    int getNumberOfModes() const { return numModes_; }
    
    void setFrequencyRange(double minFreq, double maxFreq);
    void getFrequencyRange(double& minFreq, double& maxFreq) const {
        minFreq = minFreq_; maxFreq = maxFreq_;
    }
    
    // Normalization options
    enum class NormalizationType {
        MASS,          // phi^T * M * phi = I
        MAX_VALUE,     // max(|phi|) = 1
        EUCLIDEAN      // ||phi|| = 1
    };
    
    void setNormalizationType(NormalizationType type) { normType_ = type; }
    NormalizationType getNormalizationType() const { return normType_; }
    
    // Eigenvalue solver options
    enum class SolverType {
        SUBSPACE_ITERATION,  // Subspace iteration method
        LANCZOS,             // Lanczos method
        JACOBI_DAVIDSON,     // Jacobi-Davidson method
        SHIFT_INVERT        // Shift-invert method
    };
    
    void setSolverType(SolverType type) { solverType_ = type; }
    SolverType getSolverType() const { return solverType_; }
    
    void setMaxIterations(int maxIter) { maxIterations_ = maxIter; }
    void setTolerance(double tol) { tolerance_ = tol; }
    
    // Results access
    const std::vector<double>& getEigenvalues() const { return eigenvalues_; }
    const std::vector<double>& getFrequencies() const { return frequencies_; }
    const std::vector<double>& getPeriods() const { return periods_; }
    const std::vector<std::vector<double>>& getEigenvectors() const { return eigenvectors_; }
    
    // Modal participation factors
    const std::vector<double>& getParticipationFactors() const { return participationFactors_; }
    const std::vector<double>& getEffectiveMasses() const { return effectiveMasses_; }
    
    // Get specific mode
    void getMode(int modeIndex, double& frequency, std::vector<double>& modeShape) const;
    
private:
    /**
     * @brief Solve the generalized eigenvalue problem
     * @param model FE model
     * @return true if successful
     */
    bool solveEigenProblem(FEModel* model);
    
    /**
     * @brief Solve using subspace iteration method
     * @param model FE model
     * @return true if successful
     */
    bool solveBySubspaceIteration(FEModel* model);
    
    /**
     * @brief Solve using Lanczos method
     * @param model FE model
     * @return true if successful
     */
    bool solveByLanczos(FEModel* model);
    
    /**
     * @brief Assemble mass matrix
     * @param model FE model
     */
    void assembleMassMatrix(FEModel* model);
    
    /**
     * @brief Assemble stiffness matrix
     * @param model FE model
     */
    void assembleStiffnessMatrix(FEModel* model);
    
    /**
     * @brief Normalize eigenvectors
     */
    void normalizeEigenvectors(FEModel* model);
    
    /**
     * @brief Convert eigenvalues to frequencies and periods
     */
    void convertToFrequencies();
    
    /**
     * @brief Calculate modal participation factors
     * @param model FE model
     */
    void calculateParticipationFactors(FEModel* model);
    
    /**
     * @brief Calculate effective modal masses
     * @param model FE model
     */
    void calculateEffectiveMasses(FEModel* model);
    
    /**
     * @brief Sort modes by frequency
     */
    void sortModes();
    
    /**
     * @brief Check orthogonality of modes
     * @param model FE model
     */
    void checkOrthogonality(FEModel* model);
    
    // Configuration
    int numModes_;                          ///< Number of modes to extract
    double minFreq_;                        ///< Minimum frequency (Hz)
    double maxFreq_;                        ///< Maximum frequency (Hz)
    NormalizationType normType_;            ///< Normalization type
    SolverType solverType_;                 ///< Eigenvalue solver type
    int maxIterations_;                     ///< Max iterations for iterative solvers
    double tolerance_;                      ///< Convergence tolerance
    
    // Results
    std::vector<double> eigenvalues_;       ///< Eigenvalues (omega^2)
    std::vector<double> frequencies_;       ///< Natural frequencies (Hz)
    std::vector<double> periods_;           ///< Natural periods (s)
    std::vector<std::vector<double>> eigenvectors_;  ///< Mode shapes
    
    // Modal properties
    std::vector<double> participationFactors_;  ///< Modal participation factors
    std::vector<double> effectiveMasses_;       ///< Effective modal masses
};

#endif // MODALSTEP_H
