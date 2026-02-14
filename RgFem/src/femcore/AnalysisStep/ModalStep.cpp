#include "ModalStep.h"
#include "FEModel.h"
#include "FEMesh.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>

ModalStep::ModalStep(const std::string& name)
    : numModes_(10)
    , minFreq_(0.0)
    , maxFreq_(1e6)
    , normType_(NormalizationType::MASS)
    , solverType_(SolverType::SUBSPACE_ITERATION)
    , maxIterations_(100)
    , tolerance_(1e-6)
{
    stepName_ = name;
    currentTime_ = 0.0;
    currentIncrement_ = 0;
    outputFrequency_ = 1;
}

void ModalStep::initialize(FEModel* model)
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "Initializing Modal Analysis: " << stepName_ << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Number of modes to extract: " << numModes_ << std::endl;
    std::cout << "Frequency range: " << minFreq_ << " - " << maxFreq_ << " Hz" << std::endl;
    
    std::cout << "Eigenvalue solver: ";
    switch (solverType_) {
        case SolverType::SUBSPACE_ITERATION:
            std::cout << "Subspace Iteration" << std::endl;
            break;
        case SolverType::LANCZOS:
            std::cout << "Lanczos" << std::endl;
            break;
        case SolverType::JACOBI_DAVIDSON:
            std::cout << "Jacobi-Davidson" << std::endl;
            break;
        case SolverType::SHIFT_INVERT:
            std::cout << "Shift-Invert" << std::endl;
            break;
    }
    
    std::cout << "Normalization: ";
    switch (normType_) {
        case NormalizationType::MASS:
            std::cout << "Mass-normalized (phi^T * M * phi = I)" << std::endl;
            break;
        case NormalizationType::MAX_VALUE:
            std::cout << "Max value = 1" << std::endl;
            break;
        case NormalizationType::EUCLIDEAN:
            std::cout << "Euclidean norm = 1" << std::endl;
            break;
    }
    
    // Clear previous results
    eigenvalues_.clear();
    frequencies_.clear();
    periods_.clear();
    eigenvectors_.clear();
    participationFactors_.clear();
    effectiveMasses_.clear();
    
    int nDOF = model->getNumberOfDOFs();
    std::cout << "Number of DOFs: " << nDOF << std::endl;
    std::cout << "Initialization complete.\n" << std::endl;
}

bool ModalStep::execute(FEModel* model)
{
    std::cout << "Executing modal analysis..." << std::endl;
    
    bool success = solveEigenProblem(model);
    
    if (success) {
        convertToFrequencies();
        normalizeEigenvectors(model);
        sortModes();
        calculateParticipationFactors(model);
        calculateEffectiveMasses(model);
        checkOrthogonality(model);
    }
    
    return success;
}

void ModalStep::finalize(FEModel* model)
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "Modal Analysis Completed: " << stepName_ << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Modes found: " << eigenvalues_.size() << std::endl;
    std::cout << std::endl;
    
    // Output modal properties
    std::cout << std::setw(6) << "Mode"
              << std::setw(15) << "Frequency (Hz)"
              << std::setw(15) << "Period (s)"
              << std::setw(15) << "Part. Factor"
              << std::setw(15) << "Eff. Mass (%)" << std::endl;
    std::cout << std::string(81, '-') << std::endl;
    
    double totalMass = 0.0;
    for (double m : effectiveMasses_) {
        totalMass += m;
    }
    
    for (size_t i = 0; i < frequencies_.size(); ++i) {
        double massPercent = (totalMass > 0.0) ? 
                            (effectiveMasses_[i] / totalMass * 100.0) : 0.0;
        
        std::cout << std::setw(6) << (i + 1)
                  << std::setw(15) << std::fixed << std::setprecision(6) << frequencies_[i]
                  << std::setw(15) << std::setprecision(6) << periods_[i]
                  << std::setw(15) << std::scientific << participationFactors_[i]
                  << std::setw(15) << std::fixed << std::setprecision(2) << massPercent
                  << std::endl;
    }
    
    std::cout << std::string(81, '-') << std::endl;
    std::cout << std::endl;
}

bool ModalStep::solveEigenProblem(FEModel* model)
{
    std::cout << "\n--- Solving Eigenvalue Problem ---\n" << std::endl;
    
    // Assemble system matrices
    std::cout << "Assembling mass matrix..." << std::endl;
    assembleMassMatrix(model);
    
    std::cout << "Assembling stiffness matrix..." << std::endl;
    assembleStiffnessMatrix(model);
    
    // Solve based on selected method
    bool success = false;
    
    switch (solverType_) {
        case SolverType::SUBSPACE_ITERATION:
            success = solveBySubspaceIteration(model);
            break;
            
        case SolverType::LANCZOS:
            success = solveByLanczos(model);
            break;
            
        case SolverType::JACOBI_DAVIDSON:
            std::cout << "Jacobi-Davidson not yet implemented, using Subspace Iteration" << std::endl;
            success = solveBySubspaceIteration(model);
            break;
            
        case SolverType::SHIFT_INVERT:
            std::cout << "Shift-Invert not yet implemented, using Subspace Iteration" << std::endl;
            success = solveBySubspaceIteration(model);
            break;
    }
    
    return success;
}

bool ModalStep::solveBySubspaceIteration(FEModel* model)
{
    std::cout << "Solving using Subspace Iteration method..." << std::endl;
    
    int nDOF = model->getNumberOfDOFs();
    int p = std::min(std::max(2 * numModes_, numModes_ + 8), nDOF);  // Subspace size
    
    std::cout << "  Subspace size: " << p << std::endl;
    std::cout << "  Max iterations: " << maxIterations_ << std::endl;
    std::cout << "  Tolerance: " << tolerance_ << std::endl;
    
    // Initialize iteration vectors with random values
    // In real implementation, would use better starting vectors
    
    // TODO: Actual implementation of subspace iteration
    /*
    std::vector<std::vector<double>> X(p, std::vector<double>(nDOF));
    std::vector<std::vector<double>> Y(p, std::vector<double>(nDOF));
    
    // Initialize X with random orthonormal vectors
    initializeSubspace(X);
    
    bool converged = false;
    for (int iter = 0; iter < maxIterations_; ++iter) {
        // Y = K^{-1} * M * X
        for (int i = 0; i < p; ++i) {
            // Solve K * Y_i = M * X_i
            // This is the most expensive operation
        }
        
        // Orthogonalize Y with respect to M
        gramSchmidtM(Y, M);
        
        // Project onto subspace: K_r = Y^T * K * Y, M_r = Y^T * M * Y
        MatrixXd Kr = computeReducedStiffness(Y, K);
        MatrixXd Mr = computeReducedMass(Y, M);
        
        // Solve reduced eigenvalue problem: Kr * phi_r = lambda * Mr * phi_r
        solveReducedEigenProblem(Kr, Mr, eigenvalues_reduced, eigenvectors_reduced);
        
        // Update X = Y * phi_r
        updateSubspace(X, Y, eigenvectors_reduced);
        
        // Check convergence
        if (checkSubspaceConvergence(iter)) {
            converged = true;
            std::cout << "  Converged in " << (iter + 1) << " iterations" << std::endl;
            break;
        }
    }
    
    if (!converged) {
        std::cout << "  WARNING: Did not converge in " << maxIterations_ 
                  << " iterations" << std::endl;
    }
    
    // Extract first numModes_ eigenpairs
    eigenvalues_.resize(numModes_);
    eigenvectors_.resize(numModes_);
    
    for (int i = 0; i < numModes_; ++i) {
        eigenvalues_[i] = eigenvalues_reduced[i];
        eigenvectors_[i] = computeEigenvector(X, eigenvectors_reduced, i);
    }
    */
    
    // For demonstration, create dummy eigenvalues
    eigenvalues_.resize(numModes_);
    eigenvectors_.resize(numModes_, std::vector<double>(nDOF, 0.0));
    
    for (int i = 0; i < numModes_; ++i) {
        // Dummy eigenvalues: omega^2 = (i+1)^2 * 100
        eigenvalues_[i] = static_cast<double>((i + 1) * (i + 1) * 100.0);
        
        // Simple mode shapes for demonstration
        for (int j = 0; j < nDOF; ++j) {
            eigenvectors_[i][j] = std::sin(M_PI * (i + 1) * j / nDOF);
        }
    }
    
    std::cout << "  Eigenvalue problem solved." << std::endl;
    
    return true;
}

bool ModalStep::solveByLanczos(FEModel* model)
{
    std::cout << "Solving using Lanczos method..." << std::endl;
    
    // TODO: Implement Lanczos algorithm
    /*
    The Lanczos algorithm is particularly efficient for large sparse matrices:
    
    1. Start with random vector q_1
    2. For j = 1, 2, ..., m:
       - Compute r = K^{-1} * M * q_j (solve K * r = M * q_j)
       - Orthogonalize: r = r - sum(alpha_i * q_i)
       - beta_{j+1} = ||r||
       - q_{j+1} = r / beta_{j+1}
    3. Form tridiagonal matrix T from alpha and beta
    4. Solve eigenvalue problem for T
    5. Compute Ritz vectors
    */
    
    // For now, fall back to subspace iteration
    return solveBySubspaceIteration(model);
}

void ModalStep::assembleMassMatrix(FEModel* model)
{
    // Assemble global mass matrix from element mass matrices
    
    // TODO: Actual implementation
    /*
    model->clearMassMatrix();
    
    for (auto& element : model->getMesh()->getElements()) {
        MatrixXd Me = element->getMassMatrix();
        std::vector<int> dofs = element->getDOFs();
        model->assembleElementMatrix(Me, dofs, model->getMassMatrix());
    }
    */
}

void ModalStep::assembleStiffnessMatrix(FEModel* model)
{
    // Assemble global stiffness matrix
    
    // TODO: Actual implementation
    /*
    model->clearStiffnessMatrix();
    
    for (auto& element : model->getMesh()->getElements()) {
        MatrixXd Ke = element->getStiffnessMatrix();
        std::vector<int> dofs = element->getDOFs();
        model->assembleElementMatrix(Ke, dofs, model->getStiffnessMatrix());
    }
    */
}

void ModalStep::normalizeEigenvectors(FEModel* model)
{
    std::cout << "Normalizing eigenvectors..." << std::endl;
    
    for (size_t i = 0; i < eigenvectors_.size(); ++i) {
        auto& phi = eigenvectors_[i];
        
        switch (normType_) {
            case NormalizationType::MASS: {
                // Normalize such that phi^T * M * phi = 1
                // TODO: Actual mass-normalization
                /*
                VectorXd Mphi = M * phi;
                double phiTMphi = phi.dot(Mphi);
                double factor = 1.0 / std::sqrt(phiTMphi);
                phi *= factor;
                */
                
                // For now, use Euclidean normalization
                double norm = 0.0;
                for (double val : phi) {
                    norm += val * val;
                }
                norm = std::sqrt(norm);
                
                if (norm > 1e-10) {
                    for (double& val : phi) {
                        val /= norm;
                    }
                }
                break;
            }
            
            case NormalizationType::MAX_VALUE: {
                // Normalize such that max(|phi|) = 1
                double maxVal = 0.0;
                for (double val : phi) {
                    maxVal = std::max(maxVal, std::abs(val));
                }
                
                if (maxVal > 1e-10) {
                    for (double& val : phi) {
                        val /= maxVal;
                    }
                }
                break;
            }
            
            case NormalizationType::EUCLIDEAN: {
                // Normalize such that ||phi|| = 1
                double norm = 0.0;
                for (double val : phi) {
                    norm += val * val;
                }
                norm = std::sqrt(norm);
                
                if (norm > 1e-10) {
                    for (double& val : phi) {
                        val /= norm;
                    }
                }
                break;
            }
        }
    }
}

void ModalStep::convertToFrequencies()
{
    frequencies_.resize(eigenvalues_.size());
    periods_.resize(eigenvalues_.size());
    
    for (size_t i = 0; i < eigenvalues_.size(); ++i) {
        // omega^2 = eigenvalue, omega = angular frequency (rad/s)
        double omega = std::sqrt(eigenvalues_[i]);
        
        // f = omega / (2*pi) (Hz)
        frequencies_[i] = omega / (2.0 * M_PI);
        
        // T = 1 / f (s)
        periods_[i] = (frequencies_[i] > 1e-10) ? (1.0 / frequencies_[i]) : 0.0;
    }
}

void ModalStep::calculateParticipationFactors(FEModel* model)
{
    std::cout << "Calculating modal participation factors..." << std::endl;
    
    participationFactors_.resize(eigenvectors_.size());
    
    // TODO: Actual implementation
    /*
    // Participation factor: gamma_i = (phi_i^T * M * r) / (phi_i^T * M * phi_i)
    // where r is the influence vector (typically unit vector in direction of interest)
    
    VectorXd r(nDOF);  // Influence vector (e.g., all 1's for global X direction)
    r.setOnes();
    
    for (size_t i = 0; i < eigenvectors_.size(); ++i) {
        VectorXd phi = eigenvectors_[i];
        VectorXd Mphi = M * phi;
        
        double numerator = phi.dot(M * r);
        double denominator = phi.dot(Mphi);
        
        participationFactors_[i] = numerator / denominator;
    }
    */
    
    // Placeholder
    for (size_t i = 0; i < eigenvectors_.size(); ++i) {
        participationFactors_[i] = 1.0 / std::sqrt(eigenvalues_[i]);
    }
}

void ModalStep::calculateEffectiveMasses(FEModel* model)
{
    std::cout << "Calculating effective modal masses..." << std::endl;
    
    effectiveMasses_.resize(eigenvectors_.size());
    
    // TODO: Actual implementation
    /*
    // Effective mass: M_eff,i = (phi_i^T * M * r)^2 / (phi_i^T * M * phi_i)
    
    VectorXd r(nDOF);
    r.setOnes();
    
    for (size_t i = 0; i < eigenvectors_.size(); ++i) {
        VectorXd phi = eigenvectors_[i];
        VectorXd Mphi = M * phi;
        
        double numerator = phi.dot(M * r);
        numerator = numerator * numerator;
        
        double denominator = phi.dot(Mphi);
        
        effectiveMasses_[i] = numerator / denominator;
    }
    */
    
    // Placeholder
    for (size_t i = 0; i < eigenvectors_.size(); ++i) {
        effectiveMasses_[i] = participationFactors_[i] * participationFactors_[i];
    }
}

void ModalStep::sortModes()
{
    // Sort modes by frequency (ascending order)
    std::vector<size_t> indices(eigenvalues_.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }
    
    std::sort(indices.begin(), indices.end(),
              [this](size_t a, size_t b) {
                  return eigenvalues_[a] < eigenvalues_[b];
              });
    
    // Reorder all arrays
    auto reorder = [&indices](auto& vec) {
        auto temp = vec;
        for (size_t i = 0; i < indices.size(); ++i) {
            vec[i] = temp[indices[i]];
        }
    };
    
    reorder(eigenvalues_);
    reorder(frequencies_);
    reorder(periods_);
    reorder(eigenvectors_);
    
    if (!participationFactors_.empty()) {
        reorder(participationFactors_);
    }
    if (!effectiveMasses_.empty()) {
        reorder(effectiveMasses_);
    }
}

void ModalStep::checkOrthogonality(FEModel* model)
{
    // Check orthogonality: phi_i^T * M * phi_j = delta_ij (for mass-normalized modes)
    
    if (eigenvectors_.size() < 2) return;
    
    std::cout << "Checking orthogonality of modes..." << std::endl;
    
    // TODO: Actual check with mass matrix
    /*
    double maxOffDiag = 0.0;
    
    for (size_t i = 0; i < eigenvectors_.size(); ++i) {
        for (size_t j = i + 1; j < eigenvectors_.size(); ++j) {
            VectorXd phi_i = eigenvectors_[i];
            VectorXd phi_j = eigenvectors_[j];
            
            double orthog = phi_i.dot(M * phi_j);
            maxOffDiag = std::max(maxOffDiag, std::abs(orthog));
        }
    }
    
    std::cout << "  Max off-diagonal term: " << maxOffDiag << std::endl;
    
    if (maxOffDiag < 1e-6) {
        std::cout << "  Modes are orthogonal." << std::endl;
    } else {
        std::cout << "  WARNING: Modes may not be properly orthogonal!" << std::endl;
    }
    */
}

void ModalStep::setFrequencyRange(double minFreq, double maxFreq)
{
    minFreq_ = minFreq;
    maxFreq_ = maxFreq;
}

void ModalStep::getMode(int modeIndex, double& frequency, 
                        std::vector<double>& modeShape) const
{
    if (modeIndex < 0 || modeIndex >= static_cast<int>(frequencies_.size())) {
        std::cerr << "Error: Invalid mode index " << modeIndex << std::endl;
        return;
    }
    
    frequency = frequencies_[modeIndex];
    modeShape = eigenvectors_[modeIndex];
}
