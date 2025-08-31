#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

// Forward declarations
class Node;
class Material;

/**
 * @brief Abstract base class for finite elements with geometric nonlinearity
 * Supports infinite strain formulation using appropriate strain measures
 */
class Element {
public:
    // Type definitions
    using NodePtr = std::shared_ptr<Node>;
    using MaterialPtr = std::shared_ptr<Material>;
    using Matrix = MatrixXd;
    using Vector = VectorXd;
    using SparseMatrix = SparseMatrix<double>;

protected:
    // Element properties
    int element_id_;
    std::vector<NodePtr> nodes_;
    MaterialPtr material_;
    int num_nodes_;
    int num_dof_per_node_;
    int total_dof_;
    
    // Geometric properties
    Matrix reference_coords_;      // Reference configuration (X)
    Matrix current_coords_;        // Current configuration (x)
    
    // Kinematic variables
    Matrix deformation_gradient_;  // F = dx/dX
    Matrix right_cauchy_green_;    // C = F^T * F
    Matrix green_lagrange_strain_; // E = 0.5 * (C - I)
    Matrix pk2_stress_;           // Second Piola-Kirchhoff stress
    Matrix cauchy_stress_;        // Cauchy stress
    
    // Numerical integration
    Matrix gauss_points_;
    Vector gauss_weights_;
    int num_gauss_points_;

public:
    // Constructors
    Element(int id, const std::vector<NodePtr>& nodes, MaterialPtr material);
    virtual ~Element() = default;
    
    // Pure virtual functions (must be implemented by derived classes)
    virtual void computeShapeFunctions(const Vector& xi, Vector& N) const = 0;
    virtual void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const = 0;
    virtual void setupGaussIntegration() = 0;
    virtual int getNumNodes() const = 0;
    virtual int getDofPerNode() const = 0;
    
    // Main computation functions
    virtual void updateCurrentConfiguration();
    virtual void computeKinematicQuantities();
    virtual void computeInternalForces(Vector& f_int);
    virtual void computeTangentStiffness(SparseMatrix& K_tan);
    virtual void computeStresses();
    
    // Geometric nonlinearity specific functions
    void computeDeformationGradient(const Matrix& dN_dX, Matrix& F);
    void computeGreenLagrangeStrain(const Matrix& F, Matrix& E);
    void computeRightCauchyGreen(const Matrix& F, Matrix& C);
    void computeCauchyStress(const Matrix& F, const Matrix& S, Matrix& sigma);
    
    // Tangent stiffness components
    void computeMaterialStiffness(const Vector& xi, Matrix& C_material);
    void computeGeometricStiffness(const Vector& xi, Matrix& K_geo);
    void computeInitialDisplacementStiffness(const Vector& xi, Matrix& K_sigma);
    
    // Jacobian and coordinate transformations
    void computeJacobian(const Matrix& dN_dxi, const Matrix& coords, Matrix& J);
    void computeInverseJacobian(const Matrix& J, Matrix& J_inv, double& det_J);
    void computeCartesianDerivatives(const Matrix& dN_dxi, const Matrix& J_inv, Matrix& dN_dX);
    
    // B-matrix computations for different formulations
    void computeBMatrix_Linear(const Matrix& dN_dX, Matrix& B);
    void computeBMatrix_Nonlinear(const Matrix& dN_dX, const Matrix& F, Matrix& B_L, Matrix& B_NL);
    
    // Strain and stress measures for large deformation
    void computeLogarithmicStrain(const Matrix& C, Matrix& E_log);
    void computePrincipalStrains(const Matrix& E, Vector& principal_strains);
    void computeVonMisesStress(const Matrix& stress, double& vm_stress);
    
    // Update functions
    void updateNodalCoordinates();
    void updateElementState();
    bool checkElementQuality();
    
    // Utility functions
    Matrix getDeformationGradient() const { return deformation_gradient_; }
    Matrix getGreenLagrangeStrain() const { return green_lagrange_strain_; }
    Matrix getPK2Stress() const { return pk2_stress_; }
    Matrix getCauchyStress() const { return cauchy_stress_; }
    double getElementVolume() const;
    double getReferenceVolume() const;
    Vector getElementCentroid() const;
    
    // Access functions
    int getId() const { return element_id_; }
    const std::vector<NodePtr>& getNodes() const { return nodes_; }
    MaterialPtr getMaterial() const { return material_; }
    int getTotalDof() const { return total_dof_; }
    
    // Connectivity
    std::vector<int> getGlobalDofIndices() const;
    void assembleElementMatrix(const Matrix& K_elem, SparseMatrix& K_global);
    void assembleElementVector(const Vector& f_elem, Vector& f_global);

protected:
    // Helper functions for large deformation kinematics
    Matrix computePolarDecomposition(const Matrix& F, Matrix& R, Matrix& U);
    void computeStretchTensors(const Matrix& F, Matrix& U, Matrix& V);
    Matrix computeMatrixLogarithm(const Matrix& A);
    Matrix computeMatrixExponential(const Matrix& A);
    
    // Numerical stability functions
    double computeElementSize() const;
    bool isElementInverted() const;
    void applyElementStabilization(Matrix& K_elem);
    
    // Error checking
    void validateKinematicQuantities() const;
    void checkNumericalStability() const;
};

/**
 * @brief Specific implementation for 2D/3D continuum elements
 */
class ContinuumElement : public Element {
private:
    int spatial_dimension_;
    
public:
    ContinuumElement(int id, const std::vector<NodePtr>& nodes, 
                    MaterialPtr material, int spatial_dim);
    
    // Implement pure virtual functions
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
    int getNumNodes() const override { return nodes_.size(); }
    int getDofPerNode() const override { return spatial_dimension_; }
    
    // Specialized functions for continuum mechanics
    void computePlaneStrainModifications();
    void computeAxiSymmetricModifications(const Vector& xi);
};

/**
 * @brief Enhanced element with additional features for severe deformation
 */
class EnhancedElement : public ContinuumElement {
private:
    bool use_enhanced_strain_;
    bool use_f_bar_method_;
    Matrix enhanced_parameters_;
    
public:
    EnhancedElement(int id, const std::vector<NodePtr>& nodes, 
                   MaterialPtr material, int spatial_dim,
                   bool enhanced_strain = false, bool f_bar = false);
    
    // Enhanced strain method
    void computeEnhancedStrainMatrix(const Vector& xi, Matrix& G);
    void updateEnhancedParameters(const Matrix& residual);
    
    // F-bar method for near-incompressible materials
    void applyFBarMethod(Matrix& F, double& det_F_bar);
    
    // Mixed formulation for volumetric locking
    void computeMixedFormulation(Matrix& K_uu, Matrix& K_up, Vector& f_u);
};

#endif // ELEMENT_H