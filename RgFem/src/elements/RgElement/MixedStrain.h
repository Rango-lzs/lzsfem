#ifndef ELEMENT_HIERARCHY_H
#define ELEMENT_HIERARCHY_H

#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

// Forward declarations
class Node;
class Material;

/**
 * @brief Abstract base class for all finite elements
 * Contains common functionality regardless of strain formulation
 */
class ElementBase {
public:
    using NodePtr = std::shared_ptr<Node>;
    using MaterialPtr = std::shared_ptr<Material>;
    using Matrix = MatrixXd;
    using Vector = VectorXd;
    using SparseMatrix = SparseMatrix<double>;

protected:
    // Common properties
    int element_id_;
    std::vector<NodePtr> nodes_;
    MaterialPtr material_;
    int num_nodes_;
    int num_dof_per_node_;
    int total_dof_;
    int spatial_dimension_;
    
    // Integration
    Matrix gauss_points_;
    Vector gauss_weights_;
    int num_gauss_points_;

public:
    ElementBase(int id, const std::vector<NodePtr>& nodes, MaterialPtr material, int spatial_dim);
    virtual ~ElementBase() = default;
    
    // Pure virtual functions - must be implemented by all elements
    virtual void computeInternalForces(Vector& f_int) = 0;
    virtual void computeTangentStiffness(SparseMatrix& K_tan) = 0;
    virtual void updateElementState() = 0;
    
    // Common virtual functions with default implementations
    virtual void computeShapeFunctions(const Vector& xi, Vector& N) const = 0;
    virtual void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const = 0;
    virtual void setupGaussIntegration() = 0;
    
    // Common utility functions
    void computeJacobian(const Matrix& dN_dxi, const Matrix& coords, Matrix& J);
    void computeInverseJacobian(const Matrix& J, Matrix& J_inv, double& det_J);
    std::vector<int> getGlobalDofIndices() const;
    void assembleElementMatrix(const Matrix& K_elem, SparseMatrix& K_global);
    void assembleElementVector(const Vector& f_elem, Vector& f_global);
    
    // Access functions
    int getId() const { return element_id_; }
    const std::vector<NodePtr>& getNodes() const { return nodes_; }
    MaterialPtr getMaterial() const { return material_; }
    int getTotalDof() const { return total_dof_; }
    int getSpatialDimension() const { return spatial_dimension_; }
};

/**
 * @brief Element class for small strain analysis
 * Uses linear strain-displacement relations and engineering strain
 */
class SmallStrainElement : public ElementBase {
protected:
    // Small strain specific data
    Matrix current_strain_;        // Engineering strain tensor
    Matrix current_stress_;        // Stress tensor
    Matrix constitutive_matrix_;   // Material stiffness matrix
    
    // Kinematic matrices
    Matrix B_matrix_;              // Strain-displacement matrix (constant)
    bool B_matrix_computed_;

public:
    SmallStrainElement(int id, const std::vector<NodePtr>& nodes, 
                      MaterialPtr material, int spatial_dim);
    
    // Implement pure virtual functions
    void computeInternalForces(Vector& f_int) override;
    void computeTangentStiffness(SparseMatrix& K_tan) override;
    void updateElementState() override;
    
    // Small strain specific functions
    void computeBMatrix(const Vector& xi, Matrix& B);
    void computeStrain(const Vector& u_elem, Matrix& strain);
    void computeStress(const Matrix& strain, Matrix& stress);
    void computeConstitutiveMatrix(Matrix& D);
    
    // Simple updates (no geometric nonlinearity)
    void updateStrain();
    void updateStress();
    
    // Access functions
    Matrix getCurrentStrain() const { return current_strain_; }
    Matrix getCurrentStress() const { return current_stress_; }
    double getVonMisesStress() const;
    
protected:
    // Helper functions for small strain
    void computeLinearBMatrix(const Matrix& dN_dx, Matrix& B);
    bool validateSmallStrainAssumption() const;
};

/**
 * @brief Element class for geometrically nonlinear analysis with finite strains
 * Uses nonlinear strain-displacement relations and appropriate strain measures
 */
class LargeStrainElement : public ElementBase {
protected:
    // Large strain specific data
    Matrix reference_coords_;      // Reference configuration
    Matrix current_coords_;        // Current/deformed configuration
    Matrix deformation_gradient_;  // F = dx/dX
    Matrix right_cauchy_green_;    // C = F^T * F
    Matrix green_lagrange_strain_; // E = 0.5 * (C - I)
    Matrix pk2_stress_;           // Second Piola-Kirchhoff stress
    Matrix cauchy_stress_;        // Cauchy stress
    
    // History variables for path-dependent materials
    std::vector<Matrix> stress_history_;
    std::vector<Matrix> strain_history_;

public:
    LargeStrainElement(int id, const std::vector<NodePtr>& nodes, 
                      MaterialPtr material, int spatial_dim);
    
    // Implement pure virtual functions
    void computeInternalForces(Vector& f_int) override;
    void computeTangentStiffness(SparseMatrix& K_tan) override;
    void updateElementState() override;
    
    // Large strain specific functions
    void updateCurrentConfiguration();
    void computeKinematicQuantities();
    void computeDeformationGradient(const Matrix& dN_dX, Matrix& F);
    void computeGreenLagrangeStrain(const Matrix& F, Matrix& E);
    void computeStresses();
    
    // Tangent stiffness components
    void computeMaterialStiffness(const Vector& xi, Matrix& K_mat);
    void computeGeometricStiffness(const Vector& xi, Matrix& K_geo);
    void computeInitialDisplacementStiffness(const Vector& xi, Matrix& K_sigma);
    
    // B-matrix for nonlinear analysis
    void computeNonlinearBMatrix(const Matrix& dN_dX, const Matrix& F, 
                                Matrix& B_L, Matrix& B_NL);
    
    // Strain measures for large deformation
    void computeLogarithmicStrain(const Matrix& C, Matrix& E_log);
    void computeAlmansiStrain(const Matrix& F, Matrix& e);
    
    // Stress transformations
    void computeCauchyStress(const Matrix& F, const Matrix& S, Matrix& sigma);
    void computeKirchhoffStress(const Matrix& F, const Matrix& S, Matrix& tau);
    
    // Stability and quality checks
    bool checkElementInversion() const;
    bool checkElementQuality() const;
    void applyStabilization(Matrix& K_elem);
    
    // Access functions
    Matrix getDeformationGradient() const { return deformation_gradient_; }
    Matrix getGreenLagrangeStrain() const { return green_lagrange_strain_; }
    Matrix getPK2Stress() const { return pk2_stress_; }
    Matrix getCauchyStress() const { return cauchy_stress_; }
    double getElementVolume() const;
    double getReferenceVolume() const;
    
protected:
    // Helper functions for large deformation
    void computeCartesianDerivatives(const Matrix& dN_dxi, const Matrix& J_inv, 
                                   Matrix& dN_dX, Matrix& dN_dx);
    Matrix computePolarDecomposition(const Matrix& F, Matrix& R, Matrix& U);
    double computeJacobian(const Matrix& F) const;
    void updateHistoryVariables();
};

/**
 * @brief Factory class to create appropriate element types
 */
class ElementFactory {
public:
    enum class AnalysisType {
        SMALL_STRAIN,
        LARGE_STRAIN,
        AUTO_DETECT  // Automatically switch based on strain magnitude
    };
    
    static std::shared_ptr<ElementBase> createElement(
        int id, 
        const std::vector<std::shared_ptr<Node>>& nodes,
        std::shared_ptr<Material> material,
        int spatial_dim,
        AnalysisType analysis_type = AnalysisType::SMALL_STRAIN);
};

/**
 * @brief Specific implementations for different element geometries
 */

// 2D Triangular elements
class Triangle3SmallStrain : public SmallStrainElement {
public:
    Triangle3SmallStrain(int id, const std::vector<NodePtr>& nodes, MaterialPtr material);
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
};

class Triangle3LargeStrain : public LargeStrainElement {
public:
    Triangle3LargeStrain(int id, const std::vector<NodePtr>& nodes, MaterialPtr material);
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
};

// 2D Quadrilateral elements  
class Quad4SmallStrain : public SmallStrainElement {
public:
    Quad4SmallStrain(int id, const std::vector<NodePtr>& nodes, MaterialPtr material);
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
};

class Quad4LargeStrain : public LargeStrainElement {
public:
    Quad4LargeStrain(int id, const std::vector<NodePtr>& nodes, MaterialPtr material);
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
};

// 3D Tetrahedral elements
class Tetrahedron4SmallStrain : public SmallStrainElement {
public:
    Tetrahedron4SmallStrain(int id, const std::vector<NodePtr>& nodes, MaterialPtr material);
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
};

class Tetrahedron4LargeStrain : public LargeStrainElement {
public:
    Tetrahedron4LargeStrain(int id, const std::vector<NodePtr>& nodes, MaterialPtr material);
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
};

// 3D Hexahedral elements
class Hexahedron8SmallStrain : public SmallStrainElement {
public:
    Hexahedron8SmallStrain(int id, const std::vector<NodePtr>& nodes, MaterialPtr material);
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
};

class Hexahedron8LargeStrain : public LargeStrainElement {
public:
    Hexahedron8LargeStrain(int id, const std::vector<NodePtr>& nodes, MaterialPtr material);
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
};

#endif // ELEMENT_HIERARCHY_H