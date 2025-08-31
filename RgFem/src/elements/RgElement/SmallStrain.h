#ifndef SMALL_STRAIN_ELEMENT_H
#define SMALL_STRAIN_ELEMENT_H

#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

// Forward declarations
class Node;
class Material;

/**
 * @brief Abstract base class for small strain finite elements
 * Uses linear strain-displacement relations: ε = ½(∇u + ∇u^T)
 * Assumes small displacements and small rotations
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
    int spatial_dimension_;
    
    // Current element state
    Vector nodal_displacements_;    // Current nodal displacements
    Matrix element_strain_;         // Current strain at integration points
    Matrix element_stress_;         // Current stress at integration points
    
    // Element matrices (computed once for small strain)
    Matrix stiffness_matrix_;       // Element stiffness matrix [K_e]
    Matrix mass_matrix_;           // Element mass matrix [M_e] 
    Vector internal_force_;        // Internal force vector {f_int}
    
    // Constitutive relationship
    Matrix constitutive_matrix_;   // Material stiffness matrix [D]
    
    // Numerical integration
    Matrix gauss_points_;          // Gauss integration points in natural coordinates
    Vector gauss_weights_;         // Gauss integration weights
    int num_gauss_points_;
    
    // Geometric properties
    double element_volume_;        // Element volume/area
    double element_thickness_;     // Element thickness (for 2D elements)
    
    // Flags for optimization
    bool stiffness_computed_;      // Whether stiffness matrix is up-to-date
    bool mass_computed_;          // Whether mass matrix is up-to-date

public:
    // Constructors
    Element(int id, const std::vector<NodePtr>& nodes, MaterialPtr material, int spatial_dim);
    virtual ~Element() = default;
    
    // Pure virtual functions (must be implemented by derived classes)
    virtual void computeShapeFunctions(const Vector& xi, Vector& N) const = 0;
    virtual void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const = 0;
    virtual void setupGaussIntegration() = 0;
    virtual int getNumNodes() const = 0;
    virtual int getDofPerNode() const = 0;
    virtual std::string getElementType() const = 0;
    
    // Main computation functions
    virtual void computeStiffnessMatrix();
    virtual void computeMassMatrix();
    virtual void computeInternalForces();
    virtual void computeStrainStress();
    virtual void updateElementState(const Vector& global_displacement);
    
    // Small strain specific functions
    void computeBMatrix(const Vector& xi, Matrix& B);
    void computeConstitutiveMatrix(Matrix& D);
    void computeStrain(const Matrix& B, const Vector& u_elem, Vector& strain);
    void computeStress(const Vector& strain, Vector& stress);
    
    // Coordinate transformation and Jacobian
    void computeJacobian(const Matrix& dN_dxi, Matrix& J, double& det_J);
    void computeCartesianDerivatives(const Matrix& dN_dxi, const Matrix& J_inv, Matrix& dN_dx);
    void mapToGlobalCoordinates(const Vector& xi, Vector& global_coords);
    
    // Integration and assembly
    double integrateOverElement(const std::function<double(const Vector&)>& integrand);
    void assembleElementMatrix(const Matrix& elem_matrix, SparseMatrix& global_matrix);
    void assembleElementVector(const Vector& elem_vector, Vector& global_vector);
    
    // Element quality and validation
    bool validateElementGeometry();
    double computeAspectRatio();
    double computeSkewness();
    double computeElementSize();
    
    // Stress recovery and smoothing
    void extrapolateStressToNodes(Matrix& nodal_stress);
    void computePrincipalStresses(const Vector& stress, Vector& principal_stresses, Matrix& directions);
    double computeVonMisesStress(const Vector& stress);
    double computeMaxShearStress(const Vector& stress);
    
    // Material property queries
    double getYoungsModulus() const;
    double getPoissonsRatio() const;
    double getDensity() const;
    double getThickness() const { return element_thickness_; }
    void setThickness(double thickness) { element_thickness_ = thickness; }
    
    // Access functions
    int getId() const { return element_id_; }
    const std::vector<NodePtr>& getNodes() const { return nodes_; }
    MaterialPtr getMaterial() const { return material_; }
    int getTotalDof() const { return total_dof_; }
    int getSpatialDimension() const { return spatial_dimension_; }
    
    Matrix getStiffnessMatrix() const { return stiffness_matrix_; }
    Matrix getMassMatrix() const { return mass_matrix_; }
    Vector getInternalForces() const { return internal_force_; }
    Matrix getCurrentStrain() const { return element_strain_; }
    Matrix getCurrentStress() const { return element_stress_; }
    double getVolume() const { return element_volume_; }
    
    // Connectivity and DOF mapping
    std::vector<int> getGlobalDofIndices() const;
    Vector getElementDisplacements(const Vector& global_displacement) const;
    void scatterGlobalToElement(const Vector& global_vector, Vector& element_vector) const;
    void gatherElementToGlobal(const Vector& element_vector, Vector& global_vector) const;
    
    // Output functions
    void printElementInfo() const;
    void exportElementData(const std::string& filename) const;

protected:
    // Helper functions for small strain formulation
    void computeElementVolume();
    void initializeElementMatrices();
    void updateInternalState();
    
    // Numerical integration helpers
    Vector transformGaussPoint(const Vector& xi) const;
    double getGaussWeight(int point_index) const;
    
    // B-matrix computation for different element types
    virtual void computeLinearBMatrix(const Matrix& dN_dx, Matrix& B) = 0;
    
    // Validation helpers
    bool checkStrainMagnitude() const;  // Verify small strain assumption
    bool checkRotationMagnitude() const; // Verify small rotation assumption
    void warnIfLargeDeformation() const;
};

/**
 * @brief 2D Continuum Element for small strain analysis
 */
class ContinuumElement2D : public Element {
protected:
    enum class PlaneType {
        PLANE_STRESS,    // σ_z = 0 (thin plates)
        PLANE_STRAIN,    // ε_z = 0 (long cylinders) 
        AXISYMMETRIC     // Axisymmetric problems
    };
    
    PlaneType plane_type_;

public:
    ContinuumElement2D(int id, const std::vector<NodePtr>& nodes, MaterialPtr material, PlaneType type);
    
    int getDofPerNode() const override { return 2; } // u, v displacements
    
    // 2D specific functions
    void setPlaneType(PlaneType type) { plane_type_ = type; }
    PlaneType getPlaneType() const { return plane_type_; }
    
protected:
    void computeLinearBMatrix(const Matrix& dN_dx, Matrix& B) override;
    void compute2DConstitutiveMatrix(Matrix& D);
};

/**
 * @brief 3D Continuum Element for small strain analysis
 */
class ContinuumElement3D : public Element {
public:
    ContinuumElement3D(int id, const std::vector<NodePtr>& nodes, MaterialPtr material);
    
    int getDofPerNode() const override { return 3; } // u, v, w displacements
    
protected:
    void computeLinearBMatrix(const Matrix& dN_dx, Matrix& B) override;
    void compute3DConstitutiveMatrix(Matrix& D);
};

/**
 * @brief Specific 2D element implementations
 */

// 3-node triangular element (CST - Constant Strain Triangle)
class Triangle3 : public ContinuumElement2D {
public:
    Triangle3(int id, const std::vector<NodePtr>& nodes, MaterialPtr material, PlaneType type);
    
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
    int getNumNodes() const override { return 3; }
    std::string getElementType() const override { return "Triangle3"; }
};

// 4-node quadrilateral element (Q4)
class Quadrilateral4 : public ContinuumElement2D {
public:
    Quadrilateral4(int id, const std::vector<NodePtr>& nodes, MaterialPtr material, PlaneType type);
    
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
    int getNumNodes() const override { return 4; }
    std::string getElementType() const override { return "Quadrilateral4"; }
};

// 6-node triangular element (LST - Linear Strain Triangle)
class Triangle6 : public ContinuumElement2D {
public:
    Triangle6(int id, const std::vector<NodePtr>& nodes, MaterialPtr material, PlaneType type);
    
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
    int getNumNodes() const override { return 6; }
    std::string getElementType() const override { return "Triangle6"; }
};

// 8-node quadrilateral element (Q8)
class Quadrilateral8 : public ContinuumElement2D {
public:
    Quadrilateral8(int id, const std::vector<NodePtr>& nodes, MaterialPtr material, PlaneType type);
    
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
    int getNumNodes() const override { return 8; }
    std::string getElementType() const override { return "Quadrilateral8"; }
};

/**
 * @brief Specific 3D element implementations
 */

// 4-node tetrahedral element
class Tetrahedron4 : public ContinuumElement3D {
public:
    Tetrahedron4(int id, const std::vector<NodePtr>& nodes, MaterialPtr material);
    
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
    int getNumNodes() const override { return 4; }
    std::string getElementType() const override { return "Tetrahedron4"; }
};

// 8-node hexahedral element (brick element)
class Hexahedron8 : public ContinuumElement3D {
public:
    Hexahedron8(int id, const std::vector<NodePtr>& nodes, MaterialPtr material);
    
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
    int getNumNodes() const override { return 8; }
    std::string getElementType() const override { return "Hexahedron8"; }
};

// 10-node tetrahedral element
class Tetrahedron10 : public ContinuumElement3D {
public:
    Tetrahedron10(int id, const std::vector<NodePtr>& nodes, MaterialPtr material);
    
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
    int getNumNodes() const override { return 10; }
    std::string getElementType() const override { return "Tetrahedron10"; }
};

// 20-node hexahedral element
class Hexahedron20 : public ContinuumElement3D {
public:
    Hexahedron20(int id, const std::vector<NodePtr>& nodes, MaterialPtr material);
    
    void computeShapeFunctions(const Vector& xi, Vector& N) const override;
    void computeShapeDerivatives(const Vector& xi, Matrix& dN_dxi) const override;
    void setupGaussIntegration() override;
    int getNumNodes() const override { return 20; }
    std::string getElementType() const override { return "Hexahedron20"; }
};

/**
 * @brief Element factory for creating small strain elements
 */
class ElementFactory {
public:
    enum class ElementType {
        TRIANGLE3,
        TRIANGLE6, 
        QUADRILATERAL4,
        QUADRILATERAL8,
        TETRAHEDRON4,
        TETRAHEDRON10,
        HEXAHEDRON8,
        HEXAHEDRON20
    };
    
    static std::shared_ptr<Element> createElement(
        ElementType type,
        int id,
        const std::vector<std::shared_ptr<Node>>& nodes,
        std::shared_ptr<Material> material,
        ContinuumElement2D::PlaneType plane_type = ContinuumElement2D::PlaneType::PLANE_STRESS);
        
    static std::string getElementTypeName(ElementType type);
};

#endif // SMALL_STRAIN_ELEMENT_H