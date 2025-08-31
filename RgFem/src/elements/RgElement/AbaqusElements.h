/**
 * Element Class Hierarchy with Infinite Strain Capability
 * Based on ABAQUS Element Library with Geometric Nonlinearity Extensions
 */

// ============================================================================
// ABSTRACT BASE CLASSES WITH STRAIN FORMULATION DISTINCTION
// ============================================================================

class Element {
    // Root base class for all elements
};

// Strain formulation base classes
class SmallStrainElementBase : public Element {
    // Base for small strain formulation (linear geometry)
    // Uses engineering strain: ε = ½(∇u + ∇u^T)
};

class LargeStrainElementBase : public Element {
    // Base for large strain formulation (nonlinear geometry)
    // Uses Green-Lagrange strain: E = ½(F^T F - I)
    // Supports infinite strains
};

// ============================================================================
// CONTINUUM ELEMENT BASE CLASSES
// ============================================================================

// Small strain continuum elements
class SmallStrainContinuumElement : public SmallStrainElementBase {
    // Base for small strain solid/continuum elements
};

class SmallStrainCPE_Element : public SmallStrainContinuumElement {};  // 2D Plane Strain
class SmallStrainCPS_Element : public SmallStrainContinuumElement {};  // 2D Plane Stress  
class SmallStrainCAX_Element : public SmallStrainContinuumElement {};  // 2D Axisymmetric
class SmallStrainC3D_Element : public SmallStrainContinuumElement {};  // 3D Stress

// Large strain continuum elements (infinite strain capable)
class LargeStrainContinuumElement : public LargeStrainElementBase {
    // Base for large strain solid/continuum elements
    // Handles finite deformation kinematics
};

class LargeStrainCPE_Element : public LargeStrainContinuumElement {};   // 2D Plane Strain
class LargeStrainCPS_Element : public LargeStrainContinuumElement {};   // 2D Plane Stress  
class LargeStrainCAX_Element : public LargeStrainContinuumElement {};   // 2D Axisymmetric
class LargeStrainC3D_Element : public LargeStrainContinuumElement {};   // 3D Stress

// ============================================================================
// STRUCTURAL ELEMENT BASE CLASSES
// ============================================================================

// Small strain structural elements
class SmallStrainStructuralElement : public SmallStrainElementBase {};

class SmallStrainBeamElement : public SmallStrainStructuralElement {};
class SmallStrainShellElement : public SmallStrainStructuralElement {};
class SmallStrainMembraneElement : public SmallStrainStructuralElement {};
class SmallStrainTrussElement : public SmallStrainStructuralElement {};

// Large strain structural elements (for large rotation/deformation)
class LargeStrainStructuralElement : public LargeStrainElementBase {};

class LargeStrainBeamElement : public LargeStrainStructuralElement {};   // Corotational beams
class LargeStrainShellElement : public LargeStrainStructuralElement {};  // Finite rotation shells
class LargeStrainMembraneElement : public LargeStrainStructuralElement {};
class LargeStrainTrussElement : public LargeStrainStructuralElement {};

// ============================================================================
// SPECIFIC ELEMENT IMPLEMENTATIONS - SMALL STRAIN
// ============================================================================

// ---- 2D Small Strain Elements ----
class CPE3 : public SmallStrainCPE_Element {};    // 3-node triangular plane strain
class CPE4 : public SmallStrainCPE_Element {};    // 4-node quad plane strain
class CPE4R : public SmallStrainCPE_Element {};   // 4-node reduced integration
class CPE4I : public SmallStrainCPE_Element {};   // 4-node incompatible mode
class CPE6 : public SmallStrainCPE_Element {};    // 6-node triangular
class CPE8 : public SmallStrainCPE_Element {};    // 8-node quad
class CPE8R : public SmallStrainCPE_Element {};   // 8-node reduced integration

class CPS3 : public SmallStrainCPS_Element {};    // 3-node triangular plane stress
class CPS4 : public SmallStrainCPS_Element {};    // 4-node quad plane stress
class CPS4R : public SmallStrainCPS_Element {};   // 4-node reduced integration
class CPS4I : public SmallStrainCPS_Element {};   // 4-node incompatible mode
class CPS6 : public SmallStrainCPS_Element {};    // 6-node triangular
class CPS8 : public SmallStrainCPS_Element {};    // 8-node quad
class CPS8R : public SmallStrainCPS_Element {};   // 8-node reduced integration

class CAX3 : public SmallStrainCAX_Element {};    // 3-node triangular axisymmetric
class CAX4 : public SmallStrainCAX_Element {};    // 4-node quad axisymmetric
class CAX4R : public SmallStrainCAX_Element {};   // 4-node reduced integration
class CAX4I : public SmallStrainCAX_Element {};   // 4-node incompatible mode
class CAX6 : public SmallStrainCAX_Element {};    // 6-node triangular
class CAX8 : public SmallStrainCAX_Element {};    // 8-node quad
class CAX8R : public SmallStrainCAX_Element {};   // 8-node reduced integration

// ---- 3D Small Strain Elements ----
class C3D4 : public SmallStrainC3D_Element {};    // 4-node tetrahedral
class C3D6 : public SmallStrainC3D_Element {};    // 6-node triangular prism
class C3D8 : public SmallStrainC3D_Element {};    // 8-node hexahedral
class C3D8R : public SmallStrainC3D_Element {};   // 8-node reduced integration
class C3D8I : public SmallStrainC3D_Element {};   // 8-node incompatible mode
class C3D10 : public SmallStrainC3D_Element {};   // 10-node tetrahedral
class C3D15 : public SmallStrainC3D_Element {};   // 15-node triangular prism
class C3D20 : public SmallStrainC3D_Element {};   // 20-node hexahedral
class C3D20R : public SmallStrainC3D_Element {};  // 20-node reduced integration
class C3D27 : public SmallStrainC3D_Element {};   // 27-node hexahedral

// ============================================================================
// SPECIFIC ELEMENT IMPLEMENTATIONS - LARGE STRAIN (INFINITE STRAIN CAPABLE)
// ============================================================================

// ---- 2D Large Strain Elements ----
// Using "NL" (Nonlinear) suffix to distinguish from small strain versions
class CPE3NL : public LargeStrainCPE_Element {};    // 3-node triangular plane strain - large strain
class CPE4NL : public LargeStrainCPE_Element {};    // 4-node quad plane strain - large strain
class CPE4RNL : public LargeStrainCPE_Element {};   // 4-node reduced integration - large strain
class CPE4INL : public LargeStrainCPE_Element {};   // 4-node incompatible mode - large strain
class CPE6NL : public LargeStrainCPE_Element {};    // 6-node triangular - large strain
class CPE8NL : public LargeStrainCPE_Element {};    // 8-node quad - large strain
class CPE8RNL : public LargeStrainCPE_Element {};   // 8-node reduced integration - large strain

class CPS3NL : public LargeStrainCPS_Element {};    // 3-node triangular plane stress - large strain
class CPS4NL : public LargeStrainCPS_Element {};    // 4-node quad plane stress - large strain
class CPS4RNL : public LargeStrainCPS_Element {};   // 4-node reduced integration - large strain
class CPS4INL : public LargeStrainCPS_Element {};   // 4-node incompatible mode - large strain
class CPS6NL : public LargeStrainCPS_Element {};    // 6-node triangular - large strain
class CPS8NL : public LargeStrainCPS_Element {};    // 8-node quad - large strain
class CPS8RNL : public LargeStrainCPS_Element {};   // 8-node reduced integration - large strain

class CAX3NL : public LargeStrainCAX_Element {};    // 3-node triangular axisymmetric - large strain
class CAX4NL : public LargeStrainCAX_Element {};    // 4-node quad axisymmetric - large strain
class CAX4RNL : public LargeStrainCAX_Element {};   // 4-node reduced integration - large strain
class CAX4INL : public LargeStrainCAX_Element {};   // 4-node incompatible mode - large strain
class CAX6NL : public LargeStrainCAX_Element {};    // 6-node triangular - large strain
class CAX8NL : public LargeStrainCAX_Element {};    // 8-node quad - large strain
class CAX8RNL : public LargeStrainCAX_Element {};   // 8-node reduced integration - large strain

// ---- 3D Large Strain Elements ----
class C3D4NL : public LargeStrainC3D_Element {};    // 4-node tetrahedral - large strain
class C3D6NL : public LargeStrainC3D_Element {};    // 6-node triangular prism - large strain
class C3D8NL : public LargeStrainC3D_Element {};    // 8-node hexahedral - large strain
class C3D8RNL : public LargeStrainC3D_Element {};   // 8-node reduced integration - large strain
class C3D8INL : public LargeStrainC3D_Element {};   // 8-node incompatible mode - large strain
class C3D10NL : public LargeStrainC3D_Element {};   // 10-node tetrahedral - large strain
class C3D15NL : public LargeStrainC3D_Element {};   // 15-node triangular prism - large strain
class C3D20NL : public LargeStrainC3D_Element {};   // 20-node hexahedral - large strain
class C3D20RNL : public LargeStrainC3D_Element {};  // 20-node reduced integration - large strain
class C3D27NL : public LargeStrainC3D_Element {};   // 27-node hexahedral - large strain

// ============================================================================
// ENHANCED ELEMENTS FOR SEVERE DEFORMATION
// ============================================================================

// Enhanced strain elements for better performance under large deformation
class EnhancedStrainElementBase : public LargeStrainElementBase {
    // Base for enhanced strain formulation elements
    // Uses additional internal DOFs to improve element behavior
};

class C3D8ENL : public EnhancedStrainElementBase {};   // 8-node enhanced strain hexahedral
class C3D20ENL : public EnhancedStrainElementBase {};  // 20-node enhanced strain hexahedral
class CPE4ENL : public EnhancedStrainElementBase {};   // 4-node enhanced strain plane strain
class CPS4ENL : public EnhancedStrainElementBase {};   // 4-node enhanced strain plane stress

// F-bar elements for near-incompressible materials
class FBarElementBase : public LargeStrainElementBase {
    // Base for F-bar method elements
    // Prevents volumetric locking in nearly incompressible materials
};

class C3D8FNL : public FBarElementBase {};    // 8-node F-bar hexahedral
class C3D4FNL : public FBarElementBase {};    // 4-node F-bar tetrahedral

// ============================================================================
// STRUCTURAL ELEMENTS - LARGE STRAIN
// ============================================================================

// ---- Large Strain Beam Elements ----
class B21NL : public LargeStrainBeamElement {};     // 2-node corotational beam 2D
class B22NL : public LargeStrainBeamElement {};     // 3-node corotational beam 2D
class B31NL : public LargeStrainBeamElement {};     // 2-node corotational beam 3D
class B32NL : public LargeStrainBeamElement {};     // 3-node corotational beam 3D

// ---- Large Strain Shell Elements ----
class S3NL : public LargeStrainShellElement {};     // 3-node finite rotation shell
class S4NL : public LargeStrainShellElement {};     // 4-node finite rotation shell
class S4RNL : public LargeStrainShellElement {};    // 4-node reduced integration finite rotation shell
class S8RNL : public LargeStrainShellElement {};    // 8-node finite rotation shell

// ---- Large Strain Membrane Elements ----
class M3D3NL : public LargeStrainMembraneElement {};  // 3-node triangular membrane - large strain
class M3D4NL : public LargeStrainMembraneElement {};  // 4-node quad membrane - large strain
class M3D4RNL : public LargeStrainMembraneElement {}; // 4-node reduced integration membrane - large strain

// ---- Large Strain Truss Elements ----
class T2D2NL : public LargeStrainTrussElement {};    // 2-node truss 2D - large strain
class T3D2NL : public LargeStrainTrussElement {};    // 2-node truss 3D - large strain

// ============================================================================
// ELEMENT FACTORY WITH STRAIN FORMULATION SELECTION
// ============================================================================

class ElementFactory {
public:
    enum class StrainFormulation {
        SMALL_STRAIN,           // Linear geometric analysis
        LARGE_STRAIN,           // Nonlinear geometric analysis (infinite strain capable)
        AUTO_DETECT            // Automatically switch based on deformation
    };
    
    enum class ElementType {
        // 2D Elements
        CPE3, CPE4, CPE4R, CPE4I, CPE6, CPE8, CPE8R,
        CPS3, CPS4, CPS4R, CPS4I, CPS6, CPS8, CPS8R,
        CAX3, CAX4, CAX4R, CAX4I, CAX6, CAX8, CAX8R,
        
        // 3D Elements
        C3D4, C3D6, C3D8, C3D8R, C3D8I, 
        C3D10, C3D15, C3D20, C3D20R, C3D27,
        
        // Enhanced Elements
        C3D8E, C3D20E, CPE4E, CPS4E,  // Enhanced strain
        C3D8F, C3D4F,                  // F-bar method
        
        // Structural Elements
        B21, B22, B31, B32,            // Beams
        S3, S4, S4R, S8R,              // Shells
        T2D2, T3D2                     // Trusses
    };
    
    static std::shared_ptr<Element> createElement(
        int id,
        ElementType element_type,
        const std::vector<std::shared_ptr<Node>>& nodes,
        std::shared_ptr<Material> material,
        StrainFormulation formulation = StrainFormulation::SMALL_STRAIN);
};

// ============================================================================
// SPECIAL PURPOSE ELEMENTS (Usually strain formulation independent)
// ============================================================================

class SpecialPurposeElement : public Element {};

// These elements typically don't distinguish between strain formulations
class RigidElement : public SpecialPurposeElement {};
class R2D2 : public RigidElement {};
class R3D3 : public RigidElement {};
class R3D4 : public RigidElement {};

class SpringElement : public SpecialPurposeElement {};
class SPRING1 : public SpringElement {};
class SPRING2 : public SpringElement {};

class MassElement : public SpecialPurposeElement {};
class MASS : public MassElement {};
class ROTARYI : public MassElement {};

class ConnectorElement : public SpecialPurposeElement {};
class CONN2D2 : public ConnectorElement {};
class CONN3D2 : public ConnectorElement {};

// Cohesive elements can benefit from large strain formulation
class CohesiveElement : public SpecialPurposeElement {};
class COH2D4 : public CohesiveElement {};      // Small strain cohesive
class COH2D4NL : public CohesiveElement {};    // Large strain cohesive
class COH3D6 : public CohesiveElement {};      // Small strain cohesive
class COH3D6NL : public CohesiveElement {};    // Large strain cohesive
class COH3D8 : public CohesiveElement {};      // Small strain cohesive
class COH3D8NL : public CohesiveElement {};    // Large strain cohesive