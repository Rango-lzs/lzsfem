/**
 * Element Class Hierarchy Based on ABAQUS Element Library
 * Reference: ABAQUS Element Library Documentation
 */

// ============================================================================
// ABSTRACT BASE CLASSES
// ============================================================================

class Element {
    // Root base class for all elements
};

class ContinuumElement : public Element {
    // Base for solid/continuum elements
};

class StructuralElement : public Element {
    // Base for beam, shell, truss elements
};

class SpecialPurposeElement : public Element {
    // Base for specialized elements
};

// ============================================================================
// CONTINUUM ELEMENTS (Stress/Displacement)
// ============================================================================

// ---- 2D Continuum Elements ----
class CPE_Element : public ContinuumElement {
    // 2D Plane Strain elements
};

class CPS_Element : public ContinuumElement {
    // 2D Plane Stress elements  
};

class CAX_Element : public ContinuumElement {
    // 2D Axisymmetric elements
};

// Linear 2D Elements
class CPE3 : public CPE_Element {};      // 3-node triangular plane strain
class CPE4 : public CPE_Element {};      // 4-node quadrilateral plane strain
class CPE4R : public CPE_Element {};     // 4-node reduced integration plane strain
class CPE4I : public CPE_Element {};     // 4-node incompatible mode plane strain
class CPE6 : public CPE_Element {};      // 6-node triangular plane strain
class CPE8 : public CPE_Element {};      // 8-node quadrilateral plane strain
class CPE8R : public CPE_Element {};     // 8-node reduced integration plane strain

class CPS3 : public CPS_Element {};      // 3-node triangular plane stress
class CPS4 : public CPS_Element {};      // 4-node quadrilateral plane stress
class CPS4R : public CPS_Element {};     // 4-node reduced integration plane stress
class CPS4I : public CPS_Element {};     // 4-node incompatible mode plane stress
class CPS6 : public CPS_Element {};      // 6-node triangular plane stress
class CPS8 : public CPS_Element {};      // 8-node quadrilateral plane stress
class CPS8R : public CPS_Element {};     // 8-node reduced integration plane stress

class CAX3 : public CAX_Element {};      // 3-node triangular axisymmetric
class CAX4 : public CAX_Element {};      // 4-node quadrilateral axisymmetric
class CAX4R : public CAX_Element {};     // 4-node reduced integration axisymmetric
class CAX4I : public CAX_Element {};     // 4-node incompatible mode axisymmetric
class CAX6 : public CAX_Element {};      // 6-node triangular axisymmetric
class CAX8 : public CAX_Element {};      // 8-node quadrilateral axisymmetric
class CAX8R : public CAX_Element {};     // 8-node reduced integration axisymmetric

// ---- 3D Continuum Elements ----
class C3D_Element : public ContinuumElement {
    // 3D Stress elements
};

// Linear 3D Elements
class C3D4 : public C3D_Element {};      // 4-node tetrahedral
class C3D6 : public C3D_Element {};      // 6-node triangular prism
class C3D8 : public C3D_Element {};      // 8-node hexahedral
class C3D8R : public C3D_Element {};     // 8-node reduced integration hexahedral
class C3D8I : public C3D_Element {};     // 8-node incompatible mode hexahedral
class C3D8RH : public C3D_Element {};    // 8-node hybrid hexahedral
class C3D8H : public C3D_Element {};     // 8-node hybrid hexahedral

// Quadratic 3D Elements
class C3D10 : public C3D_Element {};     // 10-node tetrahedral
class C3D10M : public C3D_Element {};    // 10-node modified tetrahedral
class C3D15 : public C3D_Element {};     // 15-node triangular prism
class C3D20 : public C3D_Element {};     // 20-node hexahedral
class C3D20R : public C3D_Element {};    // 20-node reduced integration hexahedral
class C3D20RH : public C3D_Element {};   // 20-node hybrid hexahedral
class C3D27 : public C3D_Element {};     // 27-node hexahedral

// Enhanced 3D Elements
class C3D4H : public C3D_Element {};     // 4-node hybrid tetrahedral
class C3D6H : public C3D_Element {};     // 6-node hybrid prism
class C3D10H : public C3D_Element {};    // 10-node hybrid tetrahedral

// ============================================================================
// STRUCTURAL ELEMENTS
// ============================================================================

// ---- Beam Elements ----
class BeamElement : public StructuralElement {
    // Base for beam elements
};

class B21 : public BeamElement {};       // 2-node linear beam in 2D
class B22 : public BeamElement {};       // 3-node quadratic beam in 2D
class B23 : public BeamElement {};       // 4-node cubic beam in 2D
class B31 : public BeamElement {};       // 2-node linear beam in 3D
class B32 : public BeamElement {};       // 3-node quadratic beam in 3D
class B33 : public BeamElement {};       // 4-node cubic beam in 3D

class PIPE21 : public BeamElement {};    // 2-node linear pipe
class PIPE22 : public BeamElement {};    // 3-node quadratic pipe
class PIPE31 : public BeamElement {};    // 2-node linear pipe in 3D
class PIPE32 : public BeamElement {};    // 3-node quadratic pipe in 3D

// ---- Shell Elements ----
class ShellElement : public StructuralElement {
    // Base for shell elements
};

// General Purpose Shells
class S3 : public ShellElement {};       // 3-node triangular general purpose shell
class S4 : public ShellElement {};       // 4-node general purpose shell
class S4R : public ShellElement {};      // 4-node reduced integration shell
class S4RS : public ShellElement {};     // 4-node reduced integration with drilling DOF
class S4RSW : public ShellElement {};    // 4-node reduced integration with warping
class S8R : public ShellElement {};      // 8-node reduced integration shell
class S8R5 : public ShellElement {};     // 8-node reduced integration (5 DOF)
class S9R5 : public ShellElement {};     // 9-node reduced integration (5 DOF)

// Thick Shell Elements
class SC6R : public ShellElement {};     // 6-node thick shell
class SC8R : public ShellElement {};     // 8-node thick shell

// Conventional Shell Elements
class STRI3 : public ShellElement {};    // 3-node triangular facet shell
class STRI65 : public ShellElement {};   // 6-node triangular thin shell
class S3R : public ShellElement {};      // 3-node triangular reduced integration
class S6 : public ShellElement {};       // 6-node triangular shell

// ---- Membrane Elements ----
class MembraneElement : public StructuralElement {
    // Base for membrane elements
};

class M3D3 : public MembraneElement {};  // 3-node triangular membrane
class M3D4 : public MembraneElement {};  // 4-node quadrilateral membrane
class M3D4R : public MembraneElement {}; // 4-node reduced integration membrane
class M3D6 : public MembraneElement {};  // 6-node triangular membrane
class M3D8 : public MembraneElement {};  // 8-node quadrilateral membrane
class M3D8R : public MembraneElement {}; // 8-node reduced integration membrane
class M3D9 : public MembraneElement {};  // 9-node quadrilateral membrane
class M3D9R : public MembraneElement {}; // 9-node reduced integration membrane

// ---- Truss/Tie Elements ----
class TrussElement : public StructuralElement {
    // Base for truss elements
};

class T2D2 : public TrussElement {};     // 2-node linear truss (2D)
class T2D3 : public TrussElement {};     // 3-node quadratic truss (2D)
class T3D2 : public TrussElement {};     // 2-node linear truss (3D)
class T3D3 : public TrussElement {};     // 3-node quadratic truss (3D)

// ============================================================================
// SPECIAL PURPOSE ELEMENTS
// ============================================================================

// ---- Rigid Elements ----
class RigidElement : public SpecialPurposeElement {
    // Base for rigid body elements
};

class R2D2 : public RigidElement {};     // 2-node rigid link (2D)
class R3D3 : public RigidElement {};     // 3-node rigid triangle
class R3D4 : public RigidElement {};     // 4-node rigid quadrilateral

// ---- Spring/Dashpot Elements ----
class SpringElement : public SpecialPurposeElement {
    // Base for spring elements
};

class SPRING1 : public SpringElement {}; // Grounded spring/dashpot
class SPRING2 : public SpringElement {}; // Spring/dashpot between nodes
class SPRINGA : public SpringElement {}; // Assembled spring
class DASHPOT1 : public SpringElement {};// Grounded dashpot
class DASHPOT2 : public SpringElement {};// Dashpot between nodes
class DASHPOTA : public SpringElement {};// Assembled dashpot

// ---- Mass Elements ----
class MassElement : public SpecialPurposeElement {
    // Base for mass elements
};

class MASS : public MassElement {};      // Point mass
class ROTARYI : public MassElement {};   // Rotary inertia

// ---- Connector Elements ----
class ConnectorElement : public SpecialPurposeElement {
    // Base for connector elements
};

class CONN2D2 : public ConnectorElement {};  // 2D connector
class CONN3D2 : public ConnectorElement {};  // 3D connector

// ---- Gasket Elements ----
class GasketElement : public SpecialPurposeElement {
    // Base for gasket elements
};

class GKPS4 : public GasketElement {};   // 4-node plane stress gasket
class GKPE4 : public GasketElement {};   // 4-node plane strain gasket
class GKAX4 : public GasketElement {};   // 4-node axisymmetric gasket
class GK3D8 : public GasketElement {};   // 8-node 3D gasket

// ---- Cohesive Elements ----
class CohesiveElement : public SpecialPurposeElement {
    // Base for cohesive zone elements
};

class COH2D4 : public CohesiveElement {};    // 4-node 2D cohesive
class COH3D6 : public CohesiveElement {};    // 6-node 3D cohesive (prism)
class COH3D8 : public CohesiveElement {};    // 8-node 3D cohesive (hexahedral)

// ============================================================================
// INFINITE ELEMENTS
// ============================================================================

class InfiniteElement : public SpecialPurposeElement {
    // Base for infinite elements
};

class CIN2D4 : public InfiniteElement {};    // 4-node 2D infinite
class CIN3D8 : public InfiniteElement {};    // 8-node 3D infinite

// ============================================================================
// ACOUSTIC ELEMENTS
// ============================================================================

class AcousticElement : public SpecialPurposeElement {
    // Base for acoustic elements
};

class AC2D3 : public AcousticElement {};     // 3-node 2D acoustic
class AC2D4 : public AcousticElement {};     // 4-node 2D acoustic
class AC3D4 : public AcousticElement {};     // 4-node 3D acoustic
class AC3D6 : public AcousticElement {};     // 6-node 3D acoustic
class AC3D8 : public AcousticElement {};     // 8-node 3D acoustic
class AC3D10 : public AcousticElement {};    // 10-node 3D acoustic
class AC3D15 : public AcousticElement {};    // 15-node 3D acoustic
class AC3D20 : public AcousticElement {};    // 20-node 3D acoustic

// ============================================================================
// HEAT TRANSFER ELEMENTS
// ============================================================================

class HeatTransferElement : public SpecialPurposeElement {
    // Base for heat transfer elements
};

class DC2D3 : public HeatTransferElement {}; // 3-node 2D heat transfer
class DC2D4 : public HeatTransferElement {}; // 4-node 2D heat transfer  
class DC2D6 : public HeatTransferElement {}; // 6-node 2D heat transfer
class DC2D8 : public HeatTransferElement {}; // 8-node 2D heat transfer
class DC3D4 : public HeatTransferElement {}; // 4-node 3D heat transfer
class DC3D6 : public HeatTransferElement {}; // 6-node 3D heat transfer
class DC3D8 : public HeatTransferElement {}; // 8-node 3D heat transfer
class DC3D10 : public HeatTransferElement {};// 10-node 3D heat transfer
class DC3D15 : public HeatTransferElement {};// 15-node 3D heat transfer  
class DC3D20 : public HeatTransferElement {};// 20-node 3D heat transfer

class DCAX3 : public HeatTransferElement {}; // 3-node axisymmetric heat transfer
class DCAX4 : public HeatTransferElement {}; // 4-node axisymmetric heat transfer
class DCAX6 : public HeatTransferElement {}; // 6-node axisymmetric heat transfer
class DCAX8 : public HeatTransferElement {}; // 8-node axisymmetric heat transfer

// ============================================================================
// POROUS MEDIA ELEMENTS  
// ============================================================================

class PorousElement : public SpecialPurposeElement {
    // Base for porous media elements
};

class CPE4P : public PorousElement {};   // 4-node pore pressure plane strain
class CPS4P : public PorousElement {};   // 4-node pore pressure plane stress
class CAX4P : public PorousElement {};   // 4-node pore pressure axisymmetric
class C3D8P : public PorousElement {};   // 8-node pore pressure 3D

// ============================================================================
// USER-DEFINED ELEMENTS
// ============================================================================

class UserElement : public Element {
    // Base for user-defined elements
};

class U1 : public UserElement {};        // User element with 1 node
class U2 : public UserElement {};        // User element with 2 nodes
class U3 : public UserElement {};        // User element with 3 nodes
class U4 : public UserElement {};        // User element with 4 nodes
// ... up to U20+ for different node counts