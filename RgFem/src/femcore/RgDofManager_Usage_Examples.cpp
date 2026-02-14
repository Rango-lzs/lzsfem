/*********************************************************************
 * \file   RgDofSchema_Usage_Examples.cpp
 * \brief  Usage examples for RgDofSchema system
 *
 * \date   February 2026
 *********************************************************************/

#include "RgDofSchema.h"
#include "femcore/FEModel.h"

 //=============================================================================
 // EXAMPLE 1: Basic 3D Solid Mechanics
 //=============================================================================

void Example1_Solid3D()
{
    // Create and configure DOF schema for 3D solid
    RgDofSchema schema;
    schema.SetSchema(RgDofSet::Solid3D());

    // Query schema
    int dofsPerNode = schema.GetDofsPerNode();  // Returns 3
    int uIndex = schema.GetDofIndex("u");       // Returns 0
    int vIndex = schema.GetDofIndex("v");       // Returns 1
    int wIndex = schema.GetDofIndex("w");       // Returns 2

    // Check properties
    bool is3D = schema.Is3D();                  // Returns true
    bool hasRotation = schema.HasRotationalDofs();  // Returns false

    // Print schema
    schema.PrintSchema();

    /* Output:
    =================================================
    DOF Schema Configuration
    =================================================
    DOFs per node: 3
    -------------------------------------------------
    Index  Name    Description
    -------------------------------------------------
       0   u       X displacement
       1   v       Y displacement
       2   w       Z displacement
    =================================================
    */
}

//=============================================================================
// EXAMPLE 2: 3D Beam/Shell Analysis
//=============================================================================

void Example2_BeamShell()
{
    // Create schema with rotational DOFs
    RgDofSchema schema;
    schema.SetSchema(RgDofSet::Beam3D());

    int dofsPerNode = schema.GetDofsPerNode();  // Returns 6

    // Access rotational DOFs
    int rxIndex = schema.GetDofIndex(DofType::ROTATION_X);  // Returns 3
    int ryIndex = schema.GetDofIndex(DofType::ROTATION_Y);  // Returns 4
    int rzIndex = schema.GetDofIndex(DofType::ROTATION_Z);  // Returns 5

    bool hasRotation = schema.HasRotationalDofs();  // Returns true
}

//=============================================================================
// EXAMPLE 3: Thermo-Mechanical Coupling
//=============================================================================

void Example3_ThermoMechanical()
{
    // Coupled thermo-mechanical analysis
    RgDofSchema schema;
    schema.SetSchema(RgDofSet::Coupled_ThermoMech());

    int dofsPerNode = schema.GetDofsPerNode();  // Returns 4

    // DOF indices
    int uIndex = schema.GetDofIndex("u");   // 0
    int vIndex = schema.GetDofIndex("v");   // 1
    int wIndex = schema.GetDofIndex("w");   // 2
    int tIndex = schema.GetDofIndex("T");   // 3

    bool hasStructural = schema.HasStructuralDofs();  // true
    bool hasThermal = schema.HasThermalDofs();        // true
}

//=============================================================================
// EXAMPLE 4: Using Presets
//=============================================================================

void Example4_Presets()
{
    // Quick setup using presets

    // 3D structural
    RgDofSchema* schema1 = RgDofSchemaPresets::CreateStructural3D();
    // DOFs: u, v, w

    // Beam/Shell
    RgDofSchema* schema2 = RgDofSchemaPresets::CreateBeamShell3D();
    // DOFs: u, v, w, Rx, Ry, Rz

    // Thermal
    RgDofSchema* schema3 = RgDofSchemaPresets::CreateThermal();
    // DOFs: T

    // Cleanup
    delete schema1;
    delete schema2;
    delete schema3;
}

//=============================================================================
// EXAMPLE 5: Integration with FEModel
//=============================================================================

void Example5_FEModelIntegration()
{
    FEModel model;

    // Get DOF schema from model
    RgDofSchema& schema = model.GetDofSchema();

    // Configure schema
    schema.SetSchema(RgDofSet::Solid3D());

    // Validate
    if (!schema.Validate())
    {
        RgLogError("Invalid DOF schema");
        return;
    }

    // Print configuration
    schema.PrintSchema();
}

//=============================================================================
// EXAMPLE 6: Element Compatibility Check
//=============================================================================

void Example6_ElementCompatibility()
{
    RgDofSchema schema;
    schema.SetSchema(RgDofSet::Shell());  // 6 DOFs

    // Check if solid element is compatible
    RgDofSet solidDofs = RgDofSet::Solid3D();
    bool compatible = schema.IsCompatible(solidDofs);  // true

    // Get indices for solid elements (will use subset)
    std::vector<int> indices = schema.GetIndices(solidDofs);
    // Returns: {0, 1, 2} - only displacement DOFs
}

//=============================================================================
// EXAMPLE 7: Mixed Element Mesh
//=============================================================================

void Example7_MixedElements()
{
    // Mesh with solid + shell elements
    RgDofSchema schema;
    schema.SetSchema(RgDofSet::Shell());  // Use maximum DOF set

    // Solid elements use only {u, v, w}
    RgDofSet solidDofs = RgDofSet::Solid3D();
    std::vector<int> solidIndices = schema.GetIndices(solidDofs);
    // Returns: {0, 1, 2}

    // Shell elements use all {u, v, w, Rx, Ry, Rz}
    RgDofSet shellDofs = RgDofSet::Shell();
    std::vector<int> shellIndices = schema.GetIndices(shellDofs);
    // Returns: {0, 1, 2, 3, 4, 5}
}

//=============================================================================
// EXAMPLE 8: Custom DOF Configuration
//=============================================================================

void Example8_CustomConfiguration()
{
    RgDofSchema schema;

    // Add DOFs manually for 2D + thermal
    schema.AddDof(DofType::DISPLACEMENT_X);
    schema.AddDof(DofType::DISPLACEMENT_Y);
    schema.AddDof(DofType::TEMPERATURE);

    int dofsPerNode = schema.GetDofsPerNode();  // 3
    bool is2D = schema.Is2D();                  // true
    bool hasThermal = schema.HasThermalDofs();  // true
}

//=============================================================================
// EXAMPLE 9: Boundary Condition Application
//=============================================================================

void Example9_BoundaryConditions(FEModel& model, int nodeId)
{
    RgDofSchema& schema = model.GetDofSchema();

    // Fix u = 0 at a node
    int uIndex = schema.GetDofIndex("u");
    if (uIndex >= 0)
    {
        int dofsPerNode = schema.GetDofsPerNode();
        int globalDof = nodeId * dofsPerNode + uIndex;

        // Apply BC: PrescribeDOF(globalDof, 0.0);
    }

    // Fix all DOFs (encastre)
    int dofsPerNode = schema.GetDofsPerNode();
    for (int d = 0; d < dofsPerNode; ++d)
    {
        int globalDof = nodeId * dofsPerNode + d;
        // Apply BC: PrescribeDOF(globalDof, 0.0);
    }
}

//=============================================================================
// EXAMPLE 10: Element Stiffness Assembly
//=============================================================================

void Example10_ElementAssembly(FEModel& model)
{
    RgDofSchema& schema = model.GetDofSchema();

    // For a hex8 element (8 nodes)
    int nNodes = 8;
    int dofsPerNode = schema.GetDofsPerNode();  // 3 for solid
    int totalElemDofs = nNodes * dofsPerNode;   // 24

    // Element stiffness matrix
    Matrix Ke(totalElemDofs, totalElemDofs);

    // ... compute Ke ...

    // Assemble into global system
    std::vector<int> nodeIds = { 1, 2, 3, 4, 5, 6, 7, 8 };  // Element nodes

    for (int i = 0; i < nNodes; ++i)
    {
        int nodeId = nodeIds[i];

        for (int d = 0; d < dofsPerNode; ++d)
        {
            int globalDof = nodeId * dofsPerNode + d;
            int localDof = i * dofsPerNode + d;

            // Assemble: K[globalDof][globalDof] += Ke[localDof][localDof]
        }
    }
}

//=============================================================================
// PRACTICAL: Complete Model Setup
//=============================================================================

void PracticalExample_ModelSetup()
{
    // 1. Create model
    FEModel model;

    // 2. Configure DOF schema
    RgDofSchema& schema = model.GetDofSchema();
    schema.SetSchema(RgDofSet::Solid3D());

    // 3. Validate
    if (!schema.Validate())
    {
        RgLogError("Invalid DOF configuration");
        return;
    }

    // 4. Print configuration
    schema.PrintSchema();

    // 5. Get DOF information for use throughout the model
    int dofsPerNode = schema.GetDofsPerNode();

    // 6. Now ready to:
    //    - Create elements
    //    - Apply boundary conditions
    //    - Assemble system matrices
    //    - Solve
}