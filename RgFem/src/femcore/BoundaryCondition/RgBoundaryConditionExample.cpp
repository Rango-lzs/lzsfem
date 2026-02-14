/*********************************************************************
 * \file   RgBoundaryConditionExample.cpp
 * \brief  Examples of using Rg boundary condition system
 *
 * \author Leizs
 * \date   January 2025
 *********************************************************************/

#include "RgBoundaryCondition.h"
#include "RgLoadController.h"
#include "AbaqusBoundaryParser.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include "femcore/FENodeSet.h"

//-----------------------------------------------------------------------------
// Example 1: Create boundary conditions with load curves
//-----------------------------------------------------------------------------
void Example1_LoadCurveBC(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    
    // Create a piecewise linear load curve
    RgLoadCurve* lc = new RgLoadCurve();
    lc->AddPoint(0.0, 0.0);    // t=0: value=0
    lc->AddPoint(0.5, 1.0);    // t=0.5: value=1 (ramp up)
    lc->AddPoint(1.0, 1.0);    // t=1: value=1 (hold)
    lc->AddPoint(1.5, 0.0);    // t=1.5: value=0 (ramp down)
    lc->SetInterpolation(RgLoadController::INTERP_LINEAR);
    lc->SetExtendMode(RgLoadController::EXTEND_CONSTANT);
    lc->SetName("RampLoadCurve");
    fem->AddLoadController(lc);
    
    // Create prescribed displacement BC with load curve
    FENodeSet* nodeSet = mesh.FindNodeSet("LoadedFace");
    if (nodeSet)
    {
        RgPrescribedDisplacement* bc = new RgPrescribedDisplacement(fem);
        bc->SetNodeSet(nodeSet);
        bc->SetDOF(2);           // Z direction
        bc->SetScale(0.1);       // Maximum displacement = 0.1
        bc->SetLoadController(lc);
        bc->SetName("TimeVaryingLoad");
        fem->AddBoundaryCondition(bc);
        
        // At t=0.25: displacement = 0.1 * 0.5 = 0.05
        // At t=0.5:  displacement = 0.1 * 1.0 = 0.1
        // At t=1.0:  displacement = 0.1 * 1.0 = 0.1
        // At t=1.25: displacement = 0.1 * 0.5 = 0.05
    }
}

//-----------------------------------------------------------------------------
// Example 2: Different load controller types
//-----------------------------------------------------------------------------
void Example2_ControllerTypes(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    FENodeSet* nodeSet = mesh.FindNodeSet("LoadPoint");
    
    // 2a: Linear ramp
    {
        RgLinearController* lc = new RgLinearController();
        lc->SetParameters(0.0, 1.0, 0.0, 1.0);  // Ramp from 0 to 1 over 0-1 seconds
        lc->SetName("LinearRamp");
        fem->AddLoadController(lc);
        
        RgPrescribedDisplacement* bc = new RgPrescribedDisplacement(fem);
        bc->SetNodeSet(nodeSet);
        bc->SetDOF(0);
        bc->SetScale(0.05);
        bc->SetLoadController(lc);
        bc->SetName("LinearDisplacement");
        fem->AddBoundaryCondition(bc);
    }
    
    // 2b: Step function
    {
        RgStepController* lc = new RgStepController();
        lc->SetParameters(0.5, 0.0, 1.0);  // Step from 0 to 1 at t=0.5
        lc->SetName("StepLoad");
        fem->AddLoadController(lc);
    }
    
    // 2c: Smooth step (S-curve)
    {
        RgSmoothStepController* lc = new RgSmoothStepController();
        lc->SetParameters(0.0, 1.0, 0.0, 1.0);  // Smooth transition
        lc->SetName("SmoothStep");
        fem->AddLoadController(lc);
    }
    
    // 2d: Sinusoidal loading (dynamic)
    {
        RgSineController* lc = new RgSineController();
        lc->SetParameters(
            1.0,    // Amplitude
            2.0,    // Frequency (Hz)
            0.0,    // Phase
            0.0     // Offset
        );
        lc->SetName("CyclicLoad");
        fem->AddLoadController(lc);
    }
    
    // 2e: Constant value
    {
        RgConstantController* lc = new RgConstantController();
        lc->SetValue(1.0);
        lc->SetName("Constant");
        fem->AddLoadController(lc);
    }
}

//-----------------------------------------------------------------------------
// Example 3: Smooth interpolation for better convergence
//-----------------------------------------------------------------------------
void Example3_SmoothInterpolation(FEModel* fem)
{
    // Create smooth load curve for better convergence
    RgLoadCurve* lc = new RgLoadCurve();
    
    // Add control points
    lc->AddPoint(0.0, 0.0);
    lc->AddPoint(0.1, 0.1);
    lc->AddPoint(0.5, 0.8);
    lc->AddPoint(1.0, 1.0);
    
    // Use smooth (cubic) interpolation
    lc->SetInterpolation(RgLoadController::INTERP_SMOOTH);
    lc->SetExtendMode(RgLoadController::EXTEND_CONSTANT);
    
    fem->AddLoadController(lc);
}

//-----------------------------------------------------------------------------
// Example 4: Complete workflow with Abaqus parsing and load curves
//-----------------------------------------------------------------------------
void Example4_CompleteWorkflow(const char* inputFile)
{
    FEModel* fem = new FEModel();
    
    // 1. Parse Abaqus boundary conditions
    std::vector<std::string> lines;
    // ... read file into lines ...
    
    AbaqusBoundaryParser parser(fem);
    for (size_t i = 0; i < lines.size(); ++i)
    {
        if (lines[i].find("*Boundary") == 0)
        {
            int consumed = parser.ParseBoundary(lines, i);
            i += consumed - 1;
        }
    }
    parser.CreateBoundaryConditions();
    
    // 2. Create load curves for time-varying BCs
    RgLoadCurve* loadCurve = new RgLoadCurve();
    loadCurve->AddPoint(0.0, 0.0);
    loadCurve->AddPoint(1.0, 1.0);
    loadCurve->SetName("MainLoadCurve");
    fem->AddLoadController(loadCurve);
    
    // 3. Attach load curve to prescribed BCs
    int nbc = fem->BoundaryConditions();
    for (int i = 0; i < nbc; ++i)
    {
        RgBoundaryCondition* bc = dynamic_cast<RgBoundaryCondition*>(fem->BoundaryCondition(i));
        
        // Check if it's a prescribed displacement
        RgPrescribedDisplacement* pbc = dynamic_cast<RgPrescribedDisplacement*>(bc);
        if (pbc && pbc->GetScale() != 0.0)
        {
            pbc->SetLoadController(loadCurve);
        }
    }
    
    // 4. Initialize and solve
    fem->Init();
    
    // In each time step:
    double time = 0.0;
    double dt = 0.1;
    for (int step = 0; step < 10; ++step)
    {
        time += dt;
        
        // Evaluate all load controllers
        fem->EvaluateLoadControllers(time);
        
        // Update BCs
        fem->Update();
        
        // Solve...
    }
    
    delete fem;
}

//-----------------------------------------------------------------------------
// Example 5: Cyclic loading for fatigue analysis
//-----------------------------------------------------------------------------
void Example5_CyclicLoading(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    FENodeSet* nodeSet = mesh.FindNodeSet("CyclicLoadPoint");
    
    // Create sinusoidal load with mean stress
    RgSineController* lc = new RgSineController();
    lc->SetParameters(
        0.5,    // Amplitude = 50% of max
        1.0,    // Frequency = 1 Hz
        0.0,    // Phase
        0.5     // Offset = 50% mean
    );
    // Result: value oscillates between 0 and 1
    
    lc->SetName("FatigueLoad");
    fem->AddLoadController(lc);
    
    RgPrescribedDisplacement* bc = new RgPrescribedDisplacement(fem);
    bc->SetNodeSet(nodeSet);
    bc->SetDOF(1);
    bc->SetScale(0.01);  // 1cm max
    bc->SetLoadController(lc);
    bc->SetName("CyclicDisplacement");
    fem->AddBoundaryCondition(bc);
}

//-----------------------------------------------------------------------------
// Example 6: Multi-stage loading
//-----------------------------------------------------------------------------
void Example6_MultiStageLoading(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    FENodeSet* nodeSet = mesh.FindNodeSet("LoadPoint");
    
    // Stage 1: Ramp up (0-10s)
    // Stage 2: Hold (10-20s)
    // Stage 3: Ramp down (20-30s)
    // Stage 4: Reverse (30-40s)
    
    RgLoadCurve* lc = new RgLoadCurve();
    lc->AddPoint(0.0, 0.0);     // Start
    lc->AddPoint(10.0, 1.0);    // End of ramp up
    lc->AddPoint(20.0, 1.0);    // End of hold
    lc->AddPoint(30.0, 0.0);    // End of ramp down
    lc->AddPoint(40.0, -0.5);   // Reverse loading
    
    lc->SetInterpolation(RgLoadController::INTERP_LINEAR);
    lc->SetName("MultiStageLoad");
    fem->AddLoadController(lc);
    
    RgPrescribedDisplacement* bc = new RgPrescribedDisplacement(fem);
    bc->SetNodeSet(nodeSet);
    bc->SetDOF(2);
    bc->SetScale(0.1);
    bc->SetLoadController(lc);
    bc->SetName("MultiStageBC");
    fem->AddBoundaryCondition(bc);
}

//-----------------------------------------------------------------------------
// Example 7: Using extend modes
//-----------------------------------------------------------------------------
void Example7_ExtendModes(FEModel* fem)
{
    // Curve defined for t=[0, 1]
    RgLoadCurve* lc1 = new RgLoadCurve();
    lc1->AddPoint(0.0, 0.0);
    lc1->AddPoint(1.0, 1.0);
    lc1->SetExtendMode(RgLoadController::EXTEND_CONSTANT);
    // For t<0: value=0, for t>1: value=1
    
    RgLoadCurve* lc2 = new RgLoadCurve();
    lc2->AddPoint(0.0, 0.0);
    lc2->AddPoint(1.0, 1.0);
    lc2->SetExtendMode(RgLoadController::EXTEND_EXTRAPOLATE);
    // For t<0: value extrapolated (negative)
    // For t>1: value extrapolated (>1)
    
    RgLoadCurve* lc3 = new RgLoadCurve();
    lc3->AddPoint(0.0, 0.0);
    lc3->AddPoint(0.5, 1.0);
    lc3->AddPoint(1.0, 0.0);
    lc3->SetExtendMode(RgLoadController::EXTEND_REPEAT);
    // Repeats the curve: t=1.5 -> value at t=0.5 (=1.0)
    
    fem->AddLoadController(lc1);
    fem->AddLoadController(lc2);
    fem->AddLoadController(lc3);
}

//-----------------------------------------------------------------------------
// Example 8: Programmatic BC creation (no parsing)
//-----------------------------------------------------------------------------
void Example8_ProgrammaticCreation(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    
    // Fixed BCs
    {
        FENodeSet* fixedNodes = mesh.FindNodeSet("FixedEnd");
        RgFixedBC* bc = RgBCFactory::CreateEncastre(fem, fixedNodes);
        bc->SetName("FixedSupport");
        fem->AddBoundaryCondition(bc);
    }
    
    // Symmetry BC
    {
        FENodeSet* symNodes = mesh.FindNodeSet("SymmetryPlane");
        Vector3d normal(1, 0, 0);  // YZ plane
        RgSymmetryBC* bc = RgBCFactory::CreateSymmetry(fem, symNodes, normal);
        bc->SetName("Symmetry_X");
        fem->AddBoundaryCondition(bc);
    }
    
    // Prescribed displacement with smooth load curve
    {
        FENodeSet* loadNodes = mesh.FindNodeSet("LoadedFace");
        
        RgLoadCurve* lc = new RgLoadCurve();
        lc->AddPoint(0.0, 0.0);
        lc->AddPoint(1.0, 1.0);
        lc->SetInterpolation(RgLoadController::INTERP_SMOOTH);
        fem->AddLoadController(lc);
        
        RgPrescribedDisplacement* bc = RgBCFactory::CreatePrescribed(
            fem, loadNodes, 2, 0.05);
        bc->SetLoadController(lc);
        bc->SetName("LoadedDisplacement");
        fem->AddBoundaryCondition(bc);
    }
}

//-----------------------------------------------------------------------------
// Main example
//-----------------------------------------------------------------------------
int main()
{
    FEModel fem;
    
    // Run examples
    // Example1_LoadCurveBC(&fem);
    // Example2_ControllerTypes(&fem);
    // Example3_SmoothInterpolation(&fem);
    
    return 0;
}
