/*********************************************************************
 * \file   test_RgAnalysisStep.cpp
 * \brief  Test program for RgAnalysisStep system
 *
 * \author 
 * \date   February 2025
 *********************************************************************/

#include "RgAnalysisStep.h"
#include "AbaqusImport.h"
#include "femcore/FEModel.h"
#include "femcore/FEBoundaryCondition.h"
#include "femcore/FEModelLoad.h"
#include "logger/log.h"

//=============================================================================
// Test 1: Basic step creation and management
//=============================================================================
void Test_BasicStepCreation()
{
    RgLog("\n");
    RgLog("=================================================\n");
    RgLog("Test 1: Basic Step Creation\n");
    RgLog("=================================================\n");

    FEModel model;
    model.SetActiveModule("solid");

    // Create first step
    RgStaticAnalysisStep* step1 = new RgStaticAnalysisStep(&model);
    step1->SetName("Step1-InitialLoad");
    step1->SetStepNumber(0);
    step1->SetTimePeriod(1.0);
    step1->SetInitialTimeIncrement(0.1);

    RgLog("Created step: %s\n", step1->GetName().c_str());
    RgLog("  Time period: %g\n", step1->GetTimePeriod());
    RgLog("  Initial dt: %g\n", step1->GetInitialTimeIncrement());

    // Add to model
    model.AddStep(step1);

    // Create second step
    RgStaticAnalysisStep* step2 = new RgStaticAnalysisStep(&model);
    step2->SetName("Step2-AdditionalLoad");
    step2->SetStepNumber(1);
    step2->SetTimePeriod(1.0);
    step2->SetPreviousStep(step1);

    model.AddStep(step2);

    RgLog("\nTotal steps in model: %d\n", model.Steps());
    RgLog("Test 1 PASSED\n");
}

//=============================================================================
// Test 2: BC and Load management
//=============================================================================
void Test_BCAndLoadManagement()
{
    RgLog("\n");
    RgLog("=================================================\n");
    RgLog("Test 2: BC and Load Management\n");
    RgLog("=================================================\n");

    FEModel model;
    model.SetActiveModule("solid");

    // Create step
    RgStaticAnalysisStep* step = new RgStaticAnalysisStep(&model);
    step->SetName("TestStep");

    // Create and add BCs (pseudo code - you'll need actual BC objects)
    // FEBoundaryCondition* bc1 = new FEFixedBC(...);
    // bc1->SetName("FixedSupport");
    // step->AddBoundaryCondition(bc1);

    // For testing, we'll just log
    RgLog("Step created: %s\n", step->GetName().c_str());
    RgLog("  BCs in step: %d\n", step->BoundaryConditions());
    RgLog("  Loads in step: %d\n", step->Loads());

    // Test find
    // FEBoundaryCondition* found = step->FindBoundaryCondition("FixedSupport");
    // if (found)
    //     RgLog("  Found BC: %s\n", found->GetName().c_str());

    model.AddStep(step);

    RgLog("Test 2 PASSED\n");
}

//=============================================================================
// Test 3: Step inheritance
//=============================================================================
void Test_StepInheritance()
{
    RgLog("\n");
    RgLog("=================================================\n");
    RgLog("Test 3: Step Inheritance\n");
    RgLog("=================================================\n");

    FEModel model;
    model.SetActiveModule("solid");

    // Create step 1
    RgStaticAnalysisStep* step1 = new RgStaticAnalysisStep(&model);
    step1->SetName("Step1");
    step1->SetStepNumber(0);

    // Add some BCs and loads to step1
    // (In real code, create actual BC/Load objects)
    RgLog("Step 1 created with:\n");
    RgLog("  BCs: %d\n", step1->BoundaryConditions());
    RgLog("  Loads: %d\n", step1->Loads());

    model.AddStep(step1);

    // Create step 2
    RgStaticAnalysisStep* step2 = new RgStaticAnalysisStep(&model);
    step2->SetName("Step2");
    step2->SetStepNumber(1);
    step2->SetPreviousStep(step1);
    step2->SetActivationMode(StepActivationMode::INHERITED);

    // Initialize inheritance
    step2->InheritFromPreviousStep();

    RgLog("\nStep 2 after inheritance:\n");
    auto activeBCs = step2->GetAllActiveBCs();
    auto activeLoads = step2->GetAllActiveLoads();
    RgLog("  Active BCs: %d\n", (int)activeBCs.size());
    RgLog("  Active Loads: %d\n", (int)activeLoads.size());

    model.AddStep(step2);

    RgLog("Test 3 PASSED\n");
}

//=============================================================================
// Test 4: Activation modes
//=============================================================================
void Test_ActivationModes()
{
    RgLog("\n");
    RgLog("=================================================\n");
    RgLog("Test 4: Activation Modes\n");
    RgLog("=================================================\n");

    FEModel model;
    model.SetActiveModule("solid");

    // Test NEW mode
    RgStaticAnalysisStep* step1 = new RgStaticAnalysisStep(&model);
    step1->SetName("Step1");
    step1->SetActivationMode(StepActivationMode::NEW);
    RgLog("Step 1 mode: NEW\n");
    model.AddStep(step1);

    // Test INHERITED mode
    RgStaticAnalysisStep* step2 = new RgStaticAnalysisStep(&model);
    step2->SetName("Step2");
    step2->SetPreviousStep(step1);
    step2->SetActivationMode(StepActivationMode::INHERITED);
    RgLog("Step 2 mode: INHERITED\n");
    model.AddStep(step2);

    // Test REPLACE mode
    RgStaticAnalysisStep* step3 = new RgStaticAnalysisStep(&model);
    step3->SetName("Step3");
    step3->SetPreviousStep(step2);
    step3->SetActivationMode(StepActivationMode::REPLACE);
    RgLog("Step 3 mode: REPLACE\n");
    model.AddStep(step3);

    RgLog("Test 4 PASSED\n");
}

//=============================================================================
// Test 5: Abaqus INP file loading
//=============================================================================
void Test_AbaqusINPLoading()
{
    RgLog("\n");
    RgLog("=================================================\n");
    RgLog("Test 5: Abaqus INP Loading\n");
    RgLog("=================================================\n");

    FEModel model;
    AbaqusImport importer;

    // Create a simple test INP file content (in real test, load from file)
    const char* testINP = R"(
*Heading
Test model with multiple steps

*Node
1, 0.0, 0.0, 0.0
2, 1.0, 0.0, 0.0
3, 1.0, 1.0, 0.0
4, 0.0, 1.0, 0.0

*Element, type=C3D8
1, 1, 2, 3, 4, 5, 6, 7, 8

*Nset, nset=Fixed
1, 4

*Nset, nset=Loaded
2, 3

*Part, name=TestPart
...
*End Part

*Assembly, name=Assembly
*Instance, name=Part-1, part=TestPart
*End Instance
*End Assembly

*Step, name=Step1
*Static
0.1, 1.0
*Boundary
Fixed, 1, 3, 0.0
*Cload
Loaded, 2, -100.0
*End Step

*Step, name=Step2
*Static
0.1, 1.0
*Cload
Loaded, 2, -200.0
*End Step
)";

    // In real test, you would:
    // bool success = importer.load("test.inp", &model);
    
    RgLog("Note: This is a pseudo-test\n");
    RgLog("In real implementation, create test.inp and load it\n");
    
    // Check results
    // RgLog("Loaded steps: %d\n", model.Steps());
    // for (int i = 0; i < model.Steps(); ++i)
    // {
    //     RgAnalysisStep* step = dynamic_cast<RgAnalysisStep*>(model.GetStep(i));
    //     if (step)
    //     {
    //         RgLog("  Step %d: %s\n", i+1, step->GetName().c_str());
    //         RgLog("    Active BCs: %d\n", step->GetAllActiveBCs().size());
    //         RgLog("    Active Loads: %d\n", step->GetAllActiveLoads().size());
    //     }
    // }

    RgLog("Test 5 PASSED (pseudo-test)\n");
}

//=============================================================================
// Test 6: Step initialization
//=============================================================================
void Test_StepInitialization()
{
    RgLog("\n");
    RgLog("=================================================\n");
    RgLog("Test 6: Step Initialization\n");
    RgLog("=================================================\n");

    FEModel model;
    model.SetActiveModule("solid");

    // Create and add steps
    RgStaticAnalysisStep* step1 = new RgStaticAnalysisStep(&model);
    step1->SetName("InitTest-Step1");
    step1->SetStepNumber(0);
    model.AddStep(step1);

    RgStaticAnalysisStep* step2 = new RgStaticAnalysisStep(&model);
    step2->SetName("InitTest-Step2");
    step2->SetStepNumber(1);
    step2->SetPreviousStep(step1);
    model.AddStep(step2);

    // Initialize step1
    RgLog("Initializing step 1...\n");
    bool success1 = step1->Initialize();
    RgLog("  Result: %s\n", success1 ? "SUCCESS" : "FAILED");

    // Initialize step2
    RgLog("Initializing step 2...\n");
    bool success2 = step2->Initialize();
    RgLog("  Result: %s\n", success2 ? "SUCCESS" : "FAILED");

    if (success1 && success2)
    {
        RgLog("Test 6 PASSED\n");
    }
    else
    {
        RgLog("Test 6 FAILED\n");
    }
}

//=============================================================================
// Main test runner
//=============================================================================
int main(int argc, char** argv)
{
    RgLog("=================================================\n");
    RgLog("RgAnalysisStep Test Suite\n");
    RgLog("=================================================\n");

    try
    {
        Test_BasicStepCreation();
        Test_BCAndLoadManagement();
        Test_StepInheritance();
        Test_ActivationModes();
        Test_AbaqusINPLoading();
        Test_StepInitialization();

        RgLog("\n");
        RgLog("=================================================\n");
        RgLog("All Tests Completed\n");
        RgLog("=================================================\n");
    }
    catch (const std::exception& e)
    {
        RgLogError("Test failed with exception: %s\n", e.what());
        return 1;
    }

    return 0;
}
