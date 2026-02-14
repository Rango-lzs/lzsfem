/*********************************************************************
 * \file   RgLoadExample.cpp
 * \brief  Examples of using Rg load system
 *
 * \author Leizs
 * \date   January 2025
 *********************************************************************/

#include "RgLoad.h"
#include "RgLoadController.h"
#include "AbaqusLoadParser.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include "femcore/FENodeSet.h"
#include "femcore/FEFacetSet.h"

//-----------------------------------------------------------------------------
// Example 1: Nodal force loads
//-----------------------------------------------------------------------------
void Example1_NodalForces(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    
    // 1a: Single DOF force
    {
        FENodeSet* nodeSet = mesh.FindNodeSet("LoadPoint");
        RgNodalLoad* load = RgLoadFactory::CreateNodalLoad(
            fem, nodeSet, 0, 100.0);  // 100N in X direction
        load->SetName("PointForceX");
        fem->AddModelLoad(load);
    }
    
    // 1b: Vector force
    {
        FENodeSet* nodeSet = mesh.FindNodeSet("LoadPoints");
        Vector3d force(100.0, -50.0, 0.0);  // 100N in X, -50N in Y
        RgNodalLoad* load = RgLoadFactory::CreateNodalLoad(fem, nodeSet, force);
        load->SetName("MultiPointForce");
        fem->AddModelLoad(load);
    }
    
    // 1c: Time-varying force with load curve
    {
        FENodeSet* nodeSet = mesh.FindNodeSet("DynamicLoad");
        
        // Create ramp load curve
        RgLoadCurve* lc = new RgLoadCurve();
        lc->AddPoint(0.0, 0.0);
        lc->AddPoint(1.0, 1.0);
        lc->AddPoint(2.0, 0.5);
        lc->SetInterpolation(RgLoadController::INTERP_SMOOTH);
        fem->AddLoadController(lc);
        
        RgNodalLoad* load = RgLoadFactory::CreateNodalLoad(
            fem, nodeSet, 2, 1000.0);  // 1000N max in Z
        load->SetLoadController(lc);
        load->SetName("RampedForce");
        fem->AddModelLoad(load);
        
        // At t=0.5: force = 1000 * 0.5 = 500N
    }
}

//-----------------------------------------------------------------------------
// Example 2: Surface loads (pressure and traction)
//-----------------------------------------------------------------------------
void Example2_SurfaceLoads(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    
    // 2a: Pressure load
    {
        FEFacetSet* facetSet = mesh.FindFacetSet("TopSurface");
        RgSurfaceLoad* load = RgLoadFactory::CreatePressure(
            fem, facetSet, 1.0e6);  // 1 MPa pressure
        load->SetName("Pressure");
        fem->AddModelLoad(load);
    }
    
    // 2b: Traction vector
    {
        FEFacetSet* facetSet = mesh.FindFacetSet("SideSurface");
        FESurface* surface = mesh.CreateSurface(*facetSet);
        
        Vector3d traction(100.0, 0.0, 0.0);  // 100 Pa in X
        RgSurfaceLoad* load = RgLoadFactory::CreateTraction(fem, surface, traction);
        load->SetName("ShearLoad");
        fem->AddModelLoad(load);
    }
    
    // 2c: Follower pressure (follows deformed surface)
    {
        FEFacetSet* facetSet = mesh.FindFacetSet("Balloon");
        RgSurfaceLoad* load = RgLoadFactory::CreatePressure(fem, facetSet, 0.1e6);
        load->SetFollower(true);  // Pressure normal always perpendicular to surface
        load->SetName("InflationPressure");
        fem->AddModelLoad(load);
    }
}

//-----------------------------------------------------------------------------
// Example 3: Body loads
//-----------------------------------------------------------------------------
void Example3_BodyLoads(FEModel* fem)
{
    // 3a: Gravity
    {
        Vector3d g(0.0, 0.0, -9.81);  // g in -Z direction
        RgBodyLoad* load = RgLoadFactory::CreateGravity(fem, g);
        load->SetName("Gravity");
        fem->AddModelLoad(load);
    }
    
    // 3b: Centrifugal load
    {
        Vector3d axis(0, 0, 1);        // Rotation about Z axis
        Vector3d origin(0, 0, 0);       // At origin
        double omega = 100.0;           // 100 rad/s
        
        RgBodyLoad* load = RgLoadFactory::CreateCentrifugal(
            fem, axis, origin, omega);
        load->SetName("Centrifugal");
        fem->AddModelLoad(load);
    }
    
    // 3c: Time-varying gravity (earthquake simulation)
    {
        Vector3d g(1.0, 0.0, 0.0);  // Horizontal acceleration
        RgBodyLoad* load = RgLoadFactory::CreateGravity(fem, g);
        
        // Seismic load curve
        RgSineController* lc = new RgSineController();
        lc->SetParameters(
            9.81,   // Amplitude: 1g
            2.0,    // Frequency: 2 Hz
            0.0,    // Phase
            0.0     // No offset
        );
        fem->AddLoadController(lc);
        
        load->SetLoadController(lc);
        load->SetName("Seismic");
        fem->AddModelLoad(load);
    }
}

//-----------------------------------------------------------------------------
// Example 4: Moment loads
//-----------------------------------------------------------------------------
void Example4_MomentLoads(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    
    // Applied moment about Z axis
    FENodeSet* nodeSet = mesh.FindNodeSet("MomentPoint");
    Vector3d moment(0.0, 0.0, 100.0);  // 100 N⋅m about Z
    
    RgMomentLoad* load = RgLoadFactory::CreateMoment(fem, nodeSet, moment);
    load->SetName("AppliedMoment");
    fem->AddModelLoad(load);
}

//-----------------------------------------------------------------------------
// Example 5: Parse Abaqus loads
//-----------------------------------------------------------------------------
void Example5_ParseAbaqusLoads(FEModel* fem)
{
    std::vector<std::string> inputLines = {
        "** ----------------------------------------------------------------",
        "** LOADS",
        "** ----------------------------------------------------------------",
        "** Name: Load-1 Type: Concentrated force",
        "*Cload",
        "LoadPoint, 1, 100.0",
        "LoadPoint, 2, -50.0",
        "** Name: Load-2 Type: Pressure",
        "*Dload",
        "TopFace, P, 1.0e6",
        "** Name: Load-3 Type: Gravity",
        "*Dload",
        "AllElements, GRAV, 9.81, 0.0, 0.0, -1.0",
        "** End of loads"
    };
    
    // Create parser
    AbaqusLoadParser parser(fem);
    
    // Parse all load keywords
    for (size_t i = 0; i < inputLines.size(); ++i)
    {
        std::string line = inputLines[i];
        
        if (line.find("*Cload") == 0 || line.find("*CLOAD") == 0)
        {
            int consumed = parser.ParseCload(inputLines, i);
            i += consumed - 1;
        }
        else if (line.find("*Dload") == 0 || line.find("*DLOAD") == 0)
        {
            int consumed = parser.ParseDload(inputLines, i);
            i += consumed - 1;
        }
    }
    
    // Create loads
    if (!parser.CreateLoads())
    {
        // Error handling
        return;
    }
    
    // Display parsed loads
    const auto& loadData = parser.GetLoadData();
    std::cout << "Parsed " << loadData.size() << " loads:\n";
    for (size_t i = 0; i < loadData.size(); ++i)
    {
        const auto& data = loadData[i];
        if (data.isConcentrated)
        {
            std::cout << "  Cload: " << data.targetName 
                      << ", DOF " << data.dof 
                      << ", Magnitude " << data.magnitude << "\n";
        }
        else
        {
            std::cout << "  Dload: " << data.targetName
                      << ", Type " << data.loadType
                      << ", Magnitude " << data.magnitude << "\n";
        }
    }
}

//-----------------------------------------------------------------------------
// Example 6: Multi-step loading
//-----------------------------------------------------------------------------
void Example6_MultiStepLoading(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    FENodeSet* nodeSet = mesh.FindNodeSet("LoadPoint");
    
    // Stage 1: Linear ramp (0-10s)
    // Stage 2: Hold (10-20s)
    // Stage 3: Unload (20-30s)
    
    RgLoadCurve* lc = new RgLoadCurve();
    lc->AddPoint(0.0, 0.0);
    lc->AddPoint(10.0, 1.0);
    lc->AddPoint(20.0, 1.0);
    lc->AddPoint(30.0, 0.0);
    lc->SetName("MultiStageLoadCurve");
    fem->AddLoadController(lc);
    
    RgNodalLoad* load = RgLoadFactory::CreateNodalLoad(
        fem, nodeSet, 2, 1000.0);  // 1000N max
    load->SetLoadController(lc);
    load->SetName("MultiStageLoad");
    fem->AddModelLoad(load);
}

//-----------------------------------------------------------------------------
// Example 7: Combined loads
//-----------------------------------------------------------------------------
void Example7_CombinedLoads(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    
    // Gravity (always present)
    {
        Vector3d g(0, 0, -9.81);
        RgBodyLoad* load = RgLoadFactory::CreateGravity(fem, g);
        load->SetName("Gravity");
        fem->AddModelLoad(load);
    }
    
    // Point load (time-varying)
    {
        FENodeSet* nodeSet = mesh.FindNodeSet("LoadPoint");
        
        RgLoadCurve* lc = new RgLoadCurve();
        lc->AddPoint(0.0, 0.0);
        lc->AddPoint(1.0, 1.0);
        fem->AddLoadController(lc);
        
        RgNodalLoad* load = RgLoadFactory::CreateNodalLoad(
            fem, nodeSet, 2, 500.0);
        load->SetLoadController(lc);
        load->SetName("AppliedForce");
        fem->AddModelLoad(load);
    }
    
    // Pressure (constant)
    {
        FEFacetSet* facetSet = mesh.FindFacetSet("TopFace");
        RgSurfaceLoad* load = RgLoadFactory::CreatePressure(
            fem, facetSet, 0.5e6);
        load->SetName("Pressure");
        fem->AddModelLoad(load);
    }
}

//-----------------------------------------------------------------------------
// Example 8: Initialize and update loads
//-----------------------------------------------------------------------------
void Example8_InitializeAndUpdate(FEModel* fem)
{
    // Initialize all loads
    int nloads = fem->ModelLoads();
    for (int i = 0; i < nloads; ++i)
    {
        RgLoad* load = dynamic_cast<RgLoad*>(fem->ModelLoad(i));
        if (!load->Init())
        {
            std::cerr << "Failed to initialize load: " << load->GetName() << "\n";
        }
    }
    
    // Activate loads
    for (int i = 0; i < nloads; ++i)
    {
        RgLoad* load = dynamic_cast<RgLoad*>(fem->ModelLoad(i));
        load->Activate();
    }
    
    // In time stepping loop:
    double time = 0.0;
    double dt = 0.1;
    
    for (int step = 0; step < 10; ++step)
    {
        time += dt;
        
        // Evaluate load controllers
        fem->EvaluateLoadControllers(time);
        
        // Update all loads
        for (int i = 0; i < nloads; ++i)
        {
            RgLoad* load = dynamic_cast<RgLoad*>(fem->ModelLoad(i));
            if (load->IsActive())
            {
                load->Update();
            }
        }
        
        // Solve system...
    }
}

//-----------------------------------------------------------------------------
// Example 9: Complete workflow with Abaqus input
//-----------------------------------------------------------------------------
void Example9_CompleteWorkflow(const char* inputFile)
{
    FEModel* fem = new FEModel();
    
    // 1. Read mesh and sets (assumed done)
    // ...
    
    // 2. Parse loads from Abaqus file
    std::vector<std::string> lines;
    // ... read file ...
    
    AbaqusLoadParser parser(fem);
    for (size_t i = 0; i < lines.size(); ++i)
    {
        if (lines[i].find("*Cload") == 0)
        {
            int consumed = parser.ParseCload(lines, i);
            i += consumed - 1;
        }
        else if (lines[i].find("*Dload") == 0)
        {
            int consumed = parser.ParseDload(lines, i);
            i += consumed - 1;
        }
    }
    
    parser.CreateLoads();
    
    // 3. Optionally add load curves
    int nloads = fem->ModelLoads();
    for (int i = 0; i < nloads; ++i)
    {
        RgLoad* load = dynamic_cast<RgLoad*>(fem->ModelLoad(i));
        
        // Attach load curve if needed
        // load->SetLoadController(someLoadCurve);
    }
    
    // 4. Initialize
    fem->Init();
    
    // 5. Solve
    fem->Solve();
    
    delete fem;
}

//-----------------------------------------------------------------------------
// Example 10: Cyclic loading for fatigue
//-----------------------------------------------------------------------------
void Example10_CyclicLoading(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    FENodeSet* nodeSet = mesh.FindNodeSet("FatiguePoint");
    
    // Sinusoidal load: F = F_mean + F_amp * sin(2*pi*f*t)
    RgSineController* lc = new RgSineController();
    lc->SetParameters(
        500.0,  // Amplitude: 500N
        2.0,    // Frequency: 2 Hz
        0.0,    // Phase
        1000.0  // Mean load: 1000N
    );
    // Result: Force varies from 500N to 1500N
    
    fem->AddLoadController(lc);
    
    RgNodalLoad* load = RgLoadFactory::CreateNodalLoad(
        fem, nodeSet, 2, 1.0);  // Scale is 1, all in controller
    load->SetLoadController(lc);
    load->SetName("CyclicFatigueLoad");
    fem->AddModelLoad(load);
}

//-----------------------------------------------------------------------------
int main()
{
    FEModel fem;
    
    // Run examples
    // Example1_NodalForces(&fem);
    // Example2_SurfaceLoads(&fem);
    // Example5_ParseAbaqusLoads(&fem);
    
    return 0;
}
