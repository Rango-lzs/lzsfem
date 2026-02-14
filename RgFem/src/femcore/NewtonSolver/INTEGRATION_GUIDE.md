# Solver Refactoring - Complete Implementation

## 📦 File List

### Header Files (.h)
1. **FESolver_refactored.h** - Base solver class
2. **FENewtonSolver_refactored.h** - Newton iteration base class
3. **FEStaticSolver.h** - Static analysis solver
4. **FEImplicitDynamicSolver.h** - Implicit dynamic solver
5. **FEExplicitDynamicSolver.h** - Explicit dynamic solver

### Implementation Files (.cpp)
1. **FESolver_refactored.cpp** - Base solver implementation
2. **FENewtonSolver_refactored.cpp** - Newton solver implementation
3. **FEStaticSolver.cpp** - Static solver implementation
4. **FEImplicitDynamicSolver.cpp** - Implicit dynamic implementation
5. **FEExplicitDynamicSolver.cpp** - Explicit dynamic implementation

### Documentation and Examples
1. **SOLVER_REFACTORING.md** - Design documentation
2. **SolverUsageExample.cpp** - Complete usage examples

## 🔧 Integration Steps

### Step 1: File Organization
```
YourProject/
├── femcore/
│   └── Solver/
│       ├── FESolver_refactored.h
│       ├── FESolver_refactored.cpp
│       ├── FENewtonSolver_refactored.h
│       ├── FENewtonSolver_refactored.cpp
│       ├── FEStaticSolver.h
│       ├── FEStaticSolver.cpp
│       ├── FEImplicitDynamicSolver.h
│       ├── FEImplicitDynamicSolver.cpp
│       ├── FEExplicitDynamicSolver.h
│       └── FEExplicitDynamicSolver.cpp
```

### Step 2: Update CMakeLists.txt
```cmake
# Add solver files
set(SOLVER_SOURCES
    femcore/Solver/FESolver_refactored.cpp
    femcore/Solver/FENewtonSolver_refactored.cpp
    femcore/Solver/FEStaticSolver.cpp
    femcore/Solver/FEImplicitDynamicSolver.cpp
    femcore/Solver/FEExplicitDynamicSolver.cpp
)

set(SOLVER_HEADERS
    femcore/Solver/FESolver_refactored.h
    femcore/Solver/FENewtonSolver_refactored.h
    femcore/Solver/FEStaticSolver.h
    femcore/Solver/FEImplicitDynamicSolver.h
    femcore/Solver/FEExplicitDynamicSolver.h
)

add_library(FECore ${SOLVER_SOURCES} ${SOLVER_HEADERS} ...)
```

### Step 3: Dependencies

Each solver requires:
```cpp
// Core dependencies
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include "femcore/FENode.h"
#include "femcore/FELinearSystem.h"
#include "femcore/FEGlobalMatrix.h"
#include "femcore/FEGlobalVector.h"
#include "femcore/LinearSolver.h"
#include "femcore/Domain/RgDomain.h"
#include "logger/log.h"

// For Newton solvers
#include "femcore/Solver/FELineSearch.h"
#include "femcore/Solver/FENewtonStrategy.h"
```

### Step 4: Migration from Old Code

#### Replace FESolidSolver2 for Static Analysis
```cpp
// Old code
FESolidSolver2* solver = new FESolidSolver2();

// New code
FEStaticSolver* solver = new FEStaticSolver();
solver->m_Rtol = 1e-5;
solver->m_Dtol = 1e-3;
```

#### Replace FESolidSolver2 for Dynamic Analysis
```cpp
// Old code
FESolidSolver2* solver = new FESolidSolver2();
solver->m_rhoi = 0.8;

// New code
FEImplicitDynamicSolver* solver = new FEImplicitDynamicSolver();
solver->m_timeScheme = GENERALIZED_ALPHA;
solver->m_rhoi = 0.8;
```

## 🔍 Key Implementation Notes

### 1. FESolver Base Class
**Responsibilities:**
- Equation numbering (STAGGERED only)
- DOF management
- Matrix profile building
- Solution vector management

**Key Methods:**
```cpp
bool InitEquations();           // Assign equation numbers
void BuildMatrixProfile();      // Build sparse matrix structure
int GetPartitionSize();         // For parallel solving
FENodalDofInfo GetDOFInfoFromEquation(int ieq);
```

### 2. FENewtonSolver
**Responsibilities:**
- Newton-Raphson iteration
- Convergence checking
- Line search
- QN strategies (BFGS, Broyden, JFNK)

**Key Methods:**
```cpp
bool Quasin();                  // Main iteration loop
bool QNInit();                  // Initialize iteration
bool QNUpdate();                // Update in iteration
double QNSolve();               // Solve linear system
bool CheckConvergence();        // Check tolerances
```

### 3. FEStaticSolver
**Equilibrium Equation:**
```
K * u = F_ext - F_int + F_contact + F_constraint
```

**Key Features:**
- Pure stiffness-based solution
- No inertia or damping
- Contact and constraints
- Arc-length method ready

### 4. FEImplicitDynamicSolver
**Dynamic Equilibrium:**
```
M * a + C * v + K * u = F
```

**Time Integration Schemes:**
1. **Newmark (β=0.25, γ=0.5)**
   - Average acceleration
   - Unconditionally stable
   
2. **Generalized-Alpha**
   - High-frequency damping
   - Controlled by ρ∞
   
3. **HHT-Alpha**
   - Similar to Gen-Alpha
   - α parameter

**Key Formulas:**
```cpp
// Acceleration
a(n+1) = 1/(β*dt²) * (u(n+1) - u(n)) 
       - 1/(β*dt) * v(n) 
       - (1/2β - 1) * a(n)

// Velocity
v(n+1) = v(n) + dt * ((1-γ)*a(n) + γ*a(n+1))

// Effective stiffness
K_eff = K + γ/(β*dt)*C + 1/(β*dt²)*M
```

### 5. FEExplicitDynamicSolver
**Central Difference:**
```cpp
// Accelerations
a(n) = M^-1 * (F(n) - F_int(n))

// Velocities
v(n+1/2) = v(n-1/2) + a(n) * dt

// Displacements
u(n+1) = u(n) + v(n+1/2) * dt
```

**Stability (CFL Condition):**
```
dt ≤ L_min / c
```
where:
- L_min: minimum element size
- c: wave speed = sqrt(E/ρ)

**Key Features:**
- No matrix inversion
- Lumped or consistent mass
- Automatic time step control
- CFL stability checking

## ⚠️ Important Implementation Details

### 1. Mass Matrix Assembly
For implicit dynamics, you need to implement:
```cpp
void RgDomain::MassMatrix(FELinearSystem& LS);
```

For explicit dynamics with lumped mass:
```cpp
void RgDomain::LumpedMassMatrix(std::vector<double>& M);
```

### 2. Force Calculations
All domains must implement:
```cpp
void RgDomain::InternalForces(FEGlobalVector& R);
void RgDomain::StiffnessMatrix(FELinearSystem& LS);
```

### 3. Contact Handling
Contact interfaces need:
```cpp
void FESurfacePairConstraint::LoadVector(FEGlobalVector& R);
void FESurfacePairConstraint::StiffnessMatrix(FELinearSystem& LS);
```

### 4. DOF Structure
Make sure you have these DOFs defined:
```cpp
// Displacement
DOFS::AddDOF("x");
DOFS::AddDOF("y");
DOFS::AddDOF("z");

// Velocity (for dynamics)
DOFS::AddDOF("vx");
DOFS::AddDOF("vy");
DOFS::AddDOF("vz");

// Acceleration (for explicit dynamics - optional)
DOFS::AddDOF("ax");
DOFS::AddDOF("ay");
DOFS::AddDOF("az");
```

## 🧪 Testing Recommendations

### 1. Static Solver Test
```cpp
// Cantilever beam under end load
// Expected: tip deflection = (F*L³)/(3*E*I)
```

### 2. Implicit Dynamic Test
```cpp
// Free vibration of cantilever
// Expected: frequency = (λn²/2π) * sqrt(E*I/(ρ*A*L⁴))
// λ1 = 1.875, λ2 = 4.694
```

### 3. Explicit Dynamic Test
```cpp
// Wave propagation in rod
// Expected: wave speed c = sqrt(E/ρ)
// Verify CFL stability: dt < L_min/c
```

## 🚀 Performance Optimization

### For Large Problems
1. **Use sparse matrix storage** (CSR, CSC)
2. **Enable OpenMP** for element loops
3. **Use iterative solvers** (CG, GMRES) for large systems
4. **Lumped mass** for explicit dynamics

### For Contact Problems
1. **Use augmented Lagrangian** for implicit
2. **Use penalty method** for explicit
3. **Enable adaptive time stepping**

### For Nonlinear Problems
1. **Use line search** to improve robustness
2. **Use BFGS** for quasi-Newton updates
3. **Enable auto-reformation** when diverging

## 📊 Comparison with Old Code

| Feature | Old FESolidSolver2 | New Architecture |
|---------|-------------------|------------------|
| Static | Mixed in one class | FEStaticSolver |
| Implicit Dynamic | Mixed in one class | FEImplicitDynamicSolver |
| Explicit Dynamic | Not supported | FEExplicitDynamicSolver |
| Equation Scheme | STAGGERED/BLOCK | STAGGERED only |
| Time Integration | Newmark, Gen-Alpha | Newmark, Gen-Alpha, HHT |
| Code Clarity | ~2000 lines/class | ~400-600 lines/class |
| Extensibility | Difficult | Easy |

## 🔄 Next Steps After Integration

1. **Test each solver** with simple problems
2. **Implement domain methods** (mass matrix, etc.)
3. **Implement contact methods**
4. **Add output capabilities**
5. **Performance profiling**
6. **Add more time integration schemes** (if needed)
7. **Implement adaptive time stepping**

## 📝 TODO Items in Code

Search for `TODO:` comments in the implementation:

1. **FEStaticSolver.cpp**
   - Contact force assembly
   - Constraint force assembly

2. **FEImplicitDynamicSolver.cpp**
   - Mass matrix assembly
   - Proper inertial force calculation
   - Damping matrix assembly

3. **FEExplicitDynamicSolver.cpp**
   - Element size calculation
   - Wave speed from material
   - Lumped mass assembly from elements

## 📚 References

1. **Newmark Method**: Newmark, N.M. (1959). "A method of computation for structural dynamics"
2. **Generalized-Alpha**: Chung, J., Hulbert, G.M. (1993). "A time integration algorithm for structural dynamics"
3. **HHT-Alpha**: Hilber, H.M., Hughes, T.J.R., Taylor, R.L. (1977). "Improved numerical dissipation"
4. **Explicit Dynamics**: Belytschko, T., et al. (2000). "Nonlinear Finite Elements for Continua and Structures"

---

**Last Updated**: February 2025
**Version**: 1.0
**Status**: Ready for integration and testing
