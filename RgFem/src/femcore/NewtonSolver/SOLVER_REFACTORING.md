# FE Solver Refactoring

## Overview

This refactoring separates the different analysis types into distinct solver classes with a clear inheritance hierarchy.

## Class Hierarchy

```
FESolver (Abstract base class)
├── FENewtonSolver (Newton-Raphson base for implicit methods)
│   ├── FEStaticSolver (Static/quasi-static analysis)
│   └── FEImplicitDynamicSolver (Implicit dynamic analysis)
└── FEExplicitDynamicSolver (Explicit dynamic analysis)
```

## Key Design Decisions

### 1. Separation of Concerns

**FESolver** (Base class)
- Common interface for all solvers
- Equation numbering (STAGGERED scheme only)
- DOF management
- Matrix profile building
- Basic update mechanisms

**FENewtonSolver** (Newton iteration base)
- Newton-Raphson iteration logic
- Quasi-Newton strategies (BFGS, Broyden, JFNK)
- Line search algorithms
- Convergence checking
- Stiffness matrix reformation logic
- Augmentation handling

**FEStaticSolver** (Static analysis)
- Equilibrium: K * u = F
- No inertia or damping terms
- Suitable for: stress analysis, buckling, nonlinear geometry

**FEImplicitDynamicSolver** (Implicit dynamics)
- Equilibrium: M * a + C * v + K * u = F
- Time integration schemes:
  - Newmark-beta
  - Generalized-alpha
  - HHT-alpha
- Unconditionally stable (for proper parameters)
- Suitable for: structural dynamics, vibration, impact

**FEExplicitDynamicSolver** (Explicit dynamics)
- No matrix inversion required
- Direct calculation: a = M^(-1) * (F - Fint)
- Time integration schemes:
  - Central difference
  - Forward Euler
  - Velocity Verlet
- Conditionally stable (CFL condition)
- Suitable for: wave propagation, high-speed impact, large deformations

### 2. Equation Scheme

Only STAGGERED scheme is retained:
```
| a0, b0, a1, b1, ..., an, bn |
```

This simplifies the code and is the most commonly used scheme in FE codes.

### 3. Time Integration Details

#### Implicit Methods (Newmark-beta)
```
a(n+1) = (1/β/dt²) * (u(n+1) - u(n)) - (1/β/dt) * v(n) - (1/2β - 1) * a(n)
v(n+1) = v(n) + dt * ((1-γ) * a(n) + γ * a(n+1))
```

Parameters:
- β = 0.25, γ = 0.5: Average acceleration (unconditionally stable)
- β = 1/6, γ = 0.5: Linear acceleration
- β = 0.3025, γ = 0.6: Optimal for second-order accuracy

#### Generalized-Alpha
Extends Newmark with numerical damping:
```
M * a(n+1-αm) + C * v(n+1-αf) + K * u(n+1-αf) = F(n+1-αf)
```

where:
- αm = (2*ρ∞ - 1)/(ρ∞ + 1)
- αf = ρ∞/(ρ∞ + 1)
- ρ∞ ∈ [0, 1]: spectral radius at infinity

#### Explicit Methods (Central Difference)
```
M * a(n) = F(n) - Fint(n)
v(n+1/2) = v(n-1/2) + a(n) * dt
u(n+1) = u(n) + v(n+1/2) * dt
```

Stability condition (CFL):
```
dt ≤ dt_critical = L_min / c
```

where:
- L_min: minimum element size
- c: wave speed

### 4. Key Features

#### Static Solver
- Pure equilibrium solution
- Arc-length method support (for snap-through)
- Load control and displacement control
- Contact and constraints

#### Implicit Dynamic Solver
- Multiple time integration schemes
- Automatic time step control
- Predictor-corrector schemes
- Rayleigh damping option
- Contact dynamics

#### Explicit Dynamic Solver
- Lumped or consistent mass matrix
- Automatic critical time step calculation
- CFL stability checking
- Minimal memory requirements
- No matrix factorization

## Migration Guide

### From Old FESolidSolver2 to New System

**For Static Analysis:**
```cpp
// Old
FESolidSolver2* solver = new FESolidSolver2();

// New
FEStaticSolver* solver = new FEStaticSolver();
```

**For Implicit Dynamic Analysis:**
```cpp
// Old
FESolidSolver2* solver = new FESolidSolver2();
solver->m_rhoi = 0.8;
solver->m_alphaf = ...;
solver->m_alpham = ...;

// New
FEImplicitDynamicSolver* solver = new FEImplicitDynamicSolver();
solver->m_timeScheme = GENERALIZED_ALPHA;
solver->m_rhoi = 0.8;
```

**For Explicit Dynamic Analysis:**
```cpp
// New (no direct equivalent in old system)
FEExplicitDynamicSolver* solver = new FEExplicitDynamicSolver();
solver->m_scheme = CENTRAL_DIFFERENCE;
solver->m_lumpedMass = true;
solver->m_autoTimeStep = true;
```

## Implementation Notes

### Mass Matrix Handling

**Lumped Mass (Explicit):**
```cpp
// Diagonal mass matrix
M_ii = Σ(element masses at node i)
```

**Consistent Mass (Implicit):**
```cpp
// Full mass matrix from element integration
M = Σ_e ∫(N^T * ρ * N) dV
```

### Contact Treatment

**Implicit:**
- Penalty method
- Lagrange multipliers
- Augmented Lagrangian
- Consistent tangent stiffness

**Explicit:**
- Penalty method only
- Node-to-surface contact
- Updated each time step

### Convergence Criteria

**Static/Implicit:**
- Displacement norm
- Residual norm
- Energy norm
- Component-wise checking

**Explicit:**
- No iterations required
- Stability checking only

## Performance Considerations

### Static Solver
- Sparse matrix storage
- Iterative solvers for large problems
- Multi-threading for element calculations

### Implicit Dynamic Solver
- Similar to static
- Time step control important
- Predictor improves convergence

### Explicit Dynamic Solver
- Minimal memory (no K matrix)
- Very efficient for short duration
- Excellent parallelization
- Large number of small time steps

## Testing Strategy

1. **Static:**
   - Linear elastic cantilever beam
   - Nonlinear large deformation
   - Contact patch test

2. **Implicit Dynamic:**
   - Free vibration (compare with analytical)
   - Forced vibration with damping
   - Impact with contact

3. **Explicit Dynamic:**
   - Wave propagation (compare with analytical)
   - High-velocity impact
   - CFL stability verification

## Future Extensions

1. **Adaptive time stepping** for implicit methods
2. **Subcycling** for explicit methods
3. **Coupled multi-physics** solvers
4. **Parallel domain decomposition**
5. **GPU acceleration** for explicit methods
