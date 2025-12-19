# RgGaussPoint Class Documentation

## Overview

The `RgGaussPoint` class represents a Gauss integration point used in the finite element method for numerical integration. It stores the natural coordinates and integration weight of a Gauss point, which can be used in 1D, 2D, or 3D integration schemes.

## Design Principles

The class follows these design principles:

1. **Simplicity**: The class provides a clean interface for storing and accessing Gauss point data.
2. **Flexibility**: Can handle 1D, 2D, and 3D integration schemes through a variable number of coordinates.
3. **Efficiency**: Uses std::vector for storing coordinates to allow for flexible sizing.
4. **Consistency**: Follows the naming conventions and coding style of the RgFem project.

## Class Interface

### Constructors

- `RgGaussPoint()` - Default constructor, initializes to 3D with zero coordinates and zero weight
- `RgGaussPoint(const std::vector<double>& coordinates, double weight)` - Constructor with coordinates vector and weight
- `RgGaussPoint(double r, double s, double t, double weight)` - Constructor for 3D Gauss point
- `RgGaussPoint(double r, double s, double weight)` - Constructor for 2D Gauss point
- `RgGaussPoint(double r, double weight)` - Constructor for 1D Gauss point

### Dimension Methods

- `int getDimension() const` - Get the dimension of the Gauss point (1, 2, or 3)

### Accessors

- `const std::vector<double>& getCoordinates() const` - Get all coordinates
- `double getR() const` - Get the r-coordinate (first natural coordinate)
- `double getS() const` - Get the s-coordinate (second natural coordinate)
- `double getT() const` - Get the t-coordinate (third natural coordinate)
- `double getCoordinate(int index) const` - Get coordinate by index
- `double getWeight() const` - Get the integration weight

### Mutators

- `void setCoordinates(const std::vector<double>& coordinates)` - Set coordinates from vector
- `void setCoordinates(double r, double s, double t)` - Set 3D coordinates
- `void setCoordinates(double r, double s)` - Set 2D coordinates
- `void setWeight(double weight)` - Set the integration weight

## Distinguishing Between 2D and 3D Gauss Points

In finite element analysis, distinguishing between 2D and 3D Gauss points is important for proper integration:

1. **By Dimension**: 
   - 1D Gauss points have 1 coordinate
   - 2D Gauss points have 2 coordinates (r, s)
   - 3D Gauss points have 3 coordinates (r, s, t)

2. **Usage Context**:
   - 2D Gauss points are used for planar elements (triangles, quads)
   - 3D Gauss points are used for solid elements (tetrahedra, hexahedra)

Example of distinguishing between dimensions:
```cpp
RgGaussPoint gp2d(0.577, 0.577, 1.0);     // 2D point
RgGaussPoint gp3d(0.577, 0.577, 0.577, 1.0); // 3D point

if (gp2d.getDimension() == 2) {
    // Process as 2D point
}

if (gp3d.getDimension() == 3) {
    // Process as 3D point
}
```

## Usage Example

```cpp
#include "quadrature/RgGaussPoint.h"

using namespace RgFem;

// Create a Gauss point for hexahedral elements (8-point Gauss integration)
RgGaussPoint gp3d(0.577, 0.577, 0.577, 1.0);

// Create a Gauss point for quadrilateral elements
RgGaussPoint gp2d(0.577, 0.577, 1.0);

// Access coordinates
double r = gp3d.getR();
double s = gp3d.getS();
double t = gp3d.getT();

// Access weight
double weight = gp3d.getWeight();

// Check dimension
int dim = gp3d.getDimension(); // Returns 3

// Use in integration
double integrand_value = /* some function evaluation at (r,s,t) */;
double weighted_value = integrand_value * weight;
```

## Integration with Existing Code

The `RgGaussPoint` class is designed to work seamlessly with existing finite element classes in the RgFem framework. It can be used to:

1. Store Gauss point data in element classes
2. Simplify integration loops in stiffness matrix calculations
3. Provide a uniform interface for accessing Gauss point information

In existing element classes like `RgHex8Element`, instead of separate vectors for `m_gaussR`, `m_gaussS`, `m_gaussT`, and `m_gaussW`, a vector of `RgGaussPoint` objects could be used:

```cpp
class RgHex8Element : public RgSolid3dElement {
    // Instead of:
    // std::vector<double> m_gaussR, m_gaussS, m_gaussT, m_gaussW;
    
    // Use:
    std::vector<RgGaussPoint> m_gaussPoints;
};
```

This would simplify the code and make it more maintainable.

## Future Extensions

Possible extensions to consider:

1. Adding methods to transform natural coordinates to global coordinates
2. Including Jacobian determinant information
3. Adding support for different integration schemes
4. Including material point references for constitutive modeling at Gauss points