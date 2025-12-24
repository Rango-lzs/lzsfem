#include "quadrature/RgGaussPoint.h"
#include <iostream>
#include <vector>

using namespace RgFem;

int main() {
    // Test 1D Gauss point
    RgGaussPoint gp1d(0.0, 2.0);
    std::cout << "1D Gauss point:" << std::endl;
    std::cout << "  Dimension: " << gp1d.getDimension() << std::endl;
    std::cout << "  Coordinate: " << gp1d.getR() << std::endl;
    std::cout << "  Weight: " << gp1d.getWeight() << std::endl << std::endl;

    // Test 2D Gauss point
    RgGaussPoint gp2d(-0.577, 0.577, 1.0);
    std::cout << "2D Gauss point:" << std::endl;
    std::cout << "  Dimension: " << gp2d.getDimension() << std::endl;
    std::cout << "  Coordinates: (" << gp2d.getR() << ", " << gp2d.getS() << ")" << std::endl;
    std::cout << "  Weight: " << gp2d.getWeight() << std::endl << std::endl;

    // Test 3D Gauss point
    RgGaussPoint gp3d(0.577, 0.577, 0.577, 1.0);
    std::cout << "3D Gauss point:" << std::endl;
    std::cout << "  Dimension: " << gp3d.getDimension() << std::endl;
    std::cout << "  Coordinates: (" << gp3d.getR() << ", " << gp3d.getS() << ", " << gp3d.getT() << ")" << std::endl;
    std::cout << "  Weight: " << gp3d.getWeight() << std::endl << std::endl;

    // Test general coordinate access
    std::cout << "Accessing coordinates by index:" << std::endl;
    for (size_t i = 0; i < gp3d.getCoordinates().size(); ++i) {
        std::cout << "  Coordinate " << i << ": " << gp3d.getCoordinate(static_cast<int>(i)) << std::endl;
    }

    // Test differentiation between 2D and 3D
    std::cout << "\nDifferentiation examples:" << std::endl;
    if (gp2d.getDimension() == 2) {
        std::cout << "gp2d is a 2D Gauss point" << std::endl;
    }
    
    if (gp3d.getDimension() == 3) {
        std::cout << "gp3d is a 3D Gauss point" << std::endl;
    }
    
    if (gp1d.getDimension() == 1) {
        std::cout << "gp1d is a 1D Gauss point" << std::endl;
    }

    return 0;
}