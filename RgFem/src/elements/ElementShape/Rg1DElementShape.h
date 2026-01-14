#pragma once
#include "elements/ElementShape/RgElementShape.h"
#include "elements/NaturalCoord.h"
#include "../RgElemTypeDefine.h"
#include <vector>

namespace RgFem
{
    class NaturalCoord;
}

//=============================================================================
// Base class for defining element shape classes for (3D) solid elements


class FEM_EXPORT Rg1DElementShape : public RgElementShape
{
public:
	Rg1DElementShape(ElementShape shape, int nodes) : RgElementShape(shape, nodes) {}

	//values of shape functions with size N
	virtual std::vector<double> evalH(const NaturalCoord& coord) override;

	//values of shape function derivatives with size 2,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;

	//values of shape function second derivatives with size 3,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};

// 2-node linear 1D element (line element)
class FEM_EXPORT Rg1D2NodeElementShape : public Rg1DElementShape
{
public:
    Rg1D2NodeElementShape() : Rg1DElementShape(ET_LINE2, 2) {}

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};

// 3-node quadratic 1D element
class FEM_EXPORT Rg1D3NodeElementShape : public Rg1DElementShape
{
public:
    Rg1D3NodeElementShape() : Rg1DElementShape(ET_LINE3, 3) {}

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};

// 4-node cubic 1D element
class FEM_EXPORT Rg1D4NodeElementShape : public Rg1DElementShape
{
public:
    Rg1D4NodeElementShape() : Rg1DElementShape(ET_LINE4, 4) {}

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};

// 2-node cubic Hermite 1D element (each node has value and derivative - 4 DOFs total)
class FEM_EXPORT Rg1DHermiteElementShape : public Rg1DElementShape
{
public:
    Rg1DHermiteElementShape() : Rg1DElementShape(ET_LINE4, 4) {} // Same DOF count as 4-node Lagrange

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};

