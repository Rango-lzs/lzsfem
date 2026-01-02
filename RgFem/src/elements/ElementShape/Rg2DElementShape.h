#pragma once
#include "elements/ElementShape/RgElementShape.h"
#include "elements/NaturalCoord.h"
#include "../RgElemTypeDefine.h"
#include <vector>
#include "femcore/fecore_api.h"

namespace RgFem
{
    class NaturalCoord;
}

//=============================================================================
// Base class for defining element shape classes for (3D) solid elements
namespace RgFem {

class FEM_EXPORT Rg2DElementShape : public RgElementShape
{
public:
	Rg2DElementShape(ElementShape shape, int nodes) : RgElementShape(shape, nodes) {}

	//values of shape functions with size N
	virtual std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;

	//values of shape function derivatives with size 2,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;

	//values of shape function second derivatives with size 3,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;

    ElementShape shapeType() const
    {
        return mShpType;
    }

    int nodes() const
    {
        return mNodes;
    }

private:
    ElementShape mShpType;
    int mNodes;
};

//=============================================================================
// Class for QUAD4 elements
class FEQuad4 : public Rg2DElementShape
{
public:
	FEQuad4() : Rg2DElementShape(ET_QUAD4, 4) {}

	//values of shape functions with size N
	virtual std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;

	//values of shape function derivatives with size 2,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;

	//values of shape function second derivatives with size 3,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;
};

//=============================================================================
// Class for QUAD8 elements
class FEQuad8 : public Rg2DElementShape
{
public:
	FEQuad8() : Rg2DElementShape(ET_QUAD8, 8) {}

	//values of shape functions with size N
	virtual std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;

	//values of shape function derivatives with size 2,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;

	//values of shape function second derivatives with size 3,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;
};

//=============================================================================
// Class for QUAD9 elements
class FEQuad9 : public Rg2DElementShape
{
public:
	FEQuad9() : Rg2DElementShape(ET_QUAD9, 9) {}

	//values of shape functions with size N
	virtual std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;

	//values of shape function derivatives with size 2,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;

	//values of shape function second derivatives with size 3,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;
};

//=============================================================================
// Class for TRI3 elements
class FETri3 : public Rg2DElementShape
{
public:
	FETri3() : Rg2DElementShape(ET_TRI3, 3) {}

	//values of shape functions with size N
	virtual std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;

	//values of shape function derivatives with size 2,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;

	//values of shape function second derivatives with size 3,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;
};

//=============================================================================
// Class for TRI6 elements
class FETri6 : public Rg2DElementShape
{
public:
	FETri6() : Rg2DElementShape(ET_TRI6, 6) {}

	//values of shape functions with size N
	virtual std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;

	//values of shape function derivatives with size 2,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;

	//values of shape function second derivatives with size 3,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;
};

//=============================================================================
// Class for TRI7 elements
class FETri7 : public Rg2DElementShape
{
public:
	FETri7() : Rg2DElementShape(ET_TRI7, 7) {}

	//values of shape functions with size N
	virtual std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;

	//values of shape function derivatives with size 2,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;

	//values of shape function second derivatives with size 3,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;
};

//=============================================================================
// Class for TRI10 elements
class FETri10 : public Rg2DElementShape
{
public:
	FETri10() : Rg2DElementShape(ET_TRI10, 10) {}

	//values of shape functions with size N
	virtual std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;

	//values of shape function derivatives with size 2,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;

	//values of shape function second derivatives with size 3,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;
};

} // namespace RgFem