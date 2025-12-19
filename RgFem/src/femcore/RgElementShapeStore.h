#pragma once
#include <map>
#include "fecore_enum.h"
#include "fecore_api.h"

class FEElementShape;

//-----------------------------------------------------------------------------
//! This class stores the different element shape classes
class FECORE_API RgElementShapeStore
{
public:
    //! destructor
    ~RgElementShapeStore();

    //! return the element shape store instance
    static RgElementShapeStore* GetInstance();

    //! return element shape class
    FEElementShape* GetElementShapeClass(FE_Element_Shape eshape);

    //! initialize library
    static void Initialize();

private:
    //! constructor
    RgElementShapeStore() {}
    RgElementShapeStore(const RgElementShapeStore&) {}

    //! Function to register an element shape class
    void RegisterShape(FE_Element_Shape shapeType, FEElementShape* pshape);

private:
    std::map<FE_Element_Shape, FEElementShape*> m_ShapeMap; //!< map of element shapes by shape enum
    static RgElementShapeStore* m_pThis;
};