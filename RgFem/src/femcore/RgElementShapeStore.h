#pragma once
#include <map>
#include "fecore_enum.h"
#include "fecore_api.h"
#include "elements/RgElemTypeDefine.h"

class RgElementShape;

//-----------------------------------------------------------------------------
//! This class stores the different element shape classes
class FEM_EXPORT RgElementShapeStore
{
public:
    //! destructor
    ~RgElementShapeStore();

    //! return the element shape store instance
    static RgElementShapeStore* GetInstance();

    //! return element shape class
    RgElementShape* GetElementShape(ElementShape eshape);

    //! initialize library
    static void Initialize();

private:
    //! constructor
    RgElementShapeStore() {}
    RgElementShapeStore(const RgElementShapeStore&) {}

    //! Function to register an element shape class
    void RegisterShape(ElementShape shapeType, RgElementShape* pshape);

private:
    std::map<ElementShape, RgElementShape*> m_ShapeMap; //!< map of element shapes by shape enum
    static RgElementShapeStore* m_pThis;
};