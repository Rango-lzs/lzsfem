#pragma once
#include <map>
#include "fecore_enum.h"
#include "fecore_api.h"
#include "elements/RgElemTypeDefine.h"

// Forward declarations
class FEElement;
class RgElementTraits;

//-----------------------------------------------------------------------------
//! This class stores the different element traits classes

class FEM_EXPORT RgElementTraitsStore
{
public:
    //! destructor
    ~RgElementTraitsStore();

    //! return the element traits store instance
    static RgElementTraitsStore* GetInstance();

    //! Assign a traits class to an element
    static void SetElementTraits(FEElement& el, int id);

    //! return element traits data
    RgElementTraits* GetElementTraits(int ntype);

    //! checks if the element spec is valid
    bool IsValid(const FE_Element_Spec& c);

    //! get the element spec from the type
    FE_Element_Spec GetElementSpecFromType(ElementType elemType);

    //! initialize library
    static void Initialize();

private:
    //! constructor
    RgElementTraitsStore() {}
    RgElementTraitsStore(const RgElementTraitsStore&) {}

    //! Function to register an element traits class
    void RegisterTraits(ElementType traitType, RgElementTraits* ptrait);

private:
    std::map<ElementType, RgElementTraits*> m_TraitsMap; //!< map of element traits by type enum
    static RgElementTraitsStore* m_pThis;
};