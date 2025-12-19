#pragma once
#include <map>
#include "fecore_enum.h"
#include "fecore_api.h"

// Forward declarations
class FEElement;
class FEElementTraits;

//-----------------------------------------------------------------------------
//! This class stores the different element traits classes

class FECORE_API RgElementTraitsStore
{
public:
    //! destructor
    ~RgElementTraitsStore();

    //! return the element traits store instance
    static RgElementTraitsStore* GetInstance();

    //! Assign a traits class to an element
    static void SetElementTraits(FEElement& el, int id);

    //! return element traits data
    FEElementTraits* GetElementTraits(int ntype);

    //! return the element class of a given element type
    FE_Element_Class GetElementClass(int ntype);

    //! checks if the element spec is valid
    bool IsValid(const FE_Element_Spec& c);

    //! get the element spec from the type
    FE_Element_Spec GetElementSpecFromType(FE_Element_Type elemType);

    //! initialize library
    static void Initialize();

private:
    //! constructor
    RgElementTraitsStore() {}
    RgElementTraitsStore(const RgElementTraitsStore&) {}

    //! Function to register an element traits class
    void RegisterTraits(FE_Element_Type traitType, FEElementTraits* ptrait);

private:
    std::map<FE_Element_Type, FEElementTraits*> m_TraitsMap; //!< map of element traits by type enum
    static RgElementTraitsStore* m_pThis;
};