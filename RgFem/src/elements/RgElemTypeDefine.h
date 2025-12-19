#pragma once

// 单元的类别
enum ElementCategory
{
    FE_ELEM_INVALID_CLASS,
    FE_ELEM_SOLID,
    FE_ELEM_SOLID_2D,
    FE_ELEM_SHELL,
    FE_ELEM_BEAM,
    
    FE_ELEM_TRUSS,
    FE_ELEM_SURFACE,
    FE_ELEM_DISCRETE,
    FE_ELEM_EDGE
};

//-----------------------------------------------------------------------------
// Element shapes:
// This defines the general element shape classes. This classification differs from the
// element types below, in that the latter is defined by a shape and integration rule.
// Do not change the order of these enums!
enum ElementShape
{
    ET_INVALID,
    // 3D elements
    ET_TET4,
    ET_TET10,
    ET_PENTA6,
    ET_PENTA15,
    ET_HEX8,
    ET_HEX20,
    ET_HEX27,
    ET_PYRA5,
    ET_PYRA13,

    // 2.5D elements
    ET_QUAD4,
    ET_QUAD8,
    ET_QUAD9,
    ET_TRI3,
    ET_TRI6,
    ET_TRI7,
    ET_TRI10,

    // line elements
    ET_TRUSS2,
    ET_LINE2,
    ET_DISCRETE,

    FE_ELEM_INVALID_SHAPE = 999
};

//-----------------------------------------------------------------------------
// Element types:
//  Note that these numbers are actually indices into the m_Traits array
//  of the ElementLibrary class so make sure the numbers correspond
//  with the entries into this array
//

enum ElementType
{
    // 3D soid elements
    FE_HEX8G8,
    FE_HEX8RI,
    FE_HEX8G1,
    FE_TET4G1,
    FE_TET4G4,
    FE_TET5G4,
    FE_PENTA6G6,
    FE_TET10G1,
    FE_TET10G4,
    FE_TET10G8,
    FE_TET10GL11,
    FE_TET10G4RI1,
    FE_TET10G8RI4,
    FE_TET15G4,
    FE_TET15G8,
    FE_TET15G11,
    FE_TET15G15,
    FE_TET15G15RI4,
    FE_TET20G15,
    FE_HEX20G8,
    FE_HEX20G27,
    FE_HEX27G27,
    FE_PENTA15G8,
    FE_PENTA15G21,
    FE_PYRA5G8,
    FE_PYRA13G8,

    // 2.5D surface elements
    FE_QUAD4G4,
    FE_QUAD4NI,
    FE_TRI3G1,
    FE_TRI3G3,
    FE_TRI3G7,
    FE_TRI3NI,
    FE_TRI6G3,
    FE_TRI6G4,
    FE_TRI6G7,
    //	FE_TRI6MG7,
    FE_TRI6GL7,
    FE_TRI6NI,
    FE_TRI7G3,
    FE_TRI7G4,
    FE_TRI7G7,
    FE_TRI7GL7,
    FE_TRI10G7,
    FE_TRI10G12,
    FE_QUAD8G9,
    FE_QUAD8NI,
    FE_QUAD9G9,
    FE_QUAD9NI,

    // shell elements
    FE_SHELL_QUAD4G8,
    FE_SHELL_QUAD4G12,
    FE_SHELL_QUAD8G18,
    FE_SHELL_QUAD8G27,
    FE_SHELL_TRI3G6,
    FE_SHELL_TRI3G9,
    FE_SHELL_TRI6G14,
    FE_SHELL_TRI6G21,

    // truss elements
    FE_TRUSS,

    // discrete elements
    FE_DISCRETE,

    // 2D elements
    FE2D_TRI3G1,
    FE2D_TRI6G3,
    FE2D_QUAD4G4,
    FE2D_QUAD8G9,
    FE2D_QUAD9G9,

    // line elements
    FE_LINE2G1,

    // unspecified
    FE_ELEM_INVALID_TYPE = 0xFFFF
};

// Shell formulations
enum SHELL_FORMULATION
{
    NEW_SHELL,
    OLD_SHELL,
    EAS_SHELL,
    ANS_SHELL
};

//! Helper class for creating domain classes.
struct FE_Element_Spec
{
    ElementCategory eclass; //单元的类别
    ElementShape eshape;    //单元形状
    ElementType etype;      //单元类型
    bool m_bthree_field;
    int m_shell_formulation;
    bool m_shell_norm_nodal;

    bool m_but4;
    double m_ut4_alpha;
    bool m_ut4_bdev;

    FE_Element_Spec()
    {
        eclass = FE_ELEM_INVALID_CLASS;
        eshape = FE_ELEM_INVALID_SHAPE;
        etype = FE_ELEM_INVALID_TYPE;
        m_bthree_field = false;
        m_shell_formulation = NEW_SHELL;
        m_shell_norm_nodal = true;
        m_but4 = false;
        m_ut4_alpha = 0.05;
        m_ut4_bdev = false;
    }

    bool operator==(const FE_Element_Spec& s)
    {
        if ((eclass == s.eclass) && (eshape == s.eshape) && (etype == s.etype))
            return true;
        return false;
    }
};
