set(ELEMENT_TRAITS_HEADER
    elements/ElementTraits/RgElementTraits.h
    elements/ElementTraits/RgSolidElementTraits.h
    #elements/ElementTraits/Rg2DElementTraits.h
    #elements/ElementTraits/RgShellElementTraits.h
    #elements/ElementTraits/RgSurfaceElementTraits.h
    #elements/ElementTraits/RgTrussElementTraits.h
    PARENT_SCOPE
)

set(ELEMENT_TRAITS_SOURCE
    elements/ElementTraits/RgElementTraits.cpp
    elements/ElementTraits/RgSolidElementTraits.cpp
    #elements/ElementTraits/RgShellElementTraits.cpp
    #PARENT_SCOPE
)

set(ELEMENT_HEADERS
    # Core element classes (base classes)
    elements/RgElement/RgElement.h
    elements/RgElement/RgSolidElement.h
    elements/RgElement/RgLinearSolid2dElement.h
    elements/RgElement/RgLinearSolid3dElement.h
    elements/RgElement/RgNLSolid2dElement.h
    elements/RgElement/RgNLSolid3dElement.h
    elements/RgElement/RgSolid3dElement.h
    elements/RgElement/RgLinearShellElement.h
    #elements/RgElement/RgShellElement.h
    #elements/RgElement/RgNLShellElement.h
    #elements/RgElement/RgLinearBeamElement.h
    #elements/RgElement/RgNLBeamElement.h
    #elements/RgElement/RgStructureElement.h
    #elements/RgElement/RgSurfaceElement.h
    #elements/RgElement/RgTrussElement.h

    # Specific 3D solid elements
    elements/RgElement/RgHex8Element.h
    elements/RgElement/RgHex20Element.h
    elements/RgElement/RgHex8GeomNLElement.h
    elements/RgElement/RgTet4Element.h
    elements/RgElement/RgTet10Element.h
    elements/RgElement/RgPyramid5Element.h
    elements/RgElement/RgWedge6Element.h

    # 2D elements
    #elements/RgElement/RgQuad4Element.h
    #elements/RgElement/RgQuad8Element.h
    #elements/RgElement/RgTri3Element.h
    #elements/RgElement/RgTri6Element.h

    # Shell elements
    #elements/RgElement/RgShell3Element.h
    #elements/RgElement/RgShell4Element.h

    # Beam elements
    #elements/RgElement/RgBeamElement.h
    #elements/RgElement/RgBeam2dElement.h
    #elements/RgElement/RgBeam2dGeomNLElement.h
    #elements/RgElement/RgBeam3dElement.h
    #elements/RgElement/RgBeam3dGeomNLElement.h

    # Nonlinear elements
    #elements/RgElement/RgNLQuad4Element.h
    #elements/RgElement/RgNLTri3Element.h

    # Shape classes (base to specific shapes)
    elements/ElementShape/RgElementShape.h
    elements/ElementShape/RgSolidElementShape.h
    #elements/ElementShape/Rg2DElementShape.h
    #elements/ElementShape/Rg1DElementShape.h
    elements/ElementShape/RgHex8Shape.h
    elements/ElementShape/RgHex20Shape.h
    elements/ElementShape/RgHex27Shape.h
    elements/ElementShape/RgTet4Shape.h
    elements/ElementShape/RgTet10Shape.h
    elements/ElementShape/RgPenta6Shape.h
    elements/ElementShape/RgPenta15Shape.h
    elements/ElementShape/RgPyra5Shape.h
    elements/ElementShape/RgPyra13Shape.h

    # Traits classes
    elements/ElementTraits/RgElementTraits.h
    elements/ElementTraits/RgSolidElementTraits.h
    #elements/ElementTraits/Rg2DElementTraits.h
    #elements/ElementTraits/RgShellElementTraits.h
    #elements/ElementTraits/RgSurfaceElementTraits.h
    #elements/ElementTraits/RgTrussElementTraits.h

    # Utility headers
    elements/NaturalCoord.h
    elements/RgGaussPoint.h
    elements/RgElemElemList.h
    elements/RgElementList.h
    elements/RgElementSet.h
    elements/RgElementState.h
    elements/RgElemTypeDefine.h
    ${ELEMENT_TRAITS_HEADER}
    PARENT_SCOPE
)

set(ELEMENT_SOURCES
    # Core element source files
    elements/RgElement/RgElement.cpp
    elements/RgElement/RgSolidElement.cpp
    elements/RgElement/RgLinearSolid2dElement.cpp
    elements/RgElement/RgLinearSolid3dElement.cpp
    elements/RgElement/RgNLSolid2dElement.cpp
    elements/RgElement/RgNLSolid3dElement.cpp
    elements/RgElement/RgSolid3dElement.cpp
    #elements/RgElement/RgLinearShellElement.cpp
    #elements/RgElement/RgShellElement.cpp
    #elements/RgElement/RgNLShellElement.cpp
    #elements/RgElement/RgLinearBeamElement.cpp
    #elements/RgElement/RgNLBeamElement.cpp
    #elements/RgElement/RgStructureElement.cpp
    #elements/RgElement/RgSurfaceElement.cpp
    #elements/RgElement/RgTrussElement.cpp

    # Specific 3D solid elements
    elements/RgElement/RgHex8Element.cpp
    elements/RgElement/RgHex20Element.cpp
    elements/RgElement/RgHex8GeomNLElement.cpp
    elements/RgElement/RgTet4Element.cpp
    elements/RgElement/RgTet10Element.cpp
    elements/RgElement/RgPyramid5Element.cpp
    elements/RgElement/RgWedge6Element.cpp

    # 2D elements
    #elements/RgElement/RgQuad4Element.cpp
    #elements/RgElement/RgQuad8Element.cpp
    #elements/RgElement/RgTri3Element.cpp
    #elements/RgElement/RgTri6Element.cpp

    # Shell elements
    #elements/RgElement/RgShell3Element.cpp
    #elements/RgElement/RgShell4Element.cpp

    # Beam elements
    #elements/RgElement/RgBeamElement.cpp
    #elements/RgElement/RgBeam2dElement.cpp
    #elements/RgElement/RgBeam2dGeomNLElement.cpp
    #elements/RgElement/RgBeam3dElement.cpp
    #elements/RgElement/RgBeam3dGeomNLElement.cpp

    # Nonlinear elements
    #elements/RgElement/RgNLQuad4Element.cpp
    #elements/RgElement/RgNLTri3Element.cpp

    # Shape source files
    elements/ElementShape/RgElementShape.cpp
    elements/ElementShape/RgSolidElementShape.cpp
    #elements/ElementShape/Rg2DElementShape.cpp
    #elements/ElementShape/Rg1DElementShape.cpp
    elements/ElementShape/RgHex8Shape.cpp
    elements/ElementShape/RgHex20Shape.cpp
    elements/ElementShape/RgHex27Shape.cpp
    elements/ElementShape/RgTet4Shape.cpp
    elements/ElementShape/RgTet10Shape.cpp
    elements/ElementShape/RgPenta6Shape.cpp
    elements/ElementShape/RgPenta15Shape.cpp
    elements/ElementShape/RgPyra5Shape.cpp
    elements/ElementShape/RgPyra13Shape.cpp

    # Traits source files
    elements/ElementTraits/RgElementTraits.cpp
    #elements/ElementTraits/RgShellElementTraits.cpp
    #elements/ElementTraits/Rg2DElementTraits.cpp
    #elements/ElementTraits/RgSurfaceElementTraits.cpp

    # Utility source files
    elements/NaturalCoord.cpp
    elements/RgElemElemList.cpp
    elements/RgElementList.cpp
    elements/RgElementSet.cpp
    elements/RgElementState.cpp
    elements/RgGaussPoint.cpp
    ${ELEMENT_TRAITS_SOURCE}
    PARENT_SCOPE
)