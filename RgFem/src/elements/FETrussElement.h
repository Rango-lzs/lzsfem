#pragma once
#include "elements/FEElement.h"

class FEM_EXPORT FETrussElement : public FEElement
{
public:
    FETrussElement();

    FETrussElement(const FETrussElement& el);

    FETrussElement& operator=(const FETrussElement& el);

    void Serialize(DumpStream& ar) override;

public:
    double m_a0;   // cross-sectional area
    double m_lam;  // current stretch ratio
    double m_tau;  // Kirchoff stress
    double m_L0;   // initial length
};

//-----------------------------------------------------------------------------
//! Discrete element class