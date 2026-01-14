#pragma once
#include "elements/RgElement/RgElement.h"



class FEM_EXPORT RgTrussElement : public RgElement
{
public:
    RgTrussElement();

    RgTrussElement(const RgTrussElement& el);

    RgTrussElement& operator=(const RgTrussElement& el);

    void Serialize(DumpStream& ar) override;

public:
    double m_a0;   // cross-sectional area
    double m_lam;  // current stretch ratio
    double m_tau;  // Kirchoff stress
    double m_L0;   // initial length
};
