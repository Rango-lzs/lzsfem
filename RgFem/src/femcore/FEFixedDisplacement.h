#pragma once
#include "femcore/FEFixedBC.h"

class FEFixedDisplacement : public FEFixedBC
{
    DECLARE_META_CLASS(FEFixedDisplacement, FEFixedBC);

public:
    FEFixedDisplacement();

    bool Init() override;

private:
    bool m_dofx;
    bool m_dofy;
    bool m_dofz;

    DECLARE_PARAM_LIST();
};
