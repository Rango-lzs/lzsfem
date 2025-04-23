#pragma once
#include "femcore/FEFixedBC.h"

class FEFixedRotation : public FEFixedBC
{
public:
    FEFixedRotation(FEModel* fem);

    bool Init() override;

private:
    bool m_dof[3];

    DECLARE_PARAM_LIST();
};
