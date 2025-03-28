#pragma once
#include "femcore/FEObjectBase.h"

class FEModel;
class FEModelParam;

//---------------------------------------------------------------------------------------
// Base class for evaluating model parameters
class FEM_EXPORT FEValuator : public FEObjectBase
{
public:
    FEValuator(FEModel* fem)
        : FEObjectBase(fem)
        , m_param(nullptr)
    {
    }
    virtual ~FEValuator()
    {
    }

    void SetModelParam(FEModelParam* p)
    {
        m_param = p;
    }
    FEModelParam* GetModelParam()
    {
        return m_param;
    }

private:
    FEModelParam* m_param;  //!< the model param that is using this valuator
};
