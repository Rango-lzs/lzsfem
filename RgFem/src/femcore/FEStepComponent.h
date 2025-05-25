#pragma once
#include "femcore/FEModelComponent.h"

//-----------------------------------------------------------------------------
// A Step component is a model component that can be assigned to a step.
// It adds a mechanism for activating and deactivating the component.
class FEM_EXPORT FEStepComponent : public FEModelComponent
{
public:
    FEStepComponent(FEModel* fem);
    FEStepComponent();

    //-----------------------------------------------------------------------------------
    //! This function checks if the component is active in the current step.
    bool IsActive() const;

    //-----------------------------------------------------------------------------------
    //! Activate the component.
    //! This function is called during the step initialization, right before the step is solved.
    //! This function can be used to initialize any data that could depend on the model state.
    //! Data allocation and initialization of data that does not depend on the model state should
    //! be done in Init().
    virtual void Activate();

    //-----------------------------------------------------------------------------------
    //! Deactivate the component
    virtual void Deactivate();

public:
    //! serialization
    void Serialize(DumpStream& ar);

private:
    bool m_bactive;  //!< flag indicating whether the component is active
};
