#pragma once
#include "femcore/fem_export.h"

#include <list>
#include <vector>

//-----------------------------------------------------------------------------
// Forward declarations.
class FEModel;

//-----------------------------------------------------------------------------
// callback events
#define CB_ALWAYS      0x0FFFFFFF  //!< Call for all reasons
#define CB_INIT        0x00000001  //!< Call after model initialization (i.e. FEModel::Init())
#define CB_STEP_ACTIVE 0x00000002  //!< call after step was activated (i.e.
#define CB_MAJOR_ITERS 0x00000004  //!< Call at the end of each major converged iteration
#define CB_MINOR_ITERS 0x00000008  //!< Call for each minor iteration
#define CB_SOLVED      0x00000010  //!< Call at the end of FEModel::Solve
#define CB_UPDATE_TIME                                                                                                 \
    0x00000020  //!< Call when time is updated and right before time step is solved (in FEAnalysis::Solve)
#define CB_AUGMENT          0x00000040  //!< The model is entering augmentations (called before Augment)
#define CB_STEP_SOLVED      0x00000080  //!< The step was solved
#define CB_MATRIX_REFORM    0x00000100  //!< stiffness Matrix was reformed
#define CB_REMESH           0x00000200  //!< Called after remesh
#define CB_PRE_MATRIX_SOLVE 0x00000400  //!< Called right before Matrix solve
#define CB_RESET            0x00000800  //!< Called after FEModel::Reset
#define CB_MODEL_UPDATE     0x00001000  //!< Called at the end of FEModel::Update
#define CB_TIMESTEP_SOLVED  0x00002000  //!< Called at FEAnalysis::SolveTimeStep after the solver returns.
#define CB_SERIALIZE_SAVE   0x00004000  //!< Called at the end of FEModel::Serialize when saving
#define CB_SERIALIZE_LOAD   0x00008000  //!< Called at the end of FEModel::Serialize when loading
#define CB_USER1            0x01000000  //!< can be used by users

typedef unsigned int FECORE_CB_WHEN;
typedef bool (*FECORE_CB_FNC)(FEModel*, unsigned int, void*);

//-----------------------------------------------------------------------------
// callback structure
struct FECORE_CALLBACK
{
    FECORE_CB_FNC m_pcb;     // pointer to callback function
    void* m_pd;              // pointer to user data
    FECORE_CB_WHEN m_nwhen;  // when to call function
};

//-----------------------------------------------------------------------------
// class that handles callbacks
class FEM_EXPORT CallbackHandler
{
public:
    enum CBInsertPolicy
    {
        CB_ADD_FRONT,
        CB_ADD_END
    };

public:
    CallbackHandler();
    virtual ~CallbackHandler();

    //! set callback function
    void AddCallback(FECORE_CB_FNC pcb, unsigned int nwhen, void* pd,
                     CBInsertPolicy insert = CBInsertPolicy::CB_ADD_END);

    //! call the callback function
    //! This function returns false if the run is to be aborted
    bool DoCallback(FEModel* fem, unsigned int nevent);

    //! Get the current callback reason (or zero if not inside DoCallback)
    unsigned int CurrentEvent() const;

private:
    std::list<FECORE_CALLBACK> m_pcb;  //!< pointer to callback function
    unsigned int m_event;              //!< reason for current callback (or zero)
};

