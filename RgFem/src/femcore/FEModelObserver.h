#pragma once
#include "femcore/fem_export.h"

#include <list>
#include <memory>
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


// Model 可以通过notify发送一个待特殊Event的消息给所有观察者，观察者只处理自己关注的消息
//  使用强类型枚举替代原生的宏定义
enum class FEModelEvent
{
    Init = 0x00000001,
    StepActive = 0x00000002,
    MajorIteration = 0x00000004,
    MinorIteration = 0x00000008,
    // ... 其他事件类型
    UserDefined1 = 0x01000000
};

class FEObserver
{
public:
    virtual ~FEObserver() = default;

    // 返回true表示继续执行，false表示中断流程
    virtual bool handleEvent(FEModel* model, FEModelEvent event) = 0;
};


class FESubject
{
public:
    virtual ~FESubject() = default;

    // 注册观察者（带优先级控制）
    virtual void attach(FEObserver* observer, int priority = 0) = 0;

    // 注销观察者
    virtual void detach(FEObserver* observer) = 0;

    // 通知观察者
    virtual bool notify(FEModelEvent event) = 0;
};


class FEFilteredObserver : public FEObserver
{
public:
    explicit FEFilteredObserver(FEModelEvent relevantEvents)
        : m_relevantEvents(static_cast<uint32_t>(relevantEvents))
    {
    }

    bool handleEvent(FEModel* model, FEModelEvent event) override
    {
        if (static_cast<uint32_t>(event) & m_relevantEvents)
        {
            return onEvent(model, event);
        }
        return true;
    }

protected:
    virtual bool onEvent(FEModel* model, FEModelEvent event) = 0;

private:
    uint32_t m_relevantEvents;
};

class FEModelSubject : public FESubject
{
public:
    void attach(FEObserver* observer, int priority = 0) override
    {
        m_observers.emplace_back(observer);
    }

    void detach(FEObserver* observer) override
    {
        auto it = std::find_if(m_observers.begin(), m_observers.end(),
                               [observer](const auto& curOb) { return curOb == observer; });

        if (it != m_observers.end())
        {
            m_observers.erase(it);
        }
    }

    bool notify(FEModelEvent event) override
    {
        // 按优先级降序通知
        for (auto it = m_observers.begin(); it != m_observers.end(); ++it)
        {
            if (!(*it)->handleEvent(m_model, event))
            {
                return false;
            }
        }
        return true;
    }

private:
    FEModel* m_model;
    std::vector<FEObserver*> m_observers;
};

class ConvergenceMonitor : public FEFilteredObserver
{
public:
    ConvergenceMonitor()
        : FEFilteredObserver(FEModelEvent::MajorIteration)
    {
    }

protected:
    bool onEvent(FEModel* model, FEModelEvent event) override
    {
        switch (event)
        {
            case FEModelEvent::MajorIteration:
                // logConvergence(model->getResidual());
                break;
            case FEModelEvent::MinorIteration:
                // updateProgress(model->getIteration());
                break;
        }
        return true;
    }
};

class OutputStrategy
{
public:
    virtual ~OutputStrategy() = default;
    virtual void write(const FEModelEvent& event) = 0;
};

// VTK格式输出
class VTKOutput : public OutputStrategy
{
public:
    void write(const FEModelEvent& event) override
    {
        // const auto& results = event.getResults();
        // std::string filename = "result_" + std::to_string(event.getStep()) + ".vtk";

        //// 实现VTK文件写入逻辑
        // writeVTKFile(filename, results);
    }
};

// HDF5格式输出
class HDF5Output : public OutputStrategy
{
public:
    void write(const FEModelEvent& event) override
    {
        // HDF5写入实现
    }
};


class OutputObserver : public FEFilteredObserver
{
public:
protected:
    bool onEvent(FEModel* model, FEModelEvent event) override
    {
        switch (event)
        {
            case FEModelEvent::MajorIteration:
                // logConvergence(model->getResidual());
                {
                    if (shouldWrite(event))
                    {
                        outputStrategy_->write(event);
                    }
                }
                break;
            case FEModelEvent::MinorIteration:
                // updateProgress(model->getIteration());
                break;
        }
        return true;
    }

    void setOutputStrategy(std::unique_ptr<OutputStrategy> strategy)
    {
        outputStrategy_ = std::move(strategy);
    }

private:
    bool shouldWrite(const FEModelEvent& event)
    {
        // 根据步数、时间等条件判断是否需要输出
        // return (event.getStep() % outputInterval_ == 0);
        return true;
    }

    std::unique_ptr<OutputStrategy> outputStrategy_;
    int outputInterval_ = 10;  // 每10步输出一次
};
