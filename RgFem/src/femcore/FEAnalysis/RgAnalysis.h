#pragma once
#include "femcore/FEObjectBase.h"

#include <memory>
#include <vector>

//-----------------------------------------------------------------------------
// 前置声明
class FEModel;
class FESolver;
class RgDomain;
class DumpStream;
class FEStepComponent;
class FETimeStepController;

//-----------------------------------------------------------------------------
//! 有限元分析步骤基类（职责分离、模块化设计）
class FEM_EXPORT FEAnalysis : public FEObjectBase
{
    //---------------------- 嵌套结构定义 ----------------------
    struct TimeStepData
    {
        int stepCount = 0;         // 时间步数
        double startTime = 0.0;    // 起始时间
        double endTime = 1.0;      // 终止时间
        double currentTime = 0.0;  // 当前时间
        double initialStep = 0.1;  // 初始步长
        double minStep = 1e-6;     // 最小步长
        double maxStep = 1.0;      // 最大步长
    };

    struct SolverMetrics
    {
        int totalRhsEvals = 0;       // 右端项计算次数
        int totalMatrixReforms = 0;  // 刚度矩阵重组次数
        int totalIterations = 0;     // 总迭代次数
    };

    struct OutputSettings
    {
        int plotLevel = 0;                          // 输出详细级别
        int outputStride = 1;                       // 输出间隔步长
        std::pair<int, int> plotRange{0, INT_MAX};  // 输出时间范围
        bool plotInitialState = true;               // 是否输出初始状态
        int plotHint = 0;                           // 输出模式标识
    };

    //---------------------- 核心接口 ----------------------
public:
    explicit FEAnalysis(FEModel* model);
    virtual ~FEAnalysis();

    // 生命周期管理
    bool initialize() override;
    void reset();
    void serialize(DumpStream& stream) override;

    // 执行控制
    enum class Status
    {
        Ready,
        Running,
        Converged,
        Failed
    };
    Status run();  // 主执行入口
    void abort();  // 终止分析

    // 域管理
    void addDomain(int domainId);
    void clearDomains();
    RgDomain* getDomain(size_t index) const;
    size_t domainCount() const noexcept;

    // 组件管理
    void addComponent(FEStepComponent* component);
    FEStepComponent* getComponent(size_t index) const;
    size_t componentCount() const noexcept;

    // 求解器配置
    void setSolver(FESolver* solver);
    FESolver* getSolver() const noexcept;

    // 时间步控制
    void configureTimeStepping(double start, double end, double initStep);
    void setTimeController(FETimeStepController* controller);

    // 输出配置
    void setOutputSettings(const OutputSettings& settings);
    const OutputSettings& getOutputSettings() const noexcept;

    //---------------------- 实现细节 ----------------------
private:
    bool validateModel() const;     // 模型验证
    void prepareStep();             // 步进准备
    void finalizeStep();            // 步进收尾
    void writeOutput(double time);  // 结果输出

    // 模块化子组件
    FEModel* const m_model;              // 所属模型（不可变）
    std::unique_ptr<FESolver> m_solver;  // 求解器（独占所有权）
    std::unique_ptr<FETimeStepController> m_timeController;

    // 数据分组
    TimeStepData m_timeData;
    SolverMetrics m_metrics;
    OutputSettings m_output;

    // 动态内容
    std::vector<int> m_activeDomains;                            // 激活的域ID列表
    std::vector<std::unique_ptr<FEStepComponent>> m_components;  // 步骤组件
    Status m_status = Status::Ready;                             // 当前状态
};
