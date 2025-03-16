#pragma once
#include "femcore/FEObjectBase.h"

#include <memory>
#include <vector>

//-----------------------------------------------------------------------------
// ǰ������
class FEModel;
class FESolver;
class FEDomain;
class DumpStream;
class FEStepComponent;
class FETimeStepController;

//-----------------------------------------------------------------------------
//! ����Ԫ����������ְࣨ����롢ģ�黯��ƣ�
class FEM_EXPORT FEAnalysis : public FEObjectBase
{
    //---------------------- Ƕ�׽ṹ���� ----------------------
    struct TimeStepData
    {
        int stepCount = 0;         // ʱ�䲽��
        double startTime = 0.0;    // ��ʼʱ��
        double endTime = 1.0;      // ��ֹʱ��
        double currentTime = 0.0;  // ��ǰʱ��
        double initialStep = 0.1;  // ��ʼ����
        double minStep = 1e-6;     // ��С����
        double maxStep = 1.0;      // ��󲽳�
    };

    struct SolverMetrics
    {
        int totalRhsEvals = 0;       // �Ҷ���������
        int totalMatrixReforms = 0;  // �նȾ����������
        int totalIterations = 0;     // �ܵ�������
    };

    struct OutputSettings
    {
        int plotLevel = 0;                          // �����ϸ����
        int outputStride = 1;                       // ����������
        std::pair<int, int> plotRange{0, INT_MAX};  // ���ʱ�䷶Χ
        bool plotInitialState = true;               // �Ƿ������ʼ״̬
        int plotHint = 0;                           // ���ģʽ��ʶ
    };

    //---------------------- ���Ľӿ� ----------------------
public:
    explicit FEAnalysis(FEModel* model);
    virtual ~FEAnalysis();

    // �������ڹ���
    bool initialize() override;
    void reset();
    void serialize(DumpStream& stream) override;

    // ִ�п���
    enum class Status
    {
        Ready,
        Running,
        Converged,
        Failed
    };
    Status run();  // ��ִ�����
    void abort();  // ��ֹ����

    // �����
    void addDomain(int domainId);
    void clearDomains();
    FEDomain* getDomain(size_t index) const;
    size_t domainCount() const noexcept;

    // �������
    void addComponent(FEStepComponent* component);
    FEStepComponent* getComponent(size_t index) const;
    size_t componentCount() const noexcept;

    // ���������
    void setSolver(FESolver* solver);
    FESolver* getSolver() const noexcept;

    // ʱ�䲽����
    void configureTimeStepping(double start, double end, double initStep);
    void setTimeController(FETimeStepController* controller);

    // �������
    void setOutputSettings(const OutputSettings& settings);
    const OutputSettings& getOutputSettings() const noexcept;

    //---------------------- ʵ��ϸ�� ----------------------
private:
    bool validateModel() const;     // ģ����֤
    void prepareStep();             // ����׼��
    void finalizeStep();            // ������β
    void writeOutput(double time);  // ������

    // ģ�黯�����
    FEModel* const m_model;              // ����ģ�ͣ����ɱ䣩
    std::unique_ptr<FESolver> m_solver;  // ���������ռ����Ȩ��
    std::unique_ptr<FETimeStepController> m_timeController;

    // ���ݷ���
    TimeStepData m_timeData;
    SolverMetrics m_metrics;
    OutputSettings m_output;

    // ��̬����
    std::vector<int> m_activeDomains;                            // �������ID�б�
    std::vector<std::unique_ptr<FEStepComponent>> m_components;  // �������
    Status m_status = Status::Ready;                             // ��ǰ״̬
};
