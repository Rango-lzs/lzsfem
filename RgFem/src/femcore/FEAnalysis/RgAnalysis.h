/*********************************************************************
 * \file   RgAnalysis.h
 * \brief  
 * 
 * \author Leizs
 * \date   March 2025
 *********************************************************************/

#pragma once

class AnalysisBase
{
public:
    virtual ~AnalysisBase() = default;

    // 核心接口
    virtual void initialize() = 0;      // 初始化分析
    virtual void execute() = 0;         // 执行计算
    virtual void finalize() = 0;        // 后处理
    virtual void validate() const = 0;  // 输入验证

    // 通用参数配置
    virtual void set_time_increment(double dt)
    {
        m_dt = dt;
    }
    virtual void set_output_interval(int interval)
    {
        m_output_interval = interval;
    }

protected:
    double m_dt = 0.0;          // 时间步长
    int m_output_interval = 1;  // 结果输出间隔
};

//---------------------- 派生类示例 ----------------------
class StaticAnalysis : public AnalysisBase
{
public:
    void initialize() override
    {
        // 初始化线性方程组求解器
        std::cout << "Initializing static solver...\n";
    }

    void execute() override
    {
        // 执行静力平衡迭代
        std::cout << "Solving static equilibrium...\n";
    }

    void validate() const override
    {
        // 验证静力分析参数
        if (m_dt != 0.0)
            throw std::runtime_error("Static analysis doesn't need time increment");
    }
    // ...其他接口实现
};

class ExplicitDynamicAnalysis : public AnalysisBase
{
public:
    void initialize() override
    {
        // 显式动力学特定初始化
        std::cout << "Initializing explicit time integration...\n";
    }

    void execute() override
    {
        // 显式时间步推进
        for (int step = 0; step < m_total_steps; ++step)
        {
            update_velocity();
            update_position();
            calculate_forces();
            // 条件输出
            if (step % m_output_interval == 0)
                write_results(step);
        }
    }

    void set_contact_algorithm(const std::string& algo)
    {
        m_contact_algo = algo;
    }

private:
    std::string m_contact_algo = "Penalty";
    int m_total_steps = 1000;

    void update_velocity()
    { /* 速度更新实现 */
    }
    void update_position()
    { /* 位置更新实现 */
    }
    void calculate_forces()
    { /* 内力计算实现 */
    }
    void write_results(int step)
    { /* 结果输出实现 */
    }
};