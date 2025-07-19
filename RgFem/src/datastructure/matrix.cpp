#include "Matrix.h"

#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <vector>
#include <iostream>
#include "Matrix3d.h"

// PIMPL实现细节
class Matrix::InnerData
{
public:
    Eigen::MatrixXd mat;

    InnerData() = default;
    InnerData(int nr, int nc)
        : mat(nr, nc)
    {
    }
    InnerData(const Eigen::MatrixXd& m)
        : mat(m)
    {
    }
};

// ================= 构造函数与析构函数 =================
Matrix::Matrix()
    : m_data(std::make_unique<InnerData>())
{
}
Matrix::Matrix(int nr, int nc)
    : m_data(std::make_unique<InnerData>(nr, nc))
{
}
Matrix::~Matrix() = default;  // unique_ptr自动处理析构

// 拷贝构造
Matrix::Matrix(const Matrix& m)
    : m_data(std::make_unique<InnerData>(m.m_data->mat))
{
}

// 移动构造
Matrix::Matrix(Matrix&& m) noexcept
    : m_data(std::move(m.m_data))
{
}

// 拷贝赋值
Matrix& Matrix::operator=(const Matrix& m)
{
    if (this != &m)
    {
        m_data->mat = m.m_data->mat;
    }
    return *this;
}

// 移动赋值
Matrix& Matrix::operator=(Matrix&& m) noexcept
{
    if (this != &m)
    {
        m_data = std::move(m.m_data);
    }
    return *this;
}

Matrix& Matrix::operator=(const Matrix3d& m)
{
    //m_data = std::move(m.m_matrix);
    return *this;
}

// ================= 矩阵基本操作 =================
void Matrix::resize(int nr, int nc)
{
    m_data->mat.resize(nr, nc);
    m_data->mat.setZero();
}

void Matrix::zero()
{
    m_data->mat.setZero();
}

int Matrix::rows() const
{
    return m_data->mat.rows();
}

int Matrix::columns() const
{
    return m_data->mat.cols();
}

// ================= 元素访问 =================
double* Matrix::operator[](int i)
{
    return m_data->mat.row(i).data();
}

const double* Matrix::operator[](int i) const
{
    return m_data->mat.row(i).data();
}

double& Matrix::operator()(int i, int j)
{
    return m_data->mat(i, j);
}

double Matrix::operator()(int i, int j) const
{
    return m_data->mat(i, j);
}

Matrix::operator double**()
{
    // 注意：仅当需要兼容旧代码时使用
    return const_cast<double**>(reinterpret_cast<const double**>(m_data->mat.data()));
}

// ================= 线性代数运算 =================
Matrix Matrix::transpose() const
{
    Matrix result;
    result.m_data->mat = m_data->mat.transpose();
    return result;
}

Matrix Matrix::inverse() const
{
    Matrix result;
    result.m_data->mat = m_data->mat.partialPivLu().inverse();
    return result;
}

Matrix Matrix::svd_inverse() const
{
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(m_data->mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Matrix result;
    result.m_data->mat = svd.solve(Eigen::MatrixXd::Identity(rows(), rows()));
    return result;
}

// ================= 矩阵运算 =================
Matrix Matrix::operator*(double s) const
{
    Matrix result;
    result.m_data->mat = m_data->mat * s;
    return result;
}

Matrix Matrix::operator*(const Matrix& m) const
{
    Matrix result;
    result.m_data->mat = m_data->mat * m.m_data->mat;
    return result;
}

Matrix Matrix::operator+(const Matrix& m) const
{
    Matrix result;
    result.m_data->mat = m_data->mat + m.m_data->mat;
    return result;
}

Matrix Matrix::operator-(const Matrix& m) const
{
    Matrix result;
    result.m_data->mat = m_data->mat - m.m_data->mat;
    return result;
}

Matrix& Matrix::operator*=(double s)
{
    m_data->mat *= s;
    return *this;
}

Matrix& Matrix::operator+=(const Matrix& m)
{
    m_data->mat += m.m_data->mat;
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& m)
{
    m_data->mat -= m.m_data->mat;
    return *this;
}

// ================= 求解器 =================
void Matrix::solve(std::vector<double>& x, const std::vector<double>& b)
{
    Eigen::Map<const Eigen::VectorXd> bMap(b.data(), b.size());
    Eigen::VectorXd xVec = m_data->mat.colPivHouseholderQr().solve(bMap);
    x.assign(xVec.data(), xVec.data() + xVec.size());
}

bool Matrix::lsq_solve(std::vector<double>& x, std::vector<double>& b)
{
    try
    {
        // 检查维度一致性
        if (m_data->mat.rows() != static_cast<int>(b.size()))
        {
            throw std::invalid_argument("Vector size does not match matrix rows");
        }

        // 创建 Eigen 向量映射（零拷贝）
        Eigen::Map<Eigen::VectorXd> bMap(b.data(), b.size());

        // 执行 SVD 分解并求解
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(m_data->mat, Eigen::ComputeThinU | Eigen::ComputeThinV);

        // 检查 SVD 是否成功（正确方式）
        if (svd.rank() == 0)
        {
            // 矩阵完全奇异的情况
            return false;
        }

        // 求解最小二乘问题
        Eigen::VectorXd xVec = svd.solve(bMap);

        // 检查求解状态（通过求解结果）
        if ((m_data->mat * xVec - bMap).norm() > 1e-8 * bMap.norm())
        {
            // 数值不稳定的情况
            return false;
        }

        // 复制结果到输出向量
        x.assign(xVec.data(), xVec.data() + xVec.size());
        return true;
    }
    catch (const std::exception& e)
    {
        // 处理任何异常（如内存不足）
        std::cerr << "lsq_solve error: " << e.what() << std::endl;
        return false;
    }
}

// ================= 子矩阵操作 =================
void Matrix::get(int i, int j, int rows, int cols, Matrix& A) const
{
    A.m_data->mat = m_data->mat.block(i, j, rows, cols);
}

void Matrix::fill(int i, int j, int rows, int cols, double val)
{
    m_data->mat.block(i, j, rows, cols).setConstant(val);
}

// ================= 特征值与乘法 =================
bool Matrix::eigen_vectors(Matrix& eigenVecs, std::vector<double>& eigenVals)
{
    if (rows() != columns())
        return false;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(m_data->mat);
    if (solver.info() != Eigen::Success)
        return false;

    eigenVecs.m_data->mat = solver.eigenvectors();
    const auto& evals = solver.eigenvalues();
    eigenVals.assign(evals.data(), evals.data() + evals.size());
    return true;
}

void Matrix::mult(std::vector<double>& x, std::vector<double>& y)
{
    if (y.size() < static_cast<size_t>(rows()))
    {
        y.resize(rows());
    }

    Eigen::Map<const Eigen::VectorXd> xMap(x.data(), x.size());
    Eigen::Map<Eigen::VectorXd> yMap(y.data(), y.size());
    yMap = m_data->mat * xMap;
}

// ================ 矩阵-向量乘法（成员函数） ================
std::vector<double> Matrix::operator*(const std::vector<double>& x) const
{
    // 检查输入有效性
    if (x.size() != static_cast<size_t>(columns()))
    {
        throw std::invalid_argument("Vector size does not match matrix columns");
    }

    // 创建结果向量（正确大小）
    std::vector<double> result(rows());

    // 使用Eigen进行高效计算
    Eigen::Map<const Eigen::VectorXd> xMap(x.data(), x.size());
    Eigen::Map<Eigen::VectorXd> resultMap(result.data(), result.size());
    resultMap = m_data->mat * xMap;

    return result;
}

// ================ 线性系统求解（成员函数） ================
std::vector<double> Matrix::solve(const std::vector<double>& b) const
{
    // 检查输入有效性
    if (b.size() != static_cast<size_t>(rows()))
    {
        throw std::invalid_argument("Vector size does not match matrix rows");
    }

    // 使用Eigen进行求解
    Eigen::Map<const Eigen::VectorXd> bMap(b.data(), b.size());
    Eigen::VectorXd xVec = m_data->mat.colPivHouseholderQr().solve(bMap);

    // 返回结果
    return std::vector<double>(xVec.data(), xVec.data() + xVec.size());
}

// ================= 全局函数实现 =================
std::vector<double> operator/(std::vector<double>& b, Matrix& m)
{
    return m.solve(b);
}


Matrix Matrix::OuterProduct(const std::vector<double>& a)
{
    Eigen::Map<const Eigen::VectorXd> aMap(a.data(), a.size());
    Matrix result;
    result.resize(a.size(), a.size());  // 确保大小正确
    result.m_data->mat = aMap * aMap.transpose();
    return result;
}

Matrix outer_product(std::vector<double>& a)
{
    return Matrix::OuterProduct(a);
}