#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <stdexcept>

// ================ 构造函数 ================
Matrix3d::Matrix3d()
    : m_matrix(Eigen::Matrix3d::Zero())
{
}

Matrix3d::Matrix3d(const EigenMat3d& matrix)
    : m_matrix(matrix)
{
}

Matrix3d::Matrix3d(double a)
    : m_matrix(Eigen::Matrix3d::Identity() * a)
{
}

Matrix3d::Matrix3d(double a00, double a01, double a02, double a10, double a11, double a12, double a20, double a21,
                   double a22)
{
    m_matrix << a00, a01, a02, a10, a11, a12, a20, a21, a22;
}

Matrix3d::Matrix3d(double m[3][3])
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            m_matrix(i, j) = m[i][j];
}

Matrix3d::Matrix3d(double a[9])
{
    m_matrix << a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8];
}

Matrix3d::Matrix3d(const Matrix3dd& m)
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            m_matrix(i, j) = m(i,j);
}

Matrix3d::Matrix3d(const Matrix3ds& m)
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            m_matrix(i, j) = m(i, j);
}

Matrix3d::Matrix3d(const Matrix3da& m)
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            m_matrix(i, j) = m(i, j);
}

Matrix3d::Matrix3d(const Matrix2d& m)
{
    m_matrix.setZero();
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            m_matrix(i, j) = m(i, j);
        }
    }
    m_matrix(2, 2) = 1.0;
}

Matrix3d::Matrix3d(const Vector3d& e1, const Vector3d& e2, const Vector3d& e3)
{
    m_matrix.col(0) = Eigen::Vector3d(e1.x, e1.y, e1.z);
    m_matrix.col(1) = Eigen::Vector3d(e2.x, e2.y, e2.z);
    m_matrix.col(2) = Eigen::Vector3d(e3.x, e3.y, e3.z);
}

// ================ 赋值运算符 ================
Matrix3d& Matrix3d::operator=(const Matrix3d& m)
{
    if (this != &m)
    {
        m_matrix = m.m_matrix;
    }
    return *this;
}

Matrix3d& Matrix3d::operator=(const double m[3][3])
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            m_matrix(i, j) = m[i][j];
    return *this;
}

// ================ 一元运算符 ================
Matrix3d Matrix3d::operator-() const
{
    return Matrix3d(-m_matrix);
}

// ================ 访问运算符 ================
double& Matrix3d::operator()(int i, int j)
{
    if (i < 0 || i > 2 || j < 0 || j > 2)
    {
        throw std::out_of_range("Matrix3d index out of range");
    }
    return m_matrix(i, j);
}

const double& Matrix3d::operator()(int i, int j) const
{
    if (i < 0 || i > 2 || j < 0 || j > 2)
    {
        throw std::out_of_range("Matrix3d index out of range");
    }
    return m_matrix(i, j);
}

double* Matrix3d::operator[](int i)
{
    if (i < 0 || i > 2)
    {
        throw std::out_of_range("Matrix3d row index out of range");
    }
    return m_matrix.row(i).data();
}

const double* Matrix3d::operator[](int i) const
{
    if (i < 0 || i > 2)
    {
        throw std::out_of_range("Matrix3d row index out of range");
    }
    return m_matrix.row(i).data();
}

// ================ 算术运算符 ================
Matrix3d Matrix3d::operator+(const Matrix3d& m) const
{
    return Matrix3d(m_matrix + m.m_matrix);
}

Matrix3d Matrix3d::operator-(const Matrix3d& m) const
{
    return Matrix3d(m_matrix - m.m_matrix);
}

Matrix3d Matrix3d::operator*(const Matrix3d& m) const
{
    return Matrix3d(m_matrix * m.m_matrix);
}

Matrix3d Matrix3d::operator*(double a) const
{
    return Matrix3d(m_matrix * a);
}

Matrix3d Matrix3d::operator/(double a) const
{
    if (a == 0)
    {
        throw std::invalid_argument("Matrix3d division by zero");
    }
    return Matrix3d(m_matrix / a);
}

// ================ 算术赋值运算符 ================
Matrix3d& Matrix3d::operator+=(const Matrix3d& m)
{
    m_matrix += m.m_matrix;
    return *this;
}

Matrix3d& Matrix3d::operator-=(const Matrix3d& m)
{
    m_matrix -= m.m_matrix;
    return *this;
}

Matrix3d& Matrix3d::operator*=(const Matrix3d& m)
{
    m_matrix *= m.m_matrix;
    return *this;
}

Matrix3d& Matrix3d::operator*=(double a)
{
    m_matrix *= a;
    return *this;
}

Matrix3d& Matrix3d::operator/=(double a)
{
    if (a == 0)
    {
        throw std::invalid_argument("Matrix3d division by zero");
    }
    m_matrix /= a;
    return *this;
}

// ================ 矩阵-向量乘法 ================
Vector3d Matrix3d::operator*(const Vector3d& r) const
{
    Eigen::Vector3d result = m_matrix * Eigen::Vector3d(r.x, r.y, r.z);
    return Vector3d(result.x(), result.y(), result.z());
}

// ================ 矩阵属性 ================
double Matrix3d::det() const
{
    return m_matrix.determinant();
}

double Matrix3d::trace() const
{
    return m_matrix.trace();
}

// ================ 矩阵操作 ================
void Matrix3d::zero()
{
    m_matrix.setZero();
}

void Matrix3d::unit()
{
    m_matrix.setIdentity();
}

Vector3d Matrix3d::col(int j) const
{
    if (j < 0 || j > 2)
    {
        throw std::out_of_range("Matrix3d column index out of range");
    }
    Eigen::Vector3d c = m_matrix.col(j);
    return Vector3d(c.x(), c.y(), c.z());
}

Vector3d Matrix3d::row(int i) const
{
    if (i < 0 || i > 2)
    {
        throw std::out_of_range("Matrix3d row index out of range");
    }
    Eigen::Vector3d r = m_matrix.row(i);
    return Vector3d(r.x(), r.y(), r.z());
}

void Matrix3d::setCol(int j, const Vector3d& a)
{
    if (j < 0 || j > 2)
    {
        throw std::out_of_range("Matrix3d column index out of range");
    }
    m_matrix.col(j) = Eigen::Vector3d(a.x, a.y, a.z);
}

void Matrix3d::setRow(int i, const Vector3d& a)
{
    if (i < 0 || i > 2)
    {
        throw std::out_of_range("Matrix3d row index out of range");
    }
    m_matrix.row(i) = Eigen::Vector3d(a.x, a.y, a.z);
}

Matrix3ds Matrix3d::sym() const
{
    Eigen::Matrix3d sym = 0.5 * (m_matrix + m_matrix.transpose());
    return Matrix3ds(sym(0, 0), sym(1,1), sym(2, 2), sym(0, 1), sym(1, 2), sym(0, 2));
}

Matrix3da Matrix3d::skew() const
{
    Eigen::Matrix3d skew = 0.5 * (m_matrix - m_matrix.transpose());
    return Matrix3da(skew(0, 1), skew(1, 2), skew(0, 2));
}

Matrix3d Matrix3d::inverse() const
{
    if (det() < 1e-12)
    {
        throw std::runtime_error("Matrix3d is singular, cannot compute inverse");
    }
    return Matrix3d(m_matrix.inverse());
}

bool Matrix3d::invert()
{
    if (det() < 1e-12)
    {
        return false;
    }
    m_matrix = m_matrix.inverse().eval();
    return true;
}

Matrix3d Matrix3d::transpose() const
{
    return Matrix3d(m_matrix.transpose());
}

Matrix3d Matrix3d::transinv() const
{
    return transpose().inverse();
}

void Matrix3d::skew(const Vector3d& v)
{
    m_matrix << 0, -v.z, v.y, v.z, 0, -v.x, -v.y, v.x, 0;
}

double Matrix3d::norm() const
{
    return m_matrix.lpNorm<1>();
}

double Matrix3d::dotdot(const Matrix3d& T) const
{
    return (m_matrix.array() * T.m_matrix.array()).sum();
}

// ================ 极分解 ================
void Matrix3d::right_polar(Matrix3d& R, Matrix3ds& U) const
{
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(m_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);

    // 确保旋转矩阵是正交的
    R.m_matrix = svd.matrixU() * svd.matrixV().transpose();
    if (R.det() < 0)
    {
        // 处理反射情况
        Eigen::Matrix3d V = svd.matrixV();
        V.col(2) = -V.col(2);
        R.m_matrix = svd.matrixU() * V.transpose();
    }

    // 计算右拉伸张量
    auto mat = svd.matrixV() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();
    U = Matrix3d(mat).sym();
}

void Matrix3d::left_polar(Matrix3ds& V, Matrix3d& R) const
{
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(m_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);

    // 确保旋转矩阵是正交的
    R.m_matrix = svd.matrixU() * svd.matrixV().transpose();
    if (R.det() < 0)
    {
        // 处理反射情况
        Eigen::Matrix3d U = svd.matrixU();
        U.col(2) = -U.col(2);
        R.m_matrix = U * svd.matrixV().transpose();
    }

    // 计算左拉伸张量
    auto mat = svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixU().transpose();
    V = Matrix3d(mat).sym();
}

// ================ 静态方法 ================
Matrix3d Matrix3d::identity()
{
    return Matrix3d(Eigen::Matrix3d::Identity());
}

// ================ 全局运算符 ================
Matrix3d operator*(double s, const Matrix3d& m)
{
    return m * s;
}

Matrix3d operator&(const Vector3d& a, const Vector3d& b)
{
    return Matrix3d(a.x * b.x, a.x * b.y, a.x * b.z, a.y * b.x, a.y * b.y, a.y * b.z, a.z * b.x, a.z * b.y, a.z * b.z);
}

Matrix3d skew(const Vector3d& a)
{
    return Matrix3d(0, -a.z, a.y, a.z, 0, -a.x, -a.y, a.x, 0);
}
