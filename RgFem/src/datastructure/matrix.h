#pragma once
#include "femcore/fem_export.h"
#include <memory>
#include <vector>


class Matrix3d;
class Vector3d;

class FEM_EXPORT Matrix
{
public:
    // 构造与析构
    Matrix();
    Matrix(int nr, int nc);
    Matrix(const Matrix& m);
    Matrix(Matrix&& m) noexcept;
    ~Matrix();

    // 赋值操作符
    Matrix& operator=(const Matrix& m);
    Matrix& operator=(Matrix&& m) noexcept;
    //! assignment operator
    Matrix& operator=(const Matrix3d& m);

    // 矩阵操作
    void resize(int nr, int nc);
    void zero();
    int rows() const;
    int columns() const;

    // 元素访问
    double* operator[](int i);
    const double* operator[](int i) const;
    double& operator()(int i, int j);
    double operator()(int i, int j) const;

    //将Matrix转换成double**
    operator double**();

    // 线性代数运算
    Matrix transpose() const;
    Matrix inverse() const;
    Matrix svd_inverse() const;

    // 矩阵运算
    Matrix operator*(double s) const;
    Matrix operator*(const Matrix& m) const;
    Matrix operator+(const Matrix& m) const;
    Matrix operator-(const Matrix& m) const;
    Matrix& operator*=(double s);
    Matrix& operator+=(const Matrix& m);
    Matrix& operator-=(const Matrix& m);

    // 作为成员函数的矩阵-向量乘法
    std::vector<double> operator*(const std::vector<double>& x) const;

    // 求解线性系统的成员函数（供 operator/ 使用）
    std::vector<double> solve(const std::vector<double>& b) const;

    // 求解器
    void solve(std::vector<double>& x, const std::vector<double>& b);
    bool lsq_solve(std::vector<double>& x, std::vector<double>& b);

    // 子矩阵操作
    void get(int i, int j, int rows, int cols, Matrix& A) const;
    void fill(int i, int j, int rows, int cols, double val);

    // 特征值与乘法
    bool eigen_vectors(Matrix& eigenVecs, std::vector<double>& eigenVals);
    void mult(std::vector<double>& x, std::vector<double>& y);

    static Matrix OuterProduct(const std::vector<double>& a);

private:
    class InnerData;                    // 前置声明
    std::unique_ptr<InnerData> m_data;  // PIMPL核心
};

// 全局函数声明
std::vector<double> FEM_EXPORT operator/(std::vector<double>& b, Matrix& m);

Matrix FEM_EXPORT operator*(double a, Matrix& m);

// 通过外积构造matrix m(i,j) = a(i)*a(j)
Matrix FEM_EXPORT outer_product(std::vector<double>& a);
