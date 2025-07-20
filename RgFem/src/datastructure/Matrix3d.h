#ifndef MATRIX3D_H
#define MATRIX3D_H

#include <Eigen/Dense>

//class Matrix3d;    general 3D Matrixrix of doubles 3*3
class Matrix3ds;  // symmetric 3D Matrixrix of doubles
class Matrix3da;  // anti-symmetric 3D Matrixrix of doubles
class Matrix3dd;  // diagonal Matrixrix of doubles
class Matrix2d;   // general 2D Matrixrix of doubles 2*2
class Vector3d;

using EigenMat3d = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;

class Matrix3d
{
public:
    // 构造函数
    Matrix3d();
    Matrix3d(const EigenMat3d& matrix);
    explicit Matrix3d(double a);
    Matrix3d(double a00, double a01, double a02, double a10, double a11, double a12, double a20, double a21,
             double a22);
    Matrix3d(double m[3][3]);
    Matrix3d(double a[9]);

    Matrix3d(const Matrix3dd& m);
    Matrix3d(const Matrix3ds& m);
    Matrix3d(const Matrix3da& m);

    Matrix3d(const Matrix2d& m);
    Matrix3d(const Vector3d& e1, const Vector3d& e2, const Vector3d& e3);

    // 赋值运算符
    Matrix3d& operator=(const Matrix3d& m);
    Matrix3d& operator=(const double m[3][3]);

    // 一元运算符
    Matrix3d operator-() const;

    // 访问运算符
    double& operator()(int i, int j);
    const double& operator()(int i, int j) const;
    double* operator[](int i);
    const double* operator[](int i) const;

    // 算术运算符
    Matrix3d operator+(const Matrix3d& m) const;
    Matrix3d operator-(const Matrix3d& m) const;
    Matrix3d operator*(const Matrix3d& m) const;
    Matrix3d operator*(double a) const;
    Matrix3d operator/(double a) const;

    // 算术赋值运算符
    Matrix3d& operator+=(const Matrix3d& m);
    Matrix3d& operator-=(const Matrix3d& m);
    Matrix3d& operator*=(const Matrix3d& m);
    Matrix3d& operator*=(double a);
    Matrix3d& operator/=(double a);

    // 矩阵-向量乘法
    Vector3d operator*(const Vector3d& r) const;

    // 矩阵属性
    double det() const;
    double trace() const;

    // 矩阵操作
    void zero();
    void unit();
    Vector3d col(int j) const;
    Vector3d row(int j) const;
    void setCol(int j, const Vector3d& a);
    void setRow(int i, const Vector3d& a);
    Matrix3ds sym() const;
    Matrix3da skew() const;
    Matrix3d inverse() const;
    bool invert();
    Matrix3d transpose() const;
    Matrix3d transinv() const;
    void skew(const Vector3d& v);
    double norm() const;
    double dotdot(const Matrix3d& T) const;

    // 极分解
    void right_polar(Matrix3d& R, Matrix3ds& U) const;
    void left_polar(Matrix3ds& V, Matrix3d& R) const;

    // 静态方法
    static Matrix3d identity();

private:
    EigenMat3d m_matrix;
};

// 全局运算符
Matrix3d operator*(double s, const Matrix3d& m);

// outer product for vectors
Matrix3d operator&(const Vector3d& a, const Vector3d& b);

// skew-symmetric Matrixrix of dual vector
Matrix3d skew(const Vector3d& a);

#endif  // MATRIX3D_H

#include "datastructure/Matrix3dSpecial.h"
