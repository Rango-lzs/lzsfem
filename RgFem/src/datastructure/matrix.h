#pragma once
#include "femcore/fem_export.h"
#include <memory>
#include <vector>


class Matrix3d;
class Vector3d;

class FEM_EXPORT Matrix
{
public:
    // ����������
    Matrix();
    Matrix(int nr, int nc);
    Matrix(const Matrix& m);
    Matrix(Matrix&& m) noexcept;
    ~Matrix();

    // ��ֵ������
    Matrix& operator=(const Matrix& m);
    Matrix& operator=(Matrix&& m) noexcept;
    //! assignment operator
    Matrix& operator=(const Matrix3d& m);

    // �������
    void resize(int nr, int nc);
    void zero();
    int rows() const;
    int columns() const;

    // Ԫ�ط���
    double* operator[](int i);
    const double* operator[](int i) const;
    double& operator()(int i, int j);
    double operator()(int i, int j) const;

    //��Matrixת����double**
    operator double**();

    // ���Դ�������
    Matrix transpose() const;
    Matrix inverse() const;
    Matrix svd_inverse() const;

    // ��������
    Matrix operator*(double s) const;
    Matrix operator*(const Matrix& m) const;
    Matrix operator+(const Matrix& m) const;
    Matrix operator-(const Matrix& m) const;
    Matrix& operator*=(double s);
    Matrix& operator+=(const Matrix& m);
    Matrix& operator-=(const Matrix& m);

    // ��Ϊ��Ա�����ľ���-�����˷�
    std::vector<double> operator*(const std::vector<double>& x) const;

    // �������ϵͳ�ĳ�Ա�������� operator/ ʹ�ã�
    std::vector<double> solve(const std::vector<double>& b) const;

    // �����
    void solve(std::vector<double>& x, const std::vector<double>& b);
    bool lsq_solve(std::vector<double>& x, std::vector<double>& b);

    // �Ӿ������
    void get(int i, int j, int rows, int cols, Matrix& A) const;
    void fill(int i, int j, int rows, int cols, double val);

    // ����ֵ��˷�
    bool eigen_vectors(Matrix& eigenVecs, std::vector<double>& eigenVals);
    void mult(std::vector<double>& x, std::vector<double>& y);

    static Matrix OuterProduct(const std::vector<double>& a);

private:
    class InnerData;                    // ǰ������
    std::unique_ptr<InnerData> m_data;  // PIMPL����
};

// ȫ�ֺ�������
std::vector<double> FEM_EXPORT operator/(std::vector<double>& b, Matrix& m);
// ͨ���������matrix m(i,j) = a(i)*a(j)
Matrix FEM_EXPORT outer_product(std::vector<double>& a);
