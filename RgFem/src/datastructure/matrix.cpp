#include "Matrix.h"

#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <vector>
#include <iostream>
#include "Matrix3d.h"

// PIMPLʵ��ϸ��
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

// ================= ���캯������������ =================
Matrix::Matrix()
    : m_data(std::make_unique<InnerData>())
{
}
Matrix::Matrix(int nr, int nc)
    : m_data(std::make_unique<InnerData>(nr, nc))
{
}
Matrix::~Matrix() = default;  // unique_ptr�Զ���������

// ��������
Matrix::Matrix(const Matrix& m)
    : m_data(std::make_unique<InnerData>(m.m_data->mat))
{
}

// �ƶ�����
Matrix::Matrix(Matrix&& m) noexcept
    : m_data(std::move(m.m_data))
{
}

// ������ֵ
Matrix& Matrix::operator=(const Matrix& m)
{
    if (this != &m)
    {
        m_data->mat = m.m_data->mat;
    }
    return *this;
}

// �ƶ���ֵ
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

// ================= ����������� =================
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

// ================= Ԫ�ط��� =================
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
    // ע�⣺������Ҫ���ݾɴ���ʱʹ��
    return const_cast<double**>(reinterpret_cast<const double**>(m_data->mat.data()));
}

// ================= ���Դ������� =================
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

// ================= �������� =================
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

// ================= ����� =================
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
        // ���ά��һ����
        if (m_data->mat.rows() != static_cast<int>(b.size()))
        {
            throw std::invalid_argument("Vector size does not match matrix rows");
        }

        // ���� Eigen ����ӳ�䣨�㿽����
        Eigen::Map<Eigen::VectorXd> bMap(b.data(), b.size());

        // ִ�� SVD �ֽⲢ���
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(m_data->mat, Eigen::ComputeThinU | Eigen::ComputeThinV);

        // ��� SVD �Ƿ�ɹ�����ȷ��ʽ��
        if (svd.rank() == 0)
        {
            // ������ȫ��������
            return false;
        }

        // �����С��������
        Eigen::VectorXd xVec = svd.solve(bMap);

        // ������״̬��ͨ���������
        if ((m_data->mat * xVec - bMap).norm() > 1e-8 * bMap.norm())
        {
            // ��ֵ���ȶ������
            return false;
        }

        // ���ƽ�����������
        x.assign(xVec.data(), xVec.data() + xVec.size());
        return true;
    }
    catch (const std::exception& e)
    {
        // �����κ��쳣�����ڴ治�㣩
        std::cerr << "lsq_solve error: " << e.what() << std::endl;
        return false;
    }
}

// ================= �Ӿ������ =================
void Matrix::get(int i, int j, int rows, int cols, Matrix& A) const
{
    A.m_data->mat = m_data->mat.block(i, j, rows, cols);
}

void Matrix::fill(int i, int j, int rows, int cols, double val)
{
    m_data->mat.block(i, j, rows, cols).setConstant(val);
}

// ================= ����ֵ��˷� =================
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

// ================ ����-�����˷�����Ա������ ================
std::vector<double> Matrix::operator*(const std::vector<double>& x) const
{
    // ���������Ч��
    if (x.size() != static_cast<size_t>(columns()))
    {
        throw std::invalid_argument("Vector size does not match matrix columns");
    }

    // ���������������ȷ��С��
    std::vector<double> result(rows());

    // ʹ��Eigen���и�Ч����
    Eigen::Map<const Eigen::VectorXd> xMap(x.data(), x.size());
    Eigen::Map<Eigen::VectorXd> resultMap(result.data(), result.size());
    resultMap = m_data->mat * xMap;

    return result;
}

// ================ ����ϵͳ��⣨��Ա������ ================
std::vector<double> Matrix::solve(const std::vector<double>& b) const
{
    // ���������Ч��
    if (b.size() != static_cast<size_t>(rows()))
    {
        throw std::invalid_argument("Vector size does not match matrix rows");
    }

    // ʹ��Eigen�������
    Eigen::Map<const Eigen::VectorXd> bMap(b.data(), b.size());
    Eigen::VectorXd xVec = m_data->mat.colPivHouseholderQr().solve(bMap);

    // ���ؽ��
    return std::vector<double>(xVec.data(), xVec.data() + xVec.size());
}

// ================= ȫ�ֺ���ʵ�� =================
std::vector<double> operator/(std::vector<double>& b, Matrix& m)
{
    return m.solve(b);
}


Matrix Matrix::OuterProduct(const std::vector<double>& a)
{
    Eigen::Map<const Eigen::VectorXd> aMap(a.data(), a.size());
    Matrix result;
    result.resize(a.size(), a.size());  // ȷ����С��ȷ
    result.m_data->mat = aMap * aMap.transpose();
    return result;
}

Matrix outer_product(std::vector<double>& a)
{
    return Matrix::OuterProduct(a);
}