#pragma once
#include <math.h>
#include <memory.h>
#include <vector>
#include <algorithm>
#include "Vector3d.h"
#include "femcore/fem_export.h"

class FEMesh;
class FEDofList;

double FEM_EXPORT operator*(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> FEM_EXPORT operator - (std::vector<double>& a, std::vector<double>& b);
template<typename T> void zero(std::vector<T>& a) { fill(a.begin(), a.end(), T(0)); }
template<> inline void zero<Vector3d>(std::vector<Vector3d>& a) { fill(a.begin(), a.end(), Vector3d(0,0,0)); }
template<typename T> void assign(std::vector<T>& a, const T& v) { fill(a.begin(), a.end(), v); }
void FEM_EXPORT operator+=(std::vector<double>& a, const std::vector<double>& b);
void FEM_EXPORT operator-=(std::vector<double>& a, const std::vector<double>& b);
void FEM_EXPORT operator*=(std::vector<double>& a, double b);
std::vector<double> FEM_EXPORT operator+(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> FEM_EXPORT operator*(const std::vector<double>& a, double g);
std::vector<double> FEM_EXPORT operator - (const std::vector<double>& a);

// copy vector and scale
void FEM_EXPORT vcopys(std::vector<double>& a, const std::vector<double>& b, double s);

// add scaled vector
void FEM_EXPORT vadds(std::vector<double>& a, const std::vector<double>& b, double s);
void FEM_EXPORT vsubs(std::vector<double>& a, const std::vector<double>& b, double s);

// vector subtraction: a = l - r
void FEM_EXPORT vsub(std::vector<double>& a, const std::vector<double>& l, const std::vector<double>& r);

// scale each component of a vector
void FEM_EXPORT vscale(std::vector<double>& a, const std::vector<double>& s);

// gather operation (copy mesh data to vector)
void FEM_EXPORT gather(std::vector<double>& v, FEMesh& mesh, int ndof);
void FEM_EXPORT gather(std::vector<double>& v, FEMesh& mesh, const std::vector<int>& dof);

// scatter operation (copy vector data to mesh)
void FEM_EXPORT scatter(std::vector<double>& v, FEMesh& mesh, int ndof);
void FEM_EXPORT scatter3(std::vector<double>& v, FEMesh& mesh, int ndof1, int ndof2, int ndof3);
void FEM_EXPORT scatter(std::vector<double>& v, FEMesh& mesh, const FEDofList& dofs);

// calculate l2 norm of vector
double FEM_EXPORT l2_norm(const std::vector<double>& v);
double FEM_EXPORT l2_sqrnorm(const std::vector<double>& v);
double l2_norm(double* x, int n);
