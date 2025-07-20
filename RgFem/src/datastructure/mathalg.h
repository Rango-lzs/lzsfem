#pragma once
#include "datastructure/Matrix3d.h"
#include <functional>
#include "femcore/fem_export.h"

template <class T> T weightedAverage(T* d, double* w, int n)
{
	T s = d[0] * w[0];
	for (int i = 1; i < n; ++i) s += d[i] * w[i];
	return s;
}

template <class T> T weightedAverage(T* d, double* w, int n, std::function<T(const T&)> fnc)
{
	T s = fnc(d[0]) * w[0];
	for (int i = 1; i < n; ++i) s += fnc(d[i]) * w[i];
	return s;
}

FEM_EXPORT Matrix3ds weightedAverageStructureTensor(Matrix3ds* d, double* w, int n);

// evaluate Log_p (X)
FEM_EXPORT Matrix3ds Log(const Matrix3ds& p, const Matrix3ds& X);
FEM_EXPORT Matrix3ds Exp(const Matrix3ds& p, const Matrix3ds& X);
