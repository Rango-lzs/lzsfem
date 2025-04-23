#pragma once
#include "datastructure/Vector2d.h"
#include "femcore/fem_export.h"
#include <vector>

class FEM_EXPORT PointCurve
{
	class Imp;

public:
	//! Interpolation functions
	enum INTFUNC { LINEAR = 0, STEP = 1, SMOOTH = 2, CSPLINE = 3, CPOINTS = 4, APPROX = 5, SMOOTH_STEP = 6 };

	//! Extend mode
	enum EXTMODE { CONSTANT, EXTRAPOLATE, REPEAT, REPEAT_OFFSET };

public:
	//! default constructor
	PointCurve();

	//! copy constructor
	PointCurve(const PointCurve& pc);

	//! assignment operator
	void operator = (const PointCurve& pc);

	//! destructor
	~PointCurve();

	//! call this to update internal data structures
	bool Update();

	//! adds a point to the point curve
	int Add(double x, double y);

	//! adds a point to the point curve
	int Add(const Vector2d& p);

	//! Clears the loadcurve data
	void Clear();

	//! set the x and y value of point i
	void SetPoint(int i, double x, double y);
	void SetPoint(int i, const Vector2d& p);

	//! set all points at once
	void SetPoints(const std::vector<Vector2d>& points);

	//! return all points
	std::vector<Vector2d> GetPoints() const;

	//! remove a point
	void Delete(int n);

	//! remove several points at once
	void Delete(const std::vector<int>& indexList);

	//! Set the type of interpolation
	void SetInterpolator(int fnc);

	//! return current interpolator
	int GetInterpolator() const;

	//! Set the extend mode
	void SetExtendMode(int mode);

	//! Get the extend mode
	int GetExtendMode() const;

	//! get a point
	Vector2d Point(int i) const;

	//! finds closest load point
	int FindPoint(double t, double& tval, int startIndex = 0);

	//! return nr of points
	int Points() const;

	//! see if there is a point at time t
	bool HasPoint(double t) const;

public: // operations
	
	// scale all y points by s
	void Scale(double s);

public:

	//! returns the value of the load curve at time
	double value(double x) const;

	//! returns the derivative value at time
	double derive(double x) const;

	//! returns the second derivative value at time
	double deriv2(double x) const;

	//! returns the definite integral value between a and b
	double integrate(double a, double b) const;

protected:
	double ExtendValue(double t) const;

private:
	Imp* im;
};
