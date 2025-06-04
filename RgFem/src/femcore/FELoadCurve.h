#pragma once
#include "FELoadController.h"
#include "PointCurve.h"

//-----------------------------------------------------------------------------
// Base class for load curves.
// Load curves are used to manipulate the time dependency of model parameters.
class FEM_EXPORT FELoadCurve : public FELoadController
{
public:
    // constructor
    FELoadCurve();
    FELoadCurve(const FELoadCurve& lc);

    void operator=(const FELoadCurve& lc);

    // destructor
    virtual ~FELoadCurve();

    void Serialize(DumpStream& ar) override;

    bool CopyFrom(FELoadCurve* lc);

    void Add(double time, double value);

    void Clear();

    bool Init() override;

    PointCurve& GetFunction()
    {
        return m_fnc;
    }

    int GetInterpolation() const
    {
        return m_int;
    }
    void SetInterpolation(PointCurve::INTFUNC f);

    int GetExtendMode() const
    {
        return m_ext;
    }
    void SetExtendMode(PointCurve::EXTMODE f);

    std::vector<Vector2d> GetPoints() const
    {
        return m_points;
    }

    double GetValue(double time) override;

private:
    int m_int;
    int m_ext;
    std::vector<Vector2d> m_points;

private:
    PointCurve m_fnc;  //!< function to evaluate

    DECLARE_PARAM_LIST();
};
