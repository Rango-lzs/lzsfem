/*********************************************************************
 * \file   RgLoadController.h
 * \brief  Load controller for time-varying boundary conditions and loads
 *
 * \author Leizs
 * \date   January 2025
 *********************************************************************/

#pragma once
#include "femcore/FEObjectBase.h"
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
//! Load controller base class
//! Provides time-dependent scaling factor for loads and BCs
class FEM_EXPORT RgLoadController : public FEObjectBase
{
    DECLARE_META_CLASS(RgLoadController, FEObjectBase);

public:
    //! Constructor
    RgLoadController();

    //! Destructor
    virtual ~RgLoadController();

    //! Initialize
    virtual bool Init() override;

    //! Evaluate controller at given time
    virtual void Evaluate(double time);

    //! Get current value
    double Value() const { return m_value; }

    //! Serialize
    virtual void Serialize(DumpStream& ar) override;

    //! Set controller type
    void SetType(int type) { m_type = type; }

    //! Get controller type
    int GetType() const { return m_type; }

    //! Set interpolation type
    void SetInterpolation(int interp) { m_interpolation = interp; }

    //! Get interpolation type
    int GetInterpolation() const { return m_interpolation; }

    //! Set extend mode (what happens outside time range)
    void SetExtendMode(int mode) { m_extendMode = mode; }

    //! Get extend mode
    int GetExtendMode() const { return m_extendMode; }

public:
    //! Controller types
    enum Type {
        LC_LOADCURVE = 0,    //!< Piecewise curve
        LC_LINEAR,           //!< Linear ramp
        LC_STEP,             //!< Step function
        LC_SMOOTH_STEP,      //!< Smooth step
        LC_SINE,             //!< Sinusoidal
        LC_CONSTANT          //!< Constant value
    };

    //! Interpolation types
    enum Interpolation {
        INTERP_LINEAR = 0,   //!< Linear interpolation
        INTERP_STEP,         //!< Step (piecewise constant)
        INTERP_SMOOTH,       //!< Smooth (cubic spline)
        INTERP_AKIMA         //!< Akima spline
    };

    //! Extend modes
    enum ExtendMode {
        EXTEND_CONSTANT = 0, //!< Use constant end values
        EXTEND_EXTRAPOLATE,  //!< Extrapolate linearly
        EXTEND_REPEAT,       //!< Repeat the curve
        EXTEND_CYCLE         //!< Cycle the curve
    };

protected:
    double m_value;         //!< Current value
    int m_type;             //!< Controller type
    int m_interpolation;    //!< Interpolation type
    int m_extendMode;       //!< Extend mode
};

//-----------------------------------------------------------------------------
//! Piecewise load curve
class FEM_EXPORT RgLoadCurve : public RgLoadController
{
    DECLARE_META_CLASS(RgLoadCurve, RgLoadController);

public:
    //! Constructor
    RgLoadCurve();

    //! Initialize
    bool Init() override;

    //! Evaluate at given time
    void Evaluate(double time) override;

    //! Add a point to the curve
    void AddPoint(double time, double value);

    //! Clear all points
    void Clear();

    //! Get number of points
    int Points() const { return (int)m_points.size(); }

    //! Get point
    void GetPoint(int i, double& time, double& value) const;

    //! Set point
    void SetPoint(int i, double time, double value);

    //! Serialize
    void Serialize(DumpStream& ar) override;

private:
    //! Interpolate value at time t
    double Interpolate(double t);

    //! Linear interpolation
    double LinearInterpolate(double t);

    //! Step interpolation
    double StepInterpolate(double t);

    //! Smooth interpolation (cubic)
    double SmoothInterpolate(double t);

    //! Handle extend modes
    double ExtendValue(double t);

private:
    struct Point
    {
        double time;
        double value;
        
        Point() : time(0.0), value(0.0) {}
        Point(double t, double v) : time(t), value(v) {}
    };

    std::vector<Point> m_points;  //!< Control points
};

//-----------------------------------------------------------------------------
//! Linear ramp controller
class FEM_EXPORT RgLinearController : public RgLoadController
{
    DECLARE_META_CLASS(RgLinearController, RgLoadController);

public:
    RgLinearController();

    //! Set parameters
    void SetParameters(double t0, double t1, double v0, double v1);

    //! Evaluate
    void Evaluate(double time) override;

    //! Serialize
    void Serialize(DumpStream& ar) override;

private:
    double m_t0, m_t1;  //!< Time range
    double m_v0, m_v1;  //!< Value range
};

//-----------------------------------------------------------------------------
//! Step controller
class FEM_EXPORT RgStepController : public RgLoadController
{
    DECLARE_META_CLASS(RgStepController, RgLoadController);

public:
    RgStepController();

    //! Set parameters
    void SetParameters(double t_step, double v0, double v1);

    //! Evaluate
    void Evaluate(double time) override;

    //! Serialize
    void Serialize(DumpStream& ar) override;

private:
    double m_t_step;    //!< Step time
    double m_v0, m_v1;  //!< Values before/after step
};

//-----------------------------------------------------------------------------
//! Smooth step controller
class FEM_EXPORT RgSmoothStepController : public RgLoadController
{
    DECLARE_META_CLASS(RgSmoothStepController, RgLoadController);

public:
    RgSmoothStepController();

    //! Set parameters
    void SetParameters(double t0, double t1, double v0, double v1);

    //! Evaluate
    void Evaluate(double time) override;

    //! Serialize
    void Serialize(DumpStream& ar) override;

private:
    double m_t0, m_t1;  //!< Transition time range
    double m_v0, m_v1;  //!< Value range
};

//-----------------------------------------------------------------------------
//! Sinusoidal controller
class FEM_EXPORT RgSineController : public RgLoadController
{
    DECLARE_META_CLASS(RgSineController, RgLoadController);

public:
    RgSineController();

    //! Set parameters
    void SetParameters(double amplitude, double frequency, double phase, double offset);

    //! Evaluate
    void Evaluate(double time) override;

    //! Serialize
    void Serialize(DumpStream& ar) override;

private:
    double m_amplitude; //!< Amplitude
    double m_frequency; //!< Frequency (Hz)
    double m_phase;     //!< Phase shift (radians)
    double m_offset;    //!< DC offset
};

//-----------------------------------------------------------------------------
//! Constant controller
class FEM_EXPORT RgConstantController : public RgLoadController
{
    DECLARE_META_CLASS(RgConstantController, RgLoadController);

public:
    RgConstantController();

    //! Set constant value
    void SetValue(double value);

    //! Evaluate
    void Evaluate(double time) override;

    //! Serialize
    void Serialize(DumpStream& ar) override;

private:
    double m_constantValue;
};
