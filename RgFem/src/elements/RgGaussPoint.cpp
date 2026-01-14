#include "RgGaussPoint.h"
#include <cassert>



// ============================================================================
// Constructors
// ============================================================================

RgGaussPoint::RgGaussPoint()
    : m_coordinates(3, 0.0)  // Default to 3D with zero coordinates
    , m_weight(0.0)
{
}

RgGaussPoint::RgGaussPoint(const std::vector<double>& coordinates, double weight)
    : m_coordinates(coordinates)
    , m_weight(weight)
{
}

RgGaussPoint::RgGaussPoint(double r, double s, double t, double weight)
    : m_coordinates({r, s, t})
    , m_weight(weight)
{
}

RgGaussPoint::RgGaussPoint(double r, double s, double weight)
    : m_coordinates({r, s})
    , m_weight(weight)
{
}

RgGaussPoint::RgGaussPoint(double r, double weight)
    : m_coordinates({r})
    , m_weight(weight)
{
}

// ============================================================================
// Dimension Info
// ============================================================================

int RgGaussPoint::getDimension() const
{
    return static_cast<int>(m_coordinates.size());
}

// ============================================================================
// Coordinate Accessors
// ============================================================================

const std::vector<double>& RgGaussPoint::getCoordinates() const
{
    return m_coordinates;
}

double RgGaussPoint::getR() const
{
    return m_coordinates.empty() ? 0.0 : m_coordinates[0];
}

double RgGaussPoint::getS() const
{
    return (m_coordinates.size() < 2) ? 0.0 : m_coordinates[1];
}

double RgGaussPoint::getT() const
{
    return (m_coordinates.size() < 3) ? 0.0 : m_coordinates[2];
}

double RgGaussPoint::getCoordinate(int index) const
{
    assert(index >= 0 && index < static_cast<int>(m_coordinates.size()));
    return (index < static_cast<int>(m_coordinates.size())) ? m_coordinates[index] : 0.0;
}

// ============================================================================
// Weight Accessor
// ============================================================================

double RgGaussPoint::getWeight() const
{
    return m_weight;
}

// ============================================================================
// Mutators
// ============================================================================

void RgGaussPoint::setCoordinates(const std::vector<double>& coordinates)
{
    m_coordinates = coordinates;
}

void RgGaussPoint::setCoordinates(double r, double s, double t)
{
    m_coordinates.resize(3);
    m_coordinates[0] = r;
    m_coordinates[1] = s;
    m_coordinates[2] = t;
}

void RgGaussPoint::setCoordinates(double r, double s)
{
    m_coordinates.resize(2);
    m_coordinates[0] = r;
    m_coordinates[1] = s;
}

void RgGaussPoint::setWeight(double weight)
{
    m_weight = weight;
}

