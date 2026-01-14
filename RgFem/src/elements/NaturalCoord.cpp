#include "NaturalCoord.h"
#include "RgGaussPoint.h"
#include <cassert>



// ============================================================================
// Constructors
// ============================================================================

NaturalCoord::NaturalCoord()
    : m_coordinates(3, 0.0)  // Default to 3D with zero coordinates
{
}

NaturalCoord::NaturalCoord(const std::vector<double>& coordinates)
    : m_coordinates(coordinates)
{
}

NaturalCoord::NaturalCoord(double r, double s, double t)
    : m_coordinates({r, s, t})
{
}

NaturalCoord::NaturalCoord(double r, double s)
    : m_coordinates({r, s})
{
}

NaturalCoord::NaturalCoord(double r)
    : m_coordinates({r})
{
}

NaturalCoord::NaturalCoord(const RgGaussPoint& gaussPoint)
    : m_coordinates(gaussPoint.getCoordinates())
{
}

// ============================================================================
// Dimension Info
// ============================================================================

int NaturalCoord::getDimension() const
{
    return static_cast<int>(m_coordinates.size());
}

// ============================================================================
// Coordinate Accessors
// ============================================================================

const std::vector<double>& NaturalCoord::getCoordinates() const
{
    return m_coordinates;
}

double NaturalCoord::getR() const
{
    return m_coordinates.empty() ? 0.0 : m_coordinates[0];
}

double NaturalCoord::getS() const
{
    return (m_coordinates.size() < 2) ? 0.0 : m_coordinates[1];
}

double NaturalCoord::getT() const
{
    return (m_coordinates.size() < 3) ? 0.0 : m_coordinates[2];
}

double NaturalCoord::getCoordinate(int index) const
{
    assert(index >= 0 && index < static_cast<int>(m_coordinates.size()));
    return (index < static_cast<int>(m_coordinates.size())) ? m_coordinates[index] : 0.0;
}

// ============================================================================
// Mutators
// ============================================================================

void NaturalCoord::setCoordinates(const std::vector<double>& coordinates)
{
    m_coordinates = coordinates;
}

void NaturalCoord::setCoordinates(double r, double s, double t)
{
    m_coordinates.resize(3);
    m_coordinates[0] = r;
    m_coordinates[1] = s;
    m_coordinates[2] = t;
}

void NaturalCoord::setCoordinates(double r, double s)
{
    m_coordinates.resize(2);
    m_coordinates[0] = r;
    m_coordinates[1] = s;
}

