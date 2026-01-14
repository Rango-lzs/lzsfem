#pragma once

#include <vector>

// Forward declaration
class RgGaussPoint;

/**
 * @brief Represents a natural coordinate in the finite element method
 *
 * This class stores the natural coordinates used in finite element shape functions
 * and can be used in 1D, 2D or 3D elements. It can be constructed from a RgGaussPoint
 * by extracting only the coordinate values, discarding the weight.
 */
class NaturalCoord
{
public:
    /**
     * @brief Default constructor
     */
    NaturalCoord();

    /**
     * @brief Constructor with coordinates
     * @param coordinates Natural coordinates
     */
    NaturalCoord(const std::vector<double>& coordinates);

    /**
     * @brief Constructor with individual coordinate components
     * @param r First natural coordinate
     * @param s Second natural coordinate
     * @param t Third natural coordinate
     */
    NaturalCoord(double r, double s, double t);

    /**
     * @brief Constructor for 2D natural coordinate
     * @param r First natural coordinate
     * @param s Second natural coordinate
     */
    NaturalCoord(double r, double s);

    /**
     * @brief Constructor for 1D natural coordinate
     * @param r First natural coordinate
     */
    NaturalCoord(double r);

    /**
     * @brief Constructor from RgGaussPoint, extracting only coordinates
     * @param gaussPoint The Gauss point to extract coordinates from
     */
    NaturalCoord(const RgGaussPoint& gaussPoint);

    /**
     * @brief Copy constructor
     */
    NaturalCoord(const NaturalCoord& other) = default;

    /**
     * @brief Assignment operator
     */
    NaturalCoord& operator=(const NaturalCoord& other) = default;

    /**
     * @brief Destructor
     */
    ~NaturalCoord() = default;

    // Dimension info
    /**
     * @brief Get the dimension of the natural coordinate (1, 2, or 3)
     * @return Dimension of the natural coordinate
     */
    int getDimension() const;

    // Coordinate accessors
    /**
     * @brief Get the natural coordinates
     * @return Vector of coordinates
     */
    const std::vector<double>& getCoordinates() const;

    /**
     * @brief Get the r-coordinate (first natural coordinate)
     * @return r-coordinate value
     */
    double getR() const;

    /**
     * @brief Get the s-coordinate (second natural coordinate)
     * @return s-coordinate value
     */
    double getS() const;

    /**
     * @brief Get the t-coordinate (third natural coordinate)
     * @return t-coordinate value
     */
    double getT() const;

    /**
     * @brief Get coordinate by index
     * @param index Index of coordinate (0 for r, 1 for s, 2 for t)
     * @return Coordinate value
     */
    double getCoordinate(int index) const;

    // Coordinate mutators
    /**
     * @brief Set the natural coordinates
     * @param coordinates Vector of coordinates
     */
    void setCoordinates(const std::vector<double>& coordinates);

    /**
     * @brief Set individual coordinate components
     * @param r First natural coordinate
     * @param s Second natural coordinate
     * @param t Third natural coordinate
     */
    void setCoordinates(double r, double s, double t);

    /**
     * @brief Set 2D coordinates
     * @param r First natural coordinate
     * @param s Second natural coordinate
     */
    void setCoordinates(double r, double s);

private:
    std::vector<double> m_coordinates;  ///< Natural coordinates
};
