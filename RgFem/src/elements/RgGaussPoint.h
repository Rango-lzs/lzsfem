#pragma once

#include "datastructure/Vector3d.h"
#include <vector>

namespace RgFem {

/**
 * @brief Represents a Gauss integration point in the finite element method
 * 
 * This class stores the natural coordinates and weight of a Gauss point
 * used in numerical integration. It can be used in 1D, 2D or 3D integration.
 */
class RgGaussPoint {
public:
    /**
     * @brief Default constructor
     */
    RgGaussPoint();

    /**
     * @brief Constructor with coordinates and weight
     * @param coordinates Natural coordinates of the Gauss point
     * @param weight Integration weight of the Gauss point
     */
    RgGaussPoint(const std::vector<double>& coordinates, double weight);

    /**
     * @brief Constructor with individual coordinate components and weight
     * @param r First natural coordinate
     * @param s Second natural coordinate
     * @param t Third natural coordinate
     * @param weight Integration weight
     */
    RgGaussPoint(double r, double s, double t, double weight);

    /**
     * @brief Constructor for 2D Gauss point
     * @param r First natural coordinate
     * @param s Second natural coordinate
     * @param weight Integration weight
     */
    RgGaussPoint(double r, double s, double weight);

    /**
     * @brief Constructor for 1D Gauss point
     * @param r First natural coordinate
     * @param weight Integration weight
     */
    RgGaussPoint(double r, double weight);

    /**
     * @brief Copy constructor
     */
    RgGaussPoint(const RgGaussPoint& other) = default;

    /**
     * @brief Assignment operator
     */
    RgGaussPoint& operator=(const RgGaussPoint& other) = default;

    /**
     * @brief Destructor
     */
    ~RgGaussPoint() = default;

    // Dimension info
    /**
     * @brief Get the dimension of the Gauss point (1, 2, or 3)
     * @return Dimension of the Gauss point
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

    // Weight accessor
    /**
     * @brief Get the integration weight
     * @return Integration weight
     */
    double getWeight() const;

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

    /**
     * @brief Set the integration weight
     * @param weight Integration weight
     */
    void setWeight(double weight);

private:
    std::vector<double> m_coordinates; ///< Natural coordinates of the Gauss point
    double m_weight;                  ///< Integration weight
};

} // namespace RgFem