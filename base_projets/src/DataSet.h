#ifndef DATASET_H
#define DATASET_H

#include <glm/glm.hpp>
#include <vector>
#include <string>

#include "Plane.h"

using namespace glm;

class DataSet
{
public:
    // Constructors
    DataSet(){}                        /// Empty constructor
    DataSet(const char* filename);     /// Imports a mesh from a data file
    Plane ComputeTangentPlanes();

protected:
    // Attributes
    int m_N;                                    /// Number of points in the set
    std::vector<glm::vec3> m_points;             /// Container for the vertices positions
    std::vector<Plane> m_tangentPlanes;
    std::vector<std::vector<double> > m_distances;
};


#endif // DATASET_H
