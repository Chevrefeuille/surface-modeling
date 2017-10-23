#ifndef DATASET_H
#define DATASET_H

#include <glm/glm.hpp>
#include <vector>
#include <string>

using namespace glm;

class DataSet
{
public:
    // Constructors
    DataSet(){}                        /// Empty constructor
    DataSet(const char* filename);     /// Imports a mesh from a data file

protected:
    // Attributes
    std::vector<glm::vec3> m_points;             /// Container for the vertices positions
};


#endif // DATASET_H
