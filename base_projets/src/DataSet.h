#ifndef DATASET_H
#define DATASET_H

#include <glm/glm.hpp>
#include <vector>
#include <string>
#include <tr1/unordered_map>

#include "Plane.h"

struct KeyFuncs;
struct Node;

typedef std::tr1::unordered_map<Node, double, KeyFuncs, KeyFuncs> adjacency_map;

struct Node
{
    glm::vec3 coordinates;
    adjacency_map adj;
};


// struct KeyFuncs
// {
//     size_t operator()(const Node& k)const
//     {
//         return std::tr1::hash<double>()(k.coordinates.x) ^ std::tr1::hash<double>()(k.coordinates.y) ^ std::tr1::hash<double>()(k.coordinates.z);
//     }
//
//     bool operator()(const Node& a, const Node& b)const
//     {
//             return a.coordinates.x == b.coordinates.x && a.coordinates.y == b.coordinates.y && a.coordinates.z == b.coordinates.z;
//     }
// };


class DataSet {

public:
    // Constructors
    DataSet(){}                        /// Empty constructor
    DataSet(const char* filename);     /// Imports a mesh from a data file
    Plane ComputeTangentPlanes();

    // Accessors
    int nbPoints() const {return m_N;}
    double minX() const {return min_X;};
    double minY() const {return min_Y;};
    double minZ() const {return min_Z;};
    double maxX() const {return max_X;};
    double maxY() const {return max_Y;};
    double maxZ() const {return max_Z;};

protected:
    // Attributes
    int m_N;                                    /// Number of points in the set
    int m_K;
    std::vector<glm::vec3> m_points;             /// Container for the vertices positions
    std::vector<Plane> m_tangentPlanes;
    // min/max coordinate values in each direction
    double min_X, min_Y, min_Z, max_X, max_Y, max_Z;

    std::vector<glm::vec3> ComputeNhbd(glm::vec3 x);
    glm::vec3 ComputeCentroid(std::vector<glm::vec3>);
    glm::vec3 ComputeTangent(std::vector<glm::vec3>, glm::vec3 o);

};


#endif // DATASET_H
