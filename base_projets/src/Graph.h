#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <glm/glm.hpp>
#include <tr1/functional>
#include <tr1/unordered_map>
#include <queue>
#include "Plane.h"

#define INF std::numeric_limits<double>::infinity()
struct Vertex;
typedef std::pair<double, Vertex*> ve;

struct KeyFuncs {
    size_t operator()(const Plane& p)const {
        glm::vec3 k = p.getCenter();
        return std::tr1::hash<double>()(k.x) ^ std::tr1::hash<double>()(k.y) ^ std::tr1::hash<double>()(k.z);
    }

    bool operator()(const Plane& p1, const Plane& p2)const {
        glm::vec3 a = p1.getCenter();
        glm::vec3 b = p2.getCenter();
        return a.x == b.x && a.y == b.y && a.z == b.z;
    }
};

struct Vertex {
    std::vector<ve> adj; //cost of edge, destination Vertex
    Plane plane;
    double cost; //used by Prim's algorithm
    bool isInMST; //used by Prim's algorithm
    bool isMarked; //used by DFS algorithm
    Vertex(Plane p, double _cost, bool _isInMST, bool _isMarked) : plane(p), cost(_cost), isInMST(_isInMST), isMarked(_isMarked) {}
};

struct GreaterThanByCost {
  bool operator()(const Vertex* lhs_Vertex, const Vertex* rhs_Vertex) const
  {
    return lhs_Vertex->cost > rhs_Vertex->cost;
  }
};

typedef std::tr1::unordered_map<Plane, Vertex*, KeyFuncs, KeyFuncs> vertices_map;

class Graph
{
public:
    Graph();
    vertices_map* work;
    void addVertex(const Plane& plane);
    void addEdge(const Plane& from, const Plane& to, double cost);
    void printGraph();
    void computeMSTwithPrim();
    void DFS(Vertex* curr, Vertex* prev);
    void writingPlanesIntoFile();
};

#endif // GRAPH_H
