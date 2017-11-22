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
struct vertex;
typedef std::pair<double, vertex*> ve;

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

struct vertex {
    std::vector<ve> adj; //cost of edge, destination vertex
    Plane plane;
    double cost; //used by Prim's algorithm
    vertex* prev; //used by Prim's algorithm
    bool isInMST; //used by Prim's algorithm
    bool isMarked; //used by DFS algorithm
    vertex(Plane p, double _cost, vertex* _prev, bool _isInMST, bool _isMarked) : plane(p), cost(_cost), prev(_prev), isInMST(_isInMST), isMarked(_isMarked) {}
};

struct GreaterThanByCost {
  bool operator()(const vertex* lhs_vertex, const vertex* rhs_vertex) const
  {
    return lhs_vertex->cost > rhs_vertex->cost;
  }
};

typedef std::tr1::unordered_map<Plane, vertex*, KeyFuncs, KeyFuncs> vertices_map;

class Graph
{
public:
    Graph();
    vertices_map* work;
    void addVertex(const Plane& plane);
    void addEdge(const Plane& from, const Plane& to, double cost);
    void printGraph();
    void computeMSTwithPrim();
    void DFS(vertex* curr, vertex* prev);
};

#endif // GRAPH_H
