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

#define INF std::numeric_limits<double>::infinity()
struct vertex;
typedef std::pair<double, vertex*> ve;

struct KeyFuncs {
    size_t operator()(const glm::vec3& k)const {
        return std::tr1::hash<double>()(k.x) ^ std::tr1::hash<double>()(k.y) ^ std::tr1::hash<double>()(k.z);
    }

    bool operator()(const glm::vec3& a, const glm::vec3& b)const {
        return a.x == b.x && a.y == b.y && a.z == b.z;
    }
};

struct vertex {
    std::vector<ve> adj; //cost of edge, destination vertex
    glm::vec3 coordinates;
    double cost; //useful for Prim's algorithm
    vertex* next; //useful for Prim's algorithm
    vertex(glm::vec3 c, double _cost, vertex* _next) : coordinates(c), cost(_cost), next(_next) {}
};

struct LessThanByCost
{
  bool operator()(const vertex* lhs_vertex, const vertex* rhs_vertex) const
  {
    return lhs_vertex->cost < rhs_vertex->cost;
  }
};

typedef std::tr1::unordered_map<glm::vec3, vertex*, KeyFuncs> vertices_map;

class Graph
{
public:
    Graph();
    vertices_map* work;
    void addVertex(const glm::vec3& coordinates);
    void addEdge(const glm::vec3& from, const glm::vec3& to, double cost);
    void printGraph();
    void prim(vertex& start);
};

#endif // GRAPH_H
