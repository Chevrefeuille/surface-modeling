#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <glm/glm.hpp>
#include <tr1/functional>
#include <tr1/unordered_map>
#include "Plane.h"


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
    typedef std::pair<double, vertex*> ve;
    std::vector<ve> adj; //cost of edge, destination vertex
    Plane plane;
    vertex(Plane p) : plane(p) {}
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

};

#endif // GRAPH_H
