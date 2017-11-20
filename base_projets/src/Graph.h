#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <glm/glm.hpp>



struct KeyFuncs {
    size_t operator()(const glm::vec3& k)const {
        return std::hash<double>()(k.x) ^ std::hash<double>()(k.y) ^ std::hash<double>()(k.z);
    }

    bool operator()(const glm::vec3& a, const glm::vec3& b)const {
            return a.x == b.x && a.y == b.y && a.z == b.z;
    }
};

struct vertex {
    typedef std::pair<double, vertex*> ve;
    std::vector<ve> adj; //cost of edge, destination vertex
    glm::vec3 coordinates;
    vertex(glm::vec3 c) : coordinates(c) {}
};

typedef std::map<glm::vec3, vertex*, KeyFuncs> vertices_map;

class Graph
{
public:
    vertices_map work;
    void addVertex(const glm::vec3& coordinates);
    void addEdge(const glm::vec3& from, const glm::vec3& to, double cost);
};

#endif // GRAPH_H