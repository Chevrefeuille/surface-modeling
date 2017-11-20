#include "Graph.h"

void Graph::addVertex(const glm::vec3& coordinates)
{
    vertices_map::iterator itr = work.find(coordinates);
    if (itr == work.end())
    {
        vertex *v;
        v = new vertex(coordinates);
        work[coordinates] = v;
        return;
    }
    std::cout << "\nVertex already exists!";
}

void Graph::addEdge(const glm::vec3& from, const glm::vec3& to, double cost)
{
    vertex *f = (work.find(from)->second);
    vertex *t = (work.find(to)->second);
    std::pair<double, vertex *> edge = std::make_pair(cost, t);
    f->adj.push_back(edge);
}
