#include "Graph.h"

void Graph::addVertex(const glm::vec3& coordinates)
{
    auto search = work.find(coordinates);
    if (search == work.end())
    {
        std::cout << "New Vertex (" << coordinates.x << ", " << coordinates.y << ", " << coordinates.z << ")" << std::endl;
        vertex *v;
        v = new vertex(coordinates);
        work[coordinates] = v;
        return;
    }
    std::cout << "Vertex (" << coordinates.x << ", " << coordinates.y << ", " << coordinates.z << ") already exists!" << std::endl;
}

void Graph::addEdge(const glm::vec3& from, const glm::vec3& to, double cost)
{
    vertices_map::iterator itr = work.find(from);
    if (!(itr == work.end())) {
        glm::vec3 c = itr->second->coordinates;
        //std::cout << c.y << std::endl;
        //std::cout << "found" << std::endl;
    } else {
        //std::cout << "not found" << std::endl;
    }
    //vertex *f = (work.find(from)->second);
    //vertex *t = (work.find(to)->second);
    //std::cout << "(" << f->coordinates.x << ", " << f->coordinates.y << ", " << f->coordinates.z << ")" << std::endl;
    //std::pair<double, vertex *> edge = std::make_pair(cost, t);
    //f->adj.push_back(edge);
}
