#include "Graph.h"
#include <algorithm>

Graph::Graph()
{
  work = new vertices_map();
}

void Graph::addVertex(const Plane& p)
{
    vertices_map::iterator search = (*work).find(p);
    if (search == (*work).end())
    {
        //std::cout << "New Vertex (" << coordinates.x << ", " << coordinates.y << ", " << coordinates.z << ")" << std::endl;
        vertex *v;
        v = new vertex(p);
        (*work)[p] = v;
        return;
    }
    //std::cout << "Vertex (" << coordinates.x << ", " << coordinates.y << ", " << coordinates.z << ") already exists!" << std::endl;
}

void Graph::addEdge(const Plane& from, const Plane& to, double cost)
{
    vertex *f = (*work).find(from)->second;
    vertex *t = (*work).find(to)->second;
    std::pair<double, vertex *> edge = std::make_pair(cost, t);
    if (std::find(f->adj.begin(), f->adj.end(), edge) == f->adj.end()) {
        f->adj.push_back(edge);
        //std::cout << "Adding edge" << std::endl;
    } else {
        //std::cout << "Edge already exists" << std::endl;
    }
}

void Graph::printGraph() {
    vertices_map::iterator itr;
    for (itr = (*work).begin(); itr != (*work).end(); itr++) {
        vertex *v = itr->second;
        glm::vec3 c1 = v->plane.getCenter();
        int n_neighbors = (v->adj).size();
        for (int i = 0; i < n_neighbors; i++) {
            std::pair<double, vertex*> ve = v->adj[i];
            double d = ve.first;
            glm::vec3 c2 = ve.second->plane.getCenter();
            std::cout << "(" << c1.x << ", " << c1.y << ", " << c1.z << ") ---"
            << d << "---> (" << c2.x << ", " << c2.y << ", " << c2.z << ")" << std::endl;
        }
    }
}
