#include "Graph.h"

Graph::Graph()
{
  work = new vertices_map();
}

void Graph::addVertex(const glm::vec3& coordinates)
{
    vertices_map::iterator search = (*work).find(coordinates);
    if (search == (*work).end())
    {
        //std::cout << "New Vertex (" << coordinates.x << ", " << coordinates.y << ", " << coordinates.z << ")" << std::endl;
        vertex *v;
        v = new vertex(coordinates);
        (*work)[coordinates] = v;
        return;
    }
    //std::cout << "Vertex (" << coordinates.x << ", " << coordinates.y << ", " << coordinates.z << ") already exists!" << std::endl;
}

void Graph::addEdge(const glm::vec3& from, const glm::vec3& to, double cost)
{
    vertex *f = (*work).find(from)->second;
    vertex *t = (*work).find(to)->second;
    std::pair<double, vertex *> edge = std::make_pair(cost, t);
    f->adj.push_back(edge);
}

void Graph::printGraph() {
    vertices_map::iterator itr;
    for (itr = (*work).begin(); itr != (*work).end(); itr++) {
        vertex *v = itr->second;
        glm::vec3 c1 = v->coordinates;
        int n_neighbors = (v->adj).size();
        for (int i = 0; i < n_neighbors; i++) {
            std::pair<double, vertex*> ve = v->adj[i];
            double d = ve.first;
            glm::vec3 c2 = ve.second->coordinates;
            std::cout << "(" << c1.x << ", " << c1.y << ", " << c1.z << ") ---" 
            << d << "--->" << "(" << c2.x << ", " << c2.y << ", " << c2.z << ")" << std::endl;
        }
    }
}