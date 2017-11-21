#include "Graph.h"

Graph::Graph() {
  work = new vertices_map();
}

void Graph::addVertex(const glm::vec3& coordinates) {
    vertices_map::iterator search = (*work).find(coordinates);
    if (search == (*work).end())
    {
        std::cout << "New Vertex (" << coordinates.x << ", " << coordinates.y << ", " << coordinates.z << ")" << std::endl;
        vertex *v;
        v = new vertex(coordinates, INF, NULL);
        (*work)[coordinates] = v;
        return;
    }
    std::cout << "Vertex (" << coordinates.x << ", " << coordinates.y << ", " << coordinates.z << ") already exists!" << std::endl;
}

void Graph::addEdge(const glm::vec3& from, const glm::vec3& to, double cost) {
    vertices_map::iterator itr = work->find(from);
    if (!(itr == work->end())) {
        glm::vec3 c = itr->second->coordinates;
        std::cout << c.y << std::endl;
        std::cout << "found" << std::endl;
    } else {
        std::cout << "not found" << std::endl;
    }
    vertex* f = work->find(from)->second;
    vertex* t = work->find(to)->second;
    std::cout << "(" << f->coordinates.x << ", " << f->coordinates.y << ", " << f->coordinates.z << ")" << std::endl;
    std::pair<double, vertex *> edge = std::make_pair(cost, t);
    f->adj.push_back(edge);
}

void Graph::prim(vertex& start) {
    std::priority_queue<vertex*, std::vector<vertex*>, LessThanByCost> queue;
    start.cost = 0;
    queue.push(&start);
    while (!queue.empty()) {
        vertex* t = queue.top();
        queue.pop();
        for (std::vector<ve>::iterator it = t->adj.begin() ; it != t->adj.end(); ++it) {
            vertex* u = it->second;
            if (u->cost > it->first) {
                u->cost = it->first;
                t->next = u;
                queue.push(u);
            }
        }
    }
    std::cout << "MST: " << std::endl;
    vertex* v = start.next;
    while (v != NULL) {
        std::cout << v->coordinates.x << " " << v->coordinates.y << v->coordinates.z << std::endl;
    }
}
