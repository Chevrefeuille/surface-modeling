#include "Graph.h"

Graph::Graph() {
  work = new vertices_map();
}

void Graph::addVertex(const glm::vec3& coordinates) {
    vertices_map::iterator search = (*work).find(coordinates);
    if (search == (*work).end()) {
        //std::cout << "New Vertex (" << coordinates.x << ", " << coordinates.y << ", " << coordinates.z << ")" << std::endl;
        vertex *v;
        v = new vertex(coordinates, INF, NULL, false);
        (*work)[coordinates] = v;
        return;
    }
    //std::cout << "Vertex (" << coordinates.x << ", " << coordinates.y << ", " << coordinates.z << ") already exists!" << std::endl;
}

void Graph::addEdge(const glm::vec3& from, const glm::vec3& to, double cost) {
    vertex *f = (*work).find(from)->second;
    vertex *t = (*work).find(to)->second;
    std::pair<double, vertex *> edge = std::make_pair(cost, t);
    f->adj.push_back(edge);
}

void Graph::computeMSTwithPrim() {
    vertex* start = work->begin()->second;
    std::priority_queue<vertex*, std::vector<vertex*>, GreaterThanByCost> queue;
    start->cost = 0;
    queue.push(start);
    while (!queue.empty()) {
        vertex* t = queue.top();
        //std::cout << "priority : " << t->coordinates.x << " " << t->coordinates.y << " " << t->coordinates.z << std::endl;
        //std::cout << "priority : " << t-> cost << std::endl;
        queue.pop();
        t->isInMST = true;
        for (std::vector<ve>::iterator it = t->adj.begin() ; it != t->adj.end(); ++it) {
            vertex* u = it->second;
            if (!u->isInMST && u->cost > it->first) {
                //std::cout << t->coordinates.x << " " << t->coordinates.y << " " << t->coordinates.z << std::endl;
                //std::cout << u->coordinates.x << " " << u->coordinates.y << " " << u->coordinates.z << std::endl;
                //std::cout << std::endl;
                u->cost = it->first;
                u->prev = t;
                queue.push(u);
            }
        }
    }
    std::cout << "MST: " << std::endl;
    for (vertices_map::iterator itr = (*work).begin(); itr != (*work).end(); itr++) {
        if (itr->second->prev != NULL) {
            double d = itr->second->cost;
            glm::vec3 c1 = itr->second->coordinates;
            glm::vec3 c2 = itr->second->prev->coordinates;
            std::cout << "(" << c1.x << ", " << c1.y << ", " << c1.z << ") ---"
            << d << "---> " << "(" << c2.x << ", " << c2.y << ", " << c2.z << ")" << std::endl;
        }
    }
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
            << d << "---> " << "(" << c2.x << ", " << c2.y << ", " << c2.z << ")" << std::endl;
        }
    }
}
