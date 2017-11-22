#include "Graph.h"
#include <algorithm>

Graph::Graph() {
  work = new vertices_map();
}

void Graph::addVertex(const Plane& p)
{
    vertices_map::iterator search = (*work).find(p);
    if (search == (*work).end())
    {
        //std::cout << "New Vertex (" << coordinates.x << ", " << coordinates.y << ", " << coordinates.z << ")" << std::endl;
        vertex *v;
        v = new vertex(p, INF, NULL, false, false);
        (*work)[p] = v;
        return;
    }
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
            glm::vec3 c1 = itr->second->plane.getCenter();
            glm::vec3 c2 = itr->second->prev->plane.getCenter();
            std::cout << "(" << c1.x << ", " << c1.y << ", " << c1.z << ") ---"
            << d << "---> " << "(" << c2.x << ", " << c2.y << ", " << c2.z << ")" << std::endl;
        }
    }
}

void Graph::DFS(vertex* curr, vertex* prev) {
    glm::vec3 currNormal = curr->plane.getNormal();
    if (prev != NULL &&  glm::dot(prev->plane.getNormal(), currNormal) < 0) {
        curr->plane.setNormal(prev->plane.getNormal());
    }
    curr->isMarked = true;
    glm::vec3 c1 = curr->plane.getCenter();
    std::cout << "(" << c1.x << ", " << c1.y << ", " << c1.z << ")" << std::endl;
    for (std::vector<ve>::iterator it = curr->adj.begin() ; it != curr->adj.end(); ++it) {
        vertex* u = it->second;
        if (!u->isMarked) {
            DFS(u, curr);
        }
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
