#include "Graph.h"
#include <algorithm>
#include <iostream>
#include <fstream>

Graph::Graph() {
    work = new vertices_map();
    maxZCenter = NULL;
}

Graph::~Graph() {
    delete work;
}

void Graph::addVertex(const Plane& p)
{
    vertices_map::iterator search =  work->find(p);
    if (search ==  work->end())
    {
        //std::cout << "New Vertex (" << coordinates.x << ", " << coordinates.y << ", " << coordinates.z << ")" << std::endl;
        VertexG *v;
        v = new VertexG(p, INF, false, NULL, false);
        (*work)[p] = v;

        if (maxZCenter == NULL || v->plane.getCenter().z > maxZCenter->plane.getCenter().z) {
            maxZCenter = v;
        }
    }
}

void Graph::addEdge(const Plane& from, const Plane& to, double cost)
{
    VertexG *f =  work->find(from)->second;
    VertexG *t =  work->find(to)->second;
    std::pair<double, VertexG *> edge = std::make_pair(cost, t);
    if (std::find(f->adj.begin(), f->adj.end(), edge) == f->adj.end()) {
        f->adj.push_back(edge);
        //std::cout << "Adding edge" << std::endl;
    } else {
        //std::cout << "Edge already exists" << std::endl;
    }
}

void Graph::computeMSTwithPrim() {
    for (vertices_map::iterator it = work->begin(); it != work->end(); ++it) {
        VertexG* u = it->second;
        u->isInMST = false;
        u->cost = INF;
        u->prev = NULL;
    }
    VertexG* start = work->begin()->second;
    std::priority_queue<VertexG*, std::vector<VertexG*>, GreaterThanByCost> queue;
    start->cost = 0;
    queue.push(start);
    while (!queue.empty()) {
        VertexG* t = queue.top();
        //std::cout << "priority : " << t->coordinates.x << " " << t->coordinates.y << " " << t->coordinates.z << std::endl;
        //std::cout << "priority : " << t-> cost << std::endl;
        queue.pop();
        t->isInMST = true;
        for (std::vector<ve>::iterator it = t->adj.begin() ; it != t->adj.end(); ++it) {
            VertexG* u = it->second;
            double dist = it->first;
            if (!u->isInMST && u->cost > dist) {
                //std::cout << t->plane.getCenter().x << " " << t->plane.getCenter().y << " " << t->plane.getCenter().z << std::endl;
                //std::cout << u->plane.getCenter().x << " " << u->plane.getCenter().y << " " << u->plane.getCenter().z << std::endl;
                //std::cout << it->first << std::endl;
                //std::cout << std::endl;
                u->cost = it->first;
                u->prev = t;
                queue.push(u);
            }
        }
    }
    Graph* MST = new Graph();
    for (vertices_map::iterator itr = work->begin(); itr !=  work->end(); ++itr) {
        if (itr->second->prev != NULL) {
            Plane u_plane = itr->second->plane;
            Plane v_plane = itr->second->prev->plane;
            double cost = itr->second->cost;
            MST->addVertex(u_plane);
            MST->addVertex(v_plane);
            MST->addEdge(u_plane, v_plane, cost);
            MST->addEdge(v_plane, u_plane, cost);
        }
    }
    *this = *MST;
}

void Graph::DFS(VertexG* curr, VertexG* prev) {
    glm::vec3 currNormal = curr->plane.getNormal();
    if (prev == NULL) {
        curr->plane.setNormal(abs(curr->plane.getNormal()));
    } else if (glm::dot(prev->plane.getNormal(), currNormal) < 0) {
        curr->plane.setNormal(- curr->plane.getNormal());
    }
    curr->isMarked = true;

    // glm::vec3 center = curr->plane.getCenter();
    // glm::vec3 normal = curr->plane.getNormal();
    // std::ofstream myfile ("../example.txt", std::ios_base::app);
    // if (myfile.is_open()) {
    //     myfile << center.x << " " << center.y << " " << center.z << " " << normal.x << " " << normal.y << " " << normal.z << std::endl;
    // }
    // myfile.close();

    //std::cout << "(" << c1.x << ", " << c1.y << ", " << c1.z << ")" << std::endl;
    for (std::vector<ve>::iterator it = curr->adj.begin() ; it != curr->adj.end(); ++it) {
        VertexG* u = it->second;
        if (!u->isMarked) {
            DFS(u, curr);
        }
    }
}

void Graph::printGraph() {
    std::cout << "\nGraph:" << std::endl;
    vertices_map::iterator itr;
    for (itr = work->begin(); itr !=  work->end(); itr++) {
        VertexG *v = itr->second;
        glm::vec3 c1 = v->plane.getCenter();
        int n_neighbors = (v->adj).size();
        for (int i = 0; i < n_neighbors; i++) {
            std::pair<double, VertexG*> ve = v->adj[i];
            double d = ve.first;
            glm::vec3 c2 = ve.second->plane.getCenter();
            std::cout << "(" << c1.x << ", " << c1.y << ", " << c1.z << ") ---"
            << d << "---> (" << c2.x << ", " << c2.y << ", " << c2.z << ")" << std::endl;
        }
    }
}

void Graph::writingPlanesIntoFile() {
    vertices_map::iterator itr;
    std::ofstream myfile ("../example.txt", std::ios::out);
    if (myfile.is_open()){
        for (itr = work->begin(); itr !=  work->end(); itr++) {
            glm::vec3 center = itr->second->plane.getCenter();
            glm::vec3 normal = itr->second->plane.getNormal();
            myfile << center.x << " " << center.y << " " << center.z << " " << normal.x << " " << normal.y << " " << normal.z << std::endl;
         }
    }
    myfile.close();
}
