#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <map>
#include <string>

using namespace std;

struct vertex {
    typedef pair<int, vertex*> ve;
    vector<ve> adj; //cost of edge, destination vertex
    string name;
    vertex(string s) : name(s) {}
};

class Graph
{
public:
    typedef map<string, vertex *> vmap;
    vmap work;
    void addVertex(const string&);
    void addEdge(const string& from, const string& to, double cost);
};

#endif // GRAPH_H