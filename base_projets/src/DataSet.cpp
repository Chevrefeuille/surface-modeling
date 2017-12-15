#include <cstdlib>
#include <iostream>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "DataSet.h"

struct POINT_AND_DISTANCE {
    glm::vec3 point;
	double distance;
};

bool operator<(const POINT_AND_DISTANCE& a, const POINT_AND_DISTANCE& b)
{
	return a.distance < b.distance;
}

struct PLANE_AND_DISTANCE {
    Plane plane;
	double distance;
    PLANE_AND_DISTANCE(Plane _plane, double _distance) : plane(_plane), distance(_distance) {}
};

bool operator<(const PLANE_AND_DISTANCE& a, const PLANE_AND_DISTANCE& b)
{
	return a.distance < b.distance;
}

DataSet::DataSet(const char* filename, int K, float rho) :
	m_K(K), m_rho(rho)
{

    m_graph = new Graph();
	FILE *file;
	int error;
    int nb_points;

    min_X=0.0; min_Y=0.0; min_Z=0.0; max_X=0.0; max_Y=0.0; max_Z=0.0;

	if ((file = fopen(filename, "r")) == NULL) {
		std::cout << "Unable to read : " << filename << std::endl;
	}

	// create mesh
    m_points = std::vector<glm::vec3>();

	error = fscanf(file, "%d\n", &(nb_points));
	if(error == EOF) {
		std::cout << "Unable to read : " << filename << std::endl;
	}

	m_N = nb_points;
    m_points.resize(nb_points);

	// reading points
	for (int i = 0; i < nb_points; ++i) {
        error = fscanf(file, "%f %f %f\n", &(m_points[i][0]), &(m_points[i][1]), &(m_points[i][2]));
		if(error == EOF) {
			std::cout << "Unable to read points of : " << filename << std::endl;
			exit(EXIT_FAILURE);
		}

		// Setting min/max coord
		if (m_points[i][0]<min_X) {min_X=m_points[i][0];}
		if (m_points[i][0]>max_X) {max_X=m_points[i][0];}
		if (m_points[i][1]<min_Y) {min_Y=m_points[i][1];}
		if (m_points[i][1]>max_Y) {max_Y=m_points[i][1];}
		if (m_points[i][2]<min_Z) {min_Z=m_points[i][2];}
		if (m_points[i][2]>max_Z) {max_Z=m_points[i][2];}

	}
    fclose(file);
}

DataSet::~DataSet() {
    delete m_graph;
}

void DataSet::ComputeNhbd(std::vector<glm::vec3>& Nhbd, glm::vec3 x) {
	std::vector<POINT_AND_DISTANCE> points_and_distances(m_N);
	// computing distances
	for (int i = 0; i < m_N; i++) {
		glm::vec3 x2 = m_points[i];
		double distance = glm::distance(x2, x);

		points_and_distances[i].point = x2;
		points_and_distances[i].distance = distance;

	}
	std::sort(points_and_distances.begin(), points_and_distances.end());
	// extracting k nearest neighbors
	for (int i = 1; i < m_K + 1; i++) { // skiping first element which is the point itself (at distance 0)
        Nhbd.push_back(points_and_distances[i].point);
	}
	//return Nhbd;
}

glm::vec3 DataSet::ComputeCentroid(std::vector<glm::vec3> points) {
	glm::vec3 o;
    int N = points.size();
	for (int i = 0; i < N; i++) {
		glm::vec3 x2 = points[i];
		o += x2;
	}
	o /= (m_K + 1);
	return o;
}

glm::vec3 DataSet::ComputeTangent(std::vector<glm::vec3> points, glm::vec3 o) {
	// compute covariance matrix
	Eigen::MatrixXd CV(3, 3);
	CV << 0.0,0.0,0.0,
	0.0,0.0,0.0,
	0.0,0.0,0.0;
	int N = points.size();
    for (int k = 0; k < N ; k++) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				CV(i, j) += (points[k][i] - o[i]) * (points[k][j] - o[j]);
			}
		}
	}
	// compute eigenvalues of the covariance matrix
	Eigen::EigenSolver<Eigen::MatrixXd> es(CV);
    int min_i = 0;
    float min_eigenvalue = std::abs(es.eigenvalues()[0]);
    for (int i = 1; i < 3; i++) {
        float eigenvalue = std::abs(es.eigenvalues()[i]);
        if (eigenvalue < min_eigenvalue) {
            min_eigenvalue = eigenvalue;
            min_i = i;
        }
    }

	// extract the third eigen vector
	Eigen::VectorXcd v = es.eigenvectors().col(min_i);
	glm::vec3 n(3);
	for (int i = 0; i < 3; i++) {
		n[i] = v(i).real();
	}
	return n;
}


Plane DataSet::ComputeTangentPlanes() {
	m_tangentPlanes.resize(m_N);
	for (int i = 0; i < m_N; i++) {
    	std::vector<glm::vec3> Nhbd;
        ComputeNhbd(Nhbd, m_points[i]);
        Nhbd.push_back(m_points[i]);
		glm::vec3 o = ComputeCentroid(Nhbd);
		glm::vec3 n = ComputeTangent(Nhbd, o);
		
		Plane tangentPlane(o, n);
		m_tangentPlanes[i] = tangentPlane;
	}
}


void DataSet::ComputeEMST() {
	for (int i = 0; i < m_N; i++) {
		Plane pi = m_tangentPlanes[i];
		glm::vec3 vi = pi.getCenter();
		m_graph->addVertex(pi);
		for (int j = 0; j < m_N; j++) {
			Plane pj = m_tangentPlanes[j];
			m_graph->addVertex(pj);
			glm::vec3 vj = pj.getCenter();
			double distance = glm::distance(vi, vj);
			m_graph->addEdge(pi, pj, distance);
		}
	}
    m_graph->computeMSTwithPrim();

}

void DataSet::ComputeKNeigbors(std::vector<Plane>& KNeighbors, Plane p) {

    KNeighbors.clear();
	std::vector<PLANE_AND_DISTANCE> planes_and_distances;

	// computing distances
	glm::vec3 x = p.getCenter();

	for (int i = 0; i < m_N; i++) {
		Plane p2 = m_tangentPlanes[i];
		glm::vec3 x2 = p2.getCenter();
        if (x2 != x) {
            double distance = glm::distance(x, x2);
            PLANE_AND_DISTANCE pd(p2, distance);
            planes_and_distances.push_back(pd);
        }
	}
    std::sort(planes_and_distances.begin(), planes_and_distances.end());
	// extracting k nearest neighbors
	for (int i = 0; i < m_K ; i++) {
        KNeighbors[i] = planes_and_distances[i].plane;
	}
}

void DataSet::AddKNeighborsEdges() {
    std::vector<Plane> KNeighbors (m_K);
	for (int i = 0; i < m_N; i++) {
		Plane pi = m_tangentPlanes[i];
		glm::vec3 vi = pi.getCenter();
		ComputeKNeigbors(KNeighbors, pi);
		for (int j = 0; j < m_K; j++) {
			Plane pj = KNeighbors[j];
            glm::vec3 vj = pj.getCenter();
            double distance = glm::distance(vi, vj);
            m_graph->addEdge(pi, pj, distance);
		}
	}
}

void DataSet::AssignCostOnEdges() {
    for (vertices_map::iterator it_vertices = m_graph->work->begin(); it_vertices != m_graph->work->end(); ++it_vertices) {
        VertexG* u = it_vertices->second;
        for (std::vector<ve>::iterator it_neighbors = u->adj.begin() ; it_neighbors != u->adj.end(); ++it_neighbors) {
            VertexG* v = it_neighbors->second;
            it_neighbors->first = 1 - abs(glm::dot(u->plane.getNormal(), v->plane.getNormal()));
        }
    }
}

void DataSet::AssignTangentPlanesOrientation() {
    m_graph->computeMSTwithPrim();
    m_graph->DFS(m_graph->maxZCenter, NULL);
	m_graph->writingPlanesIntoFile();
}
