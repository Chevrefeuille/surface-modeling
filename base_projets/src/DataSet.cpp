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
};

bool operator<(const PLANE_AND_DISTANCE& a, const PLANE_AND_DISTANCE& b)
{
	return a.distance < b.distance;
}

DataSet::DataSet(const char* filename) :
	m_K(2)
{
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

	// for (int i = 0; i < m_N; i++) {
	// 	std::cout << "row:" << i << std::endl;
	// 	for (int j = 0; j < m_N; j++) {
	// 		std::cout << m_distances[i][j] << std::endl;
	// 	}
	// }

    fclose(file);
}


std::vector<glm::vec3> DataSet::ComputeNhbd(glm::vec3 x) {
	std::vector<POINT_AND_DISTANCE> points_and_distances(m_N);
	// std::cout << x[0] << ", " << x[1] << ", " << x[2] << std::endl;
	// computing distances
	for (int i = 0; i < m_N; i++) {
		glm::vec3 x2 = m_points[i];
		double distance = sqrt((x2[0] - x[0]) * (x2[0] - x[0]) +
		(x2[1] - x[1]) * (x2[1] - x[1]) +
		(x2[2] - x[2]) * (x2[2] - x[2]));

		//std::cout << x2[0] << ", " << x2[1] << ", " << x2[2] << " ==> " << distance << std::endl;
		points_and_distances[i].point = x2;
		points_and_distances[i].distance = distance;

	}
	std::sort(points_and_distances.begin(), points_and_distances.end());
	// extracting k nearest neighbors
	std::vector<glm::vec3> Nhbd(m_K);
	for (int i = 1; i < m_K + 1; i++) { // skiping first element which is the point itself (at distance 0)
		Nhbd[i-1] = points_and_distances[i].point;
		// std::cout << "-" << Nhbd[i-1][0] << ", " << Nhbd[i-1][1] << ", " << Nhbd[i-1][2] << " ==> " << points_and_distances[i].distance << std::endl;
	}
	return Nhbd;
}

glm::vec3 DataSet::ComputeCentroid(std::vector<glm::vec3> points) {
	glm::vec3 o;
	for (int i = 0; i < m_K + 1; i++) {
		glm::vec3 x2 = points[i];
		// std::cout << x2[0] << ", " << x2[1] << ", " << x2[2] << std::endl;
		o += x2;
	}
	o /= (m_K + 1);
	//std::cout << o[0] << ", " << o[1] << ", " << o[2]  << std::endl;
	return o;
}

glm::vec3 DataSet::ComputeTangent(std::vector<glm::vec3> points, glm::vec3 o) {
	// compute covariance matrix
	Eigen::MatrixXd CV(3, 3);
	CV << 0.0,0.0,0.0,
	0.0,0.0,0.0,
	0.0,0.0,0.0;
    for (int k = 0; k < m_N; k++) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				CV(i, j) += (points[k][i] - o[i]) * (points[k][j] - o[j]);
			}
		}
	}
	//std::cout << CV << std::endl;
	// compute eigenvalues of the covariance matrix
	Eigen::EigenSolver<Eigen::MatrixXd> es(CV);
	//std::cout << "The eigenvalues of CV are:" << std::endl << es.eigenvalues() << std::endl;
	//std::cout << "The matrix of eigenvectors, V, is:" << std::endl << es.eigenvectors() << std::endl << std::endl;

	// extract the third eigen vector
	Eigen::VectorXcd v = es.eigenvectors().col(2);
	glm::vec3 n(3);
	for (int i = 0; i < 3; i++) {
		n[i] = v(i).real();
	}
	return n;
}


Plane DataSet::ComputeTangentPlanes() {
	m_tangentPlanes.resize(m_N);
	std::vector<glm::vec3> Nhbd(m_K + 1);
	for (int i = 0; i < m_N; i++) {
		// std::cout << m_points[i][0] << ", " << m_points[i][1] << ", " << m_points[i][2] << " => " << std::endl;
		Nhbd = ComputeNhbd(m_points[i]);
        Nhbd.push_back(m_points[i]);
		glm::vec3 o = ComputeCentroid(Nhbd);
		glm::vec3 n = ComputeTangent(Nhbd, o);
		// std::cout << "Centroid : " << o[0] << ", " << o[1] << ", " << o[2] << std::endl;
		// std::cout << "Tangent : " << n[0] << ", " << n[1] << ", " << n[2] << std::endl;
		// std::cout << "-------" << std::endl;
		Plane tangentPlane(o, n);
		m_tangentPlanes[i] = tangentPlane;
	}
}


void DataSet::ComputeEMST() {
	for (int i = 0; i < m_N; i++) {
		Plane pi = m_tangentPlanes[i];
		glm::vec3 vi = pi.getCenter();
		m_graph.addVertex(pi);
		for (int j = 0; j < m_N; j++) {
			Plane pj = m_tangentPlanes[j];
			m_graph.addVertex(pj);
			glm::vec3 vj = pj.getCenter();
			double distance = sqrt((vj[0] - vi[0]) * (vj[0] - vi[0]) +
			(vj[1] - vi[1]) * (vj[1] - vi[1]) +
			(vj[2] - vi[2]) * (vj[2] - vi[2]));
			m_graph.addEdge(pi, pj, distance);
			//std::cout << distance << std::endl;
		}
	}
    m_graph.printGraph();
    m_graph.computeMSTwithPrim();
}

std::vector<Plane> DataSet::ComputeKNeigbors(Plane p) {
	std::vector<PLANE_AND_DISTANCE> planes_and_distances(m_N - 1);

	// computing distances
	glm::vec3 x = p.getCenter();
	//std::cout << x[0] << ", " << x[1] << ", " << x[2] << std::endl;
    int k = 0;
	for (int i = 0; i < m_N; i++) {
		Plane p2 = m_tangentPlanes[i];
		glm::vec3 x2 = p2.getCenter();
        if (x2 != x) {
            double distance = sqrt((x2[0] - x[0]) * (x2[0] - x[0]) +
    		(x2[1] - x[1]) * (x2[1] - x[1]) +
    		(x2[2] - x[2]) * (x2[2] - x[2]));

    		//std::cout << x2[0] << ", " << x2[1] << ", " << x2[2] << " ==> " << distance << std::endl;
    		planes_and_distances[k].plane = p2;
            //std::cout << "->" << planes_and_distances[i].plane.getCenter()[0] << ", " << planes_and_distances[i].plane.getCenter()[1] << ", " << planes_and_distances[i].plane.getCenter()[2] << " ==> " << distance << std::endl;
    		planes_and_distances[k].distance = distance;
            k++;
        }
	}
    std::sort(planes_and_distances.begin(), planes_and_distances.end());
	// extracting k nearest neighbors
	std::vector<Plane> KNeigbors(m_K);
	for (int i = 0; i < m_K ; i++) {
		KNeigbors[i] = planes_and_distances[i].plane;
		glm::vec3 x2 = KNeigbors[i].getCenter();
		//std::cout << "*" << x2[0] << ", " << x2[1] << ", " << x2[2] << " ==> " << planes_and_distances[i].distance << std::endl;
	}
	return KNeigbors;
}

void DataSet::AddKNeighborsEdges() {
	std::vector<Plane> KNeighbors(m_K);
	for (int i = 0; i < m_N; i++) {
		Plane pi = m_tangentPlanes[i];
		glm::vec3 vi = pi.getCenter();
        //std::cout << "point " << vi.x << ", " << vi.y << ", " << vi.z << " neigbors :" << std::endl;
		KNeighbors = ComputeKNeigbors(pi);
		for (int j = 0; j < m_K; j++) {
			Plane pj = KNeighbors[j];
            glm::vec3 vj = pj.getCenter();
            //std::cout << vj.x << ", " << vj.y << ", " << vj.z << std::endl;
			m_graph.addEdge(pi, pj, 0);
		}
	}
	m_graph.printGraph();
}

void DataSet::AssignCostOnEdges() {
    for (vertices_map::iterator it_vertices = m_graph.work->begin(); it_vertices != m_graph.work->end(); ++it_vertices) {
        VertexG* u = it_vertices->second;
        for (std::vector<ve>::iterator it_neighbors = u->adj.begin() ; it_neighbors != u->adj.end(); ++it_neighbors) {
            VertexG* v = it_neighbors->second;
            it_neighbors->first = 1 - abs(glm::dot(u->plane.getNormal(), v->plane.getNormal()));
        }
    }
    m_graph.printGraph();
}

void DataSet::AssignTangentPlanesOrientation() {
    m_graph.computeMSTwithPrim();
    m_graph.DFS(m_graph.work->begin()->second, NULL);
}
