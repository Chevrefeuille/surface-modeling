#include <cstdlib>
#include <iostream>
#include <math.h>
#include <algorithm>

#include <DataSet.h>

using namespace std;

struct POINT_AND_DISTANCE {    
    glm::vec3 point;  
	double distance;  
};

bool operator<(const POINT_AND_DISTANCE& a, const POINT_AND_DISTANCE& b)
{
	return a.distance < b.distance;
}

DataSet::DataSet(const char* filename) :
	m_K(5)
{
	FILE *file;
	int error;
    int nb_points;

	if ((file = fopen(filename, "r")) == NULL) {
		std::cout << "Unable to read : " << filename << std::endl;
	}

	// create mesh
    m_points = vector<vec3>();

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
	}

	glm::vec3 x(0,0,0);
	ComputeNhbd(x);
	
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
	// computing distances
	for (int i = 0; i < m_N; i++) {
		glm::vec3 x2 = m_points[i];
		double distance = sqrt((x2[0] - x[0]) * (x2[0] - x[0]) + 
								(x2[1] - x[1]) * (x2[1] - x[1]) +
								(x2[2] - x[2]) * (x2[2] - x[2]));

								
		points_and_distances[i].point = x2;
		points_and_distances[i].distance = distance;
		
	}
	std::sort(points_and_distances.begin(), points_and_distances.end());
	std::vector<glm::vec3> Nhbd(m_K);		
	for (int i = 0; i < m_K; i++) {
		Nhbd[i] = points_and_distances[i].point;
	}
	return Nhbd;
}


Plane DataSet::ComputeTangentPlanes() {
}
