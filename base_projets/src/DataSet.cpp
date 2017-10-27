#include <cstdlib>
#include <iostream>
#include <math.h>
#include <algorithm>

#include <DataSet.h>

using namespace std;

DataSet::DataSet(const char* filename)
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

	m_distances = vector<vector<double> > ();
	m_distances.resize(nb_points);

	// computing distances
	for (int i = 0; i < m_N; i++) {
		m_distances[i].resize(nb_points);
		glm::vec3 x1 = m_points[i];
		for (int j = 0; j < m_N; j++) {
			glm::vec3 x2 = m_points[j];
			double distance = sqrt((x2[0] - x1[0]) * (x2[0] - x1[0]) + 
								   (x2[1] - x1[1]) * (x2[1] - x1[1]) +
								   (x2[2] - x1[2]) * (x2[2] - x1[2]));
			m_distances[i][j] = distance;
		}
		std::sort(m_distances[i].begin(), m_distances[i].end());   
	}

    fclose(file);
}

Plane DataSet::ComputeTangentPlanes() {
}
