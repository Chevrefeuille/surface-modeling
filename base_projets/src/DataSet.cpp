#include <cstdlib>
#include <iostream>
#include <math.h>

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

    m_points.resize(nb_points);

	// reading points
	for (int i = 0; i < nb_points; ++i) {
        error = fscanf(file, "%f %f %f\n", &(m_points[i][0]), &(m_points[i][1]), &(m_points[i][2]));
		if(error == EOF) {
			std::cout << "Unable to read points of : " << filename << std::endl;
			exit(EXIT_FAILURE);
		}
	}

    fclose(file);
}
