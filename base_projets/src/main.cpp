#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <glm/glm.hpp>


#include "DataSet.h"

using namespace glm;
using namespace std;


int main()
{
    // Mesh creation
    DataSet ds("../data/sphere.data");
    ds.ComputeTangentPlanes();

    return EXIT_SUCCESS;
}
