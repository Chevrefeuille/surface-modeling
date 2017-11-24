#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <glm/glm.hpp>

#include "DataSet.h"
#include "Mesh.h"
#include "ImplicitFunction.h"


using namespace glm;
using namespace std;


int main()
{
    // Data Set creation
    DataSet ds("../data/sphere.data");
    ds.ComputeTangentPlanes();

    ds.ComputeEMST();
    ds.AddKNeighborsEdges();
    ds.AssignCostOnEdges();
    ds.AssignTangentPlanesOrientation();

    //ImplicitFunction f();

    // ** Mesh creation from data set and iso function ** //
	// double minX = ds.minX(); double minY = ds.minY(); double minZ = ds.minZ();
	// double maxX = ds.maxX(); double maxY = ds.maxY(); double maxZ = ds.maxZ();
	// const unsigned int resX=100;
	// const unsigned int resY=100;
	// const unsigned int resZ=100;
	// minX-=2*(maxX-minX)/resX; minY-=2*(maxY-minY)/resY; minZ-=2*(maxZ-minZ)/resZ;
	// maxX+=2*(maxX-minX)/resX; maxY+=2*(maxY-minY)/resY; maxZ+=2*(maxZ-minZ)/resZ;

    //Mesh m(f, minX, maxX, minY, maxY, minZ, maxZ, resX, resY, resZ);
    //Mesh m2 = m.postProcess();

    return EXIT_SUCCESS;
}
