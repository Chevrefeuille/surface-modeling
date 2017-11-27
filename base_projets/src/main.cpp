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
   //ds.AddKNeighborsEdges();
   //ds.AssignCostOnEdges();
   //ds.AssignTangentPlanesOrientation();

    // SphereFunction f(vec3(0,0,0), 1);
    // printf("valIn = %lf\n", f.Eval(vec3(0.99,0,0)));
    // printf("valExt = %lf\n", f.Eval(vec3(1.01,0,0)));
    // printf("valOn = %lf\n", f.Eval(vec3(1,0,0)));


    // /** Mesh creation from data set and iso function **/
    // //double minX = ds.minX(); double minY = ds.minY(); double minZ = ds.minZ();
    // //double maxX = ds.maxX(); double maxY = ds.maxY(); double maxZ = ds.maxZ();
    // double minX = -1.; double minY = -1.; double minZ = -1.;
    // double maxX = 1.; double maxY = 1.; double maxZ = 1.;
    // const unsigned int resX=100;
    // const unsigned int resY=100;
    // const unsigned int resZ=100;
    // minX-=2.*(maxX-minX)/resX; minY-=2.*(maxY-minY)/resY; minZ-=2.*(maxZ-minZ)/resZ;
    // maxX+=2.*(maxX-minX)/resX; maxY+=2.*(maxY-minY)/resY; maxZ+=2.*(maxZ-minZ)/resZ;

    // Mesh m(f, minX, maxX, minY, maxY, minZ, maxZ, resX, resY, resZ);
    // printf("Mesh Created with %i point positions and %i faces\n", m.NbVertices(), m.NbFaces());
    // Mesh m2 = m.postProcess(0.001);
    // printf("Edge Collapsing done\n");

    return EXIT_SUCCESS;
}
