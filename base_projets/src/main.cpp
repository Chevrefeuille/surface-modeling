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
    DataSet ds("../data/test.data");
    ds.ComputeTangentPlanes();

    // Mesh creation from data set and iso function
    ImplicitFunction f();
    //Mesh m(ds, f);
    //m.postProcess();

    return EXIT_SUCCESS;
}
