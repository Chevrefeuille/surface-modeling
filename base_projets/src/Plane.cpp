#include <iostream>
#include <math.h>

#include <Plane.h>

using namespace std;

Plane::Plane(glm::vec3 center, glm::vec3 normal) :
    m_center(center), m_normal(normal)
{
}
