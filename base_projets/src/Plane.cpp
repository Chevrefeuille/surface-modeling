#include <iostream>
#include <math.h>

#include <Plane.h>

using namespace std;

Plane::Plane(glm::vec3 center, glm::vec3 normal) :
    m_center(center), m_normal(normal)
{
}

glm::vec3 Plane::getCenter() const {
    return m_center;
}

glm::vec3 Plane::setNormal(const glm::vec3& normal) {
    m_normal.x = normal.x;
    m_normal.y = normal.y;
    m_normal.z = normal.z;
}

glm::vec3 Plane::getNormal() const {
    return m_normal;
}
