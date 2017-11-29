#include <iostream>
#include <math.h>
#include <sstream>


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

std::string Plane::printPlane() const {
    std::string s;
    std::stringstream out;
    out << "\nCenter: " << m_center.x << " " << m_center.y << " " << m_center.z << "\nNormal: " << m_normal.x << " " << m_normal.y << " " << m_normal.z;
    s = out.str();
    return s;
}
