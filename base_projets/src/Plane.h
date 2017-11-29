#ifndef PLANE_H
#define PLANE_H

#include <glm/glm.hpp>
#include <vector>

using namespace glm;

class Plane
{
public:
    // Constructors
    Plane(){}                                      /// Empty constructor
    Plane(glm::vec3 center, glm::vec3 normal);     /// Create plane from center and normal
    glm::vec3 getNormal() const;
    glm::vec3 setNormal(const glm::vec3& normal);
    glm::vec3 getCenter() const;
    std::string printPlane() const;

protected:
    // Attributes
    glm::vec3 m_center;             /// Center of the plane
    glm::vec3 m_normal;             /// Unit normal of the plane
};

#endif //PLANE_H
