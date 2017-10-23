#ifndef OBJECT_H
#define OBJECT_H

#include <GL/glew.h>
#include <glm/glm.hpp>

#include <AbstractMesh.h>

#define DISPLAY_DEBUG_INFO  true


class Object
{
public:
    Object();
    ~Object();

    void Draw(const glm::mat4& projection_matrix, const glm::mat4& view_matrix, const GLuint PmatrixID, const GLuint VmatrixID) const;

    void GenBuffers();
    void UpdateGeometryBuffers();
    void UpdateElementsBuffer();
    void UpdateBuffers();

    void SetMesh(AbstractMesh *mesh);
    void SetShader(const GLuint programID);
    void UpdateAttributeLocations();

public:

    AbstractMesh* m_mesh;
    unsigned int m_id;

    GLuint m_programID;
    GLuint m_positionID;
    GLuint m_normalID;
    GLuint m_colorID;

    GLuint m_vertexBufferID;
    GLuint m_normalBufferID;
    GLuint m_colorBufferID;
    GLuint m_elementBufferID;

    static unsigned int s_id;
};

#endif
