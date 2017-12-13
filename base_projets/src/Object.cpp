#include "Object.h"

#include <iostream>

#include <GL/glew.h>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/ext.hpp>

#include "shader.h"

using namespace std;
using namespace glm;


unsigned int Object::s_id = 0;

Object::Object():
    m_id(s_id)
{
    s_id++;
}

Object::~Object()
{
    glDeleteBuffers(1, &m_vertexBufferID);
    glDeleteBuffers(1, &m_normalBufferID);
    glDeleteBuffers(1, &m_colorBufferID);
    glDeleteBuffers(1, &m_elementBufferID);
}

void Object::SetMesh(AbstractMesh *mesh)
{
    m_mesh = mesh;
    UpdateBuffers();
}

void Object::UpdateGeometryBuffers()
{
    glBindBuffer(GL_ARRAY_BUFFER, m_vertexBufferID);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vec3) * m_mesh->NbVertices(), m_mesh->positions_array().data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, m_normalBufferID);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vec3) * m_mesh->NbVertices(), m_mesh->normals_array().data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, m_colorBufferID);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vec3) * m_mesh->NbVertices(), m_mesh->colors_array().data(), GL_STATIC_DRAW);
}

void Object::UpdateElementsBuffer()
{
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_elementBufferID);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(glm::uint) * 3 * m_mesh->NbFaces(), m_mesh->faces_array().data(), GL_STATIC_DRAW);
}

void Object::UpdateBuffers()
{
    UpdateGeometryBuffers();
    UpdateElementsBuffer();
}




void Object::GenBuffers()
{
    glGenBuffers(1, &m_vertexBufferID);
    if(DISPLAY_DEBUG_INFO)
        cout << "  vertexBufferID = " << m_vertexBufferID << endl;

    glGenBuffers(1, &m_normalBufferID);
    if(DISPLAY_DEBUG_INFO)
        cout << "  normalBufferID = " << m_normalBufferID << endl;

    glGenBuffers(1, &m_colorBufferID);
    if(DISPLAY_DEBUG_INFO)
        cout << "  m_colorBufferID = " << m_colorBufferID << endl;

    glGenBuffers(1, &m_elementBufferID);
    if(DISPLAY_DEBUG_INFO)
        cout << "  elementBufferID = " << m_elementBufferID << endl;
}

void Object::SetShader(const GLuint programID)
{
    m_programID = programID;
    if(DISPLAY_DEBUG_INFO)
        cout << "  programID = " << m_programID << endl;

    UpdateAttributeLocations();
}

void Object::UpdateAttributeLocations()
{
    m_positionID = glGetAttribLocation(m_programID, "in_position");
    if(DISPLAY_DEBUG_INFO)
        cout << "positionID = " << m_positionID << endl;

    m_normalID = glGetAttribLocation(m_programID, "in_normal");
    if(DISPLAY_DEBUG_INFO)
        cout << "normalID = " << m_normalID << endl;

    m_colorID = glGetAttribLocation(m_programID, "in_color");
    if(DISPLAY_DEBUG_INFO)
        cout << "colorID = " << m_colorID << endl;
}


void Object::Draw(const mat4& projection_matrix, const mat4& view_matrix, const GLuint PmatrixID, const GLuint VmatrixID, bool draw_points) const
{
    // Shader program setting
    glUseProgram(m_programID);


    // Matrix transmission
    glUniformMatrix4fv(PmatrixID, 1, GL_FALSE, value_ptr(projection_matrix));
    glUniformMatrix4fv(VmatrixID, 1, GL_FALSE, value_ptr(view_matrix));


    // Pointer settings
    glEnableVertexAttribArray(m_positionID);
    glBindBuffer(GL_ARRAY_BUFFER, m_vertexBufferID);
    glVertexAttribPointer(
                m_positionID,
                3,
                GL_FLOAT,
                GL_FALSE,
                0,
                (void*)0
                );


    glEnableVertexAttribArray(m_normalID);
    glBindBuffer(GL_ARRAY_BUFFER, m_normalBufferID);
    glVertexAttribPointer(
                m_normalID,
                3,
                GL_FLOAT,
                GL_TRUE,
                0,
                (void*)0
                );


    glEnableVertexAttribArray(m_colorID);
    glBindBuffer(GL_ARRAY_BUFFER, m_colorBufferID);
    glVertexAttribPointer(
                m_colorID,
                3,
                GL_FLOAT,
                GL_FALSE,
                0,
                (void*)0
                );


    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_elementBufferID);

    if (draw_points) {
        // Draw points
        glDrawElements(
                    GL_POINTS,               // mode
                    m_mesh->NbFaces()*3,        // count
                    GL_UNSIGNED_INT,            // type
                    (void*)0                    // offset
                    );
    } else {
        // Draw triangles
        glDrawElements(
                    GL_TRIANGLES,               // mode
                    m_mesh->NbFaces()*3,        // count
                    GL_UNSIGNED_INT,            // type
                    (void*)0                    // offset
                    );
    }


    glDisableVertexAttribArray(m_positionID);
    glDisableVertexAttribArray(m_normalID);
    glDisableVertexAttribArray(m_colorID);
}
