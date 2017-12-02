#include <MeshHE.h>
#include <Mesh.h>

#include <iostream>

using namespace glm;
using namespace std;


void AbstractMesh::ClearRessources()
{
    m_positions.clear();
    m_normals.clear();
    m_colors.clear();
    m_indices.clear();
}


//***************
// Vertices color manipulation

void AbstractMesh::ColorFill(const glm::vec3 color)
{
    m_colors.resize(m_positions.size(), color);
}


void AbstractMesh::ColorFromNormals()
{
    m_colors.clear();
    m_colors.reserve(m_positions.size());
    for(glm::uint i=0; i<m_positions.size(); i++)
    {
        m_colors.push_back(vec3(0.5) + 0.5f*m_normals[i]);
    }
}

void AbstractMesh::ColorFromPositions()
{
    m_colors.clear();
    m_colors.reserve(m_positions.size());

    vector<vec3> bb = ComputeBB();

    for(glm::uint i=0; i<m_positions.size(); i++)
    {
        m_colors.push_back((m_positions[i] - bb[0]) / (bb[1] - bb[0]));
    }
}


//***************
// I/O

void AbstractMesh::write_obj(const char* filename) const
{
    FILE *file;

    if((file=fopen(filename,"w"))==NULL)
    {
        std::cout << "Unable to open : " << filename << std::endl;
    }

    fprintf(file, "OFF\n%d %d 0\n", m_positions.size(), m_indices.size()/3);

    for(unsigned int i = 0; i < m_positions.size(); i++)
    {
        vec3 p = m_positions[i];
        fprintf(file,"%f %f %f\n", p[0], p[1], p[2]);
    }

    for(unsigned int i = 0; i < m_indices.size(); i+=3)
    {   
        fprintf(file, "3 %i %i %i\n", m_indices[i], m_indices[i+1], m_indices[i+2]);
    }
}



//***************
// OpenGL utilities

const vector<vec3>& AbstractMesh::positions_array() const
{
    return m_positions;
}


const vector<vec3>& AbstractMesh::normals_array() const
{
    return m_normals;
}


const vector<vec3>& AbstractMesh::colors_array() const
{
    return m_colors;
}


const vector<glm::uint>& AbstractMesh::faces_array() const
{
    return m_indices;
}




//***************
// Geometric utilities


vector<vec3> AbstractMesh::ComputeBB() const
{
    vector<vec3> output;

    output.push_back(m_positions[0]);
    output.push_back(m_positions[0]);

    for(int i=1; i<m_positions.size(); ++i)
    {
        vec3 v = m_positions[i];

        output[0] = glm::min(output[0], v);
        output[1] = glm::max(output[1], v);

    }

    return output;
}


void AbstractMesh::Normalize()
{
    vector<vec3> bb = ComputeBB();

    vec3 centre = (bb[0] + bb[1])*0.5f;
    float radius = glm::max(glm::max(bb[1].x - bb[0].x, bb[1].y - bb[0].y), bb[1].z - bb[0].z);

    for(int i=0; i<m_positions.size(); ++i)
    {
        m_positions[i] = (m_positions[i] - centre) / radius;
    }
}


void AbstractMesh::ComputeNormals()
{
    m_normals.resize(m_positions.size());

    for(unsigned int i=0; i<m_normals.size(); i++)
    {
        m_normals[i] = vec3(0.0);
    }

    for(unsigned int i=0; i<m_indices.size(); i+=3)
    {
        unsigned int i0 = m_indices[i+0];
        unsigned int i1 = m_indices[i+1];
        unsigned int i2 = m_indices[i+2];

        vec3 p0 = m_positions[i0];
        vec3 p1 = m_positions[i1];
        vec3 p2 = m_positions[i2];

        vec3 d01 = glm::normalize(p1 - p0);
        vec3 d12 = glm::normalize(p2 - p1);
        vec3 d20 = glm::normalize(p0 - p2);

        vec3 faceNormal = glm::normalize(glm::cross(d01, -d20));

        float alpha0 = asin(length(glm::cross(d01, -d20)));
        float alpha1 = asin(length(glm::cross(d12, -d01)));
        float alpha2 = asin(length(glm::cross(d20, -d12)));

        if(glm::isnan(alpha0))
            alpha0 = 1.0f;
        if(glm::isnan(alpha1))
            alpha1 = 1.0f;
        if(glm::isnan(alpha2))
            alpha2 = 1.0f;

        m_normals[i0] += faceNormal * alpha0;
        m_normals[i1] += faceNormal * alpha1;
        m_normals[i2] += faceNormal * alpha2;
    }

    for(unsigned int i=0; i<m_normals.size(); i++)
    {
        float l = length(m_normals[i]);

        if(l < 1e-5)
        {
            cerr << "Error : normal of length 0. Normal of indice " << i << endl;
            continue;
        }

        m_normals[i] /= l;
    }
}

