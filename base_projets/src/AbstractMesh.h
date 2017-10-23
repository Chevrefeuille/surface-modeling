#ifndef ABSTRACT_MESH_HE_H
#define ABSTRACT_MESH_HE_H

#include <glm/glm.hpp>

#include <vector>

/**
 * @brief The MeshHE class.
 * Implements the half edge data structure for triangular meshes.
 */
class AbstractMesh
{
 public:

    AbstractMesh(){}                                    /// Standard constructor

    ~AbstractMesh(){ClearRessources();}

    virtual void ClearRessources();                             /// Simple ressources de-allocation

    // Vertices color manipulation
    void ColorFill(const glm::vec3 color);
    void ColorFromNormals();
    void ColorFromPositions();

    std::vector< glm::vec3 > computeBB() const ;        /// Computes the bounding box of this mesh (usefull for normalization)

    // OpenGL utilities
    virtual const std::vector<glm::vec3> &positions_array() const;      /// Generates a contiguous representation for vertices positions
    virtual const std::vector<glm::vec3> &normals_array() const;        /// Generates a contiguous representation for vertices normals
    virtual const std::vector<glm::vec3> &colors_array() const;         /// Generates a contiguous representation for vertices colors
    virtual const std::vector<glm::uint> &faces_array() const;


    // Sizes
    virtual glm::uint NbVertices() const { return m_positions.size();}
    virtual glm::uint NbFaces()    const { return m_indices.size()/3;}


    // Accessors
    const std::vector< glm::vec3 >& Vertices() const {return m_positions;}
    const std::vector< glm::vec3 >& Normals()  const {return m_normals;}
    const std::vector< glm::vec3 >& Colors()  const {return m_colors;}
    const std::vector< glm::uint >& Faces()    const {return m_indices;}
    std::vector< glm::vec3 >& Vertices() {return m_positions;}
    std::vector< glm::vec3 >& Normals()  {return m_normals;}
    std::vector< glm::vec3 >& Colors()  {return m_colors;}
    std::vector< glm::uint >& Faces()    {return m_indices;}


    // I/O
    virtual void display() const {}                 /// Displays some info about this mesh in the console
    void write_obj(const char* filename) const;


    // Geometric utilities
    std::vector< glm::vec3 > ComputeBB() const ;    /// Computes the bounding box of this mesh (usefull for normalization)
    void Normalize();                               /// Normalises and centers this mesh
    virtual void ComputeNormals();                  /// Computes new normals for each vertex of this mesh


protected:

    // Attributes
    std::vector<glm::vec3> m_positions;             /// Container for the vertices positions
    std::vector<glm::vec3> m_normals;               /// Container for the vertices normals
    std::vector<glm::vec3> m_colors;                /// Container for the vertices colors
    std::vector<glm::uint> m_indices;               /// Container for the faces indices
};

#endif // ABSTRACT_MESH_HE_H

