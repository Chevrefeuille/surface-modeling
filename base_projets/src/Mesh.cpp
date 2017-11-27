#include <cstdlib>
#include <iostream>
#include <math.h>
#include "omp.h"
#include "MeshHE.h"

#include <Mesh.h>
#include <ImplicitFunction.h>
#include <DataSet.h>

using namespace glm;
using namespace std;

/**
 * Réalise nb_itérations itérations de dichotomie, afin de trouver l'isovalue de function sur [p0, p1]
 */
vec3 findRoot(const ImplicitFunction& function, const float isoValue, const vec3& p0, const vec3& p1, unsigned nb_iter = 10)
{
    vec3 p00 = p0;
    vec3 p10 = p1;
    if(function.Eval(p0) > function.Eval(p1))
    {
        swap(p00, p10);
    }

    vec3 p = 0.5f*(p00+p10);
    for(unsigned int iter = 0; iter < nb_iter; iter++)
    {
        if(function.Eval(p) > isoValue)
        {
            p10 = 0.5f * (p00 + p10);
        }
        else
        {
            p00 = 0.5f * (p00 + p10);
        }

        p = 0.5f*(p00+p10);
    }

    return p;
}

/**
 * Calcule les points approchant l'isovalue de function sur le tétrahèdre défini par les 4 points de p[]
 */
void Mesh::ProcessTetrahedron2(const ImplicitFunction& function, const vec3 p[]) {

	// Vecteur donnant les signes des images par f de chacun des points du tétrahèdre
	bool b[4] = {function.Eval(p[0]) > 0.0, function.Eval(p[1]) > 0.0, function.Eval(p[2]) > 0.0, function.Eval(p[3]) > 0.0};

	unsigned int N = m_positions.size();

	// ** Distinction des cas selon le signe des images par f des 4 points du tétrahèdre : ** //

	// 1. 4 positifs (resp 4 négatifs)
	if((!b[0] && !b[1] && !b[2] && !b[3]) || (b[0] && b[1] && b[2] && b[3])) { return; }

	// 2. 3 positifs et 1 négatif (resp 3 négatifs et 1 positif)
	for(unsigned int i=0; i<4; i++) {
		if((b[i] && !b[(i+1)%4] && !b[(i+2)%4] && !b[(i+3)%4]) || (!b[i] && b[(i+1)%4] && b[(i+2)%4] && b[(i+3)%4])) {
			// Le point i est celui au signe différent des 3 autres
			// On forme un triangle
			vec3 p0 = findRoot(function, 0.0, p[i], p[(i+1)%4]);
			vec3 p1 = findRoot(function, 0.0, p[i], p[(i+2)%4]);
			vec3 p2 = findRoot(function, 0.0, p[i], p[(i+3)%4]);
            m_positions.push_back(p0);
            m_positions.push_back(p1);
            m_positions.push_back(p2);

			// Choix orientation des normales : -
			vec3 n0 = glm::normalize(-function.EvalDev(p0));
			vec3 n1 = glm::normalize(-function.EvalDev(p1));
			vec3 n2 = glm::normalize(-function.EvalDev(p2));
			m_normals.push_back(n0);
            m_normals.push_back(n1);
            m_normals.push_back(n2);

			// Indices dans ordre trigo
			if(dot(cross(p1-p0, p2-p0), n0+n1+n2)>0) {
				m_indices.push_back(N);
                m_indices.push_back(N+1);
                m_indices.push_back(N+2);
			}
			else {
                m_indices.push_back(N);
                m_indices.push_back(N+2);
                m_indices.push_back(N+1);
			}
			return;
		}
	}

	// 3. 2 positifs et 2 négatifs
	// Sélection des différentes combinaisons de points : on cherche pi et pj tq leur image soit de même signe
	for(unsigned int i=0; i<3; i++) {
		for(unsigned int j=i+1; j<4; j++) {
			unsigned int k = (i+1)%4;
			if(k == j) {
				k += 1;
				k %= 4;
			}

			unsigned int l = (k+1)%4;
			if(l == i) {
				l += 1;
				l %= 4;
			}
			if(l == j) {
				l += 1;
				l %= 4;
			}

			if(i == k || j == k || i == l || j == l || k == l)
				cerr << "indices are wrong ! " << endl;

			// Si les points pi et pj ont leur image de même signe
			// On forme un quadrilatère
			if(b[i] && b[j] && !b[k] && !b[l] || !b[i] && !b[j] && b[k] && b[l]) {
				vec3 p0 = findRoot(function, 0.0, p[i], p[k]);
				vec3 p1 = findRoot(function, 0.0, p[i], p[l]);
				vec3 p2 = findRoot(function, 0.0, p[j], p[k]);
				vec3 p3 = findRoot(function, 0.0, p[j], p[l]);
                m_positions.push_back(p0);
                m_positions.push_back(p1);
                m_positions.push_back(p2);
                m_positions.push_back(p3);

				vec3 n0 = glm::normalize(-function.EvalDev(p0));
				vec3 n1 = glm::normalize(-function.EvalDev(p1));
				vec3 n2 = glm::normalize(-function.EvalDev(p2));
				vec3 n3 = glm::normalize(-function.EvalDev(p3));
                m_normals.push_back(n0);
                m_normals.push_back(n1);
                m_normals.push_back(n2);
                m_normals.push_back(n3);

				if(dot(cross(p2-p0, p3-p0), n0+n3+n2)>0) {
                    m_indices.push_back(N);
                    m_indices.push_back(N+2);
                    m_indices.push_back(N+3);

                    m_indices.push_back(N);
                    m_indices.push_back(N+3);
                    m_indices.push_back(N+1);
				} else {
                    m_indices.push_back(N);
                    m_indices.push_back(N+3);
                    m_indices.push_back(N+2);

                    m_indices.push_back(N);
                    m_indices.push_back(N+1);
                    m_indices.push_back(N+3);
				}

				return;
			}
		}
	}
	cerr << "no solution found in marching tetrahedron !!" << endl;

}

/**
 * Construit un mesh défini par son iso fonction, les coord limites du marching cubes et sa résolution
 */
Mesh::Mesh(const ImplicitFunction& function, double minX, double maxX, double minY, double maxY,
		double minZ, double maxZ, const unsigned int resX, const unsigned int resY,
		const unsigned int resZ) {
	m_positions = vector<vec3>();
	m_normals  = vector<vec3>();
	m_indices    = vector<unsigned int>();

	// Boucles correspondant au marching cubes (découpage en cubes)
	for(unsigned int i=0; i < resX; i++) {
        double x0 = double(i  )/resX * (maxX - minX) + minX;
        double x1 = double(i+1)/resX * (maxX - minX) + minX;
		for(unsigned int j=0; j < resY; j++) {
            double y0 = double(j  )/resY * (maxY - minY) + minY;
            double y1 = double(j+1)/resY * (maxY - minY) + minY;
			for(unsigned int k=0; k < resZ; k++) {
                double z0 = double(k  )/resZ * (maxZ - minZ) + minZ;
                double z1 = double(k+1)/resZ * (maxZ - minZ) + minZ;

				// Les 8 points du cube défini par (i,j,k)
				vec3 p000 = vec3(x0, y0, z0);
				vec3 p001 = vec3(x1, y0, z0);
				vec3 p010 = vec3(x0, y1, z0);
				vec3 p011 = vec3(x1, y1, z0);
				vec3 p100 = vec3(x0, y0, z1);
				vec3 p101 = vec3(x1, y0, z1);
				vec3 p110 = vec3(x0, y1, z1);
				vec3 p111 = vec3(x1, y1, z1);

				// Découpage du cube en 6 tétrahèdres (définis chacun par 4 points)
				vec3 p0[4] =  {p000, p001, p011, p111};
				vec3 p1[4] =  {p000, p011, p010, p111};
				vec3 p2[4] =  {p000, p010, p110, p111};
				vec3 p3[4] =  {p000, p110, p100, p111};
				vec3 p4[4] =  {p000, p100, p101, p111};
				vec3 p5[4] =  {p000, p101, p001, p111};

				// Marching tetrahedra
				this->ProcessTetrahedron2(function, p0);
                this->ProcessTetrahedron2(function, p1);
                this->ProcessTetrahedron2(function, p2);
                this->ProcessTetrahedron2(function, p3);
                this->ProcessTetrahedron2(function, p4);
                this->ProcessTetrahedron2(function, p5);
			}
		}
	}
}

double inscribedCercleRadiusOfTriangle(const vec3& p1, const vec3& p2, const vec3& p3) {
    double areaTriangle = length(cross(p2-p1, p3-p1))/2;
    double halfPerimeter = length(p2-p1) + length(p3-p2) + length(p3-p1);
    halfPerimeter/=2;
    return areaTriangle/halfPerimeter;
}

/**
 * Edge collapsing
 * On fait pointer les anciens index sur le nouvel index
 */
void Mesh::collapseEdge(unsigned int oldIndex1, unsigned int oldIndex2, unsigned int newIndex) {
    unsigned int n = m_indices.size() / 3;
    for (int i = 0; i < n; i++) {
        // Si le triangle contient les 2 sommets à unifier:
        if (m_indices[i] == oldIndex1 || m_indices[i+1] == oldIndex1 || m_indices[i+2] == oldIndex1) {
            if (m_indices[i] == oldIndex2 || m_indices[i+1] == oldIndex2 || m_indices[i+2] == oldIndex2) {
                m_indices.erase(m_indices.begin()+i, m_indices.begin()+i+3);
                n-=3;
            }
        } else {
            if (m_indices[i] == oldIndex1 || m_indices[i] == oldIndex2) {
                m_indices[i] = newIndex;
            }
            if (m_indices[i + 1] == oldIndex1 || m_indices[i + 1] == oldIndex2) {
                m_indices[i + 1] = newIndex;
            }
            if (m_indices[i + 2] == oldIndex1 || m_indices[i + 2] == oldIndex2) {
                m_indices[i + 2] = newIndex;
            }
        }

    }
}

Mesh Mesh::postProcess(const double epsilon) {
    // Supprime les points en double et refait la correspondance avec les indices
    //this->RemoveDouble();
    //MeshHE m_he(*this);
    //return m_he;

    unsigned int n = m_indices.size()/3;

    double inscribedCercleRadius;
    double edgeAspectRatio[3*n];
    unsigned int orderedIndexEdgeAspectRatio[n];
    double* orderedEdgeAspectRatio2 = new double[n];
    vec3 p1, p2, p3;

/*
    for (int i = 0; i < 3*n; i++) {
        // 0 : p2 - p1
        // 1 : p3 - p2
        // 2 : p3 - p1
        p1 = m_positions[m_indices[3*i+0]];
        p2 = m_positions[m_indices[3*i+1]];
        p3 = m_positions[m_indices[3*i+2]];
        inscribedCercleRadius = inscribedCercleRadiusOfTriangle(p1, p2, p3);;
        edgeAspectRatio[3*i+0] = inscribedCercleRadius * length(p2-p1);
        edgeAspectRatio[3*i+1] = inscribedCercleRadius * length(p3-p2);
        edgeAspectRatio[3*i+2] = inscribedCercleRadius * length(p3-p1);
    }

    // Indicage edgeAspectRatio par ordre croissant
    unsigned int minIndex = 0; double minRatio = edgeAspectRatio[0];
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 1; j < 3*n; j++) {
            if (edgeAspectRatio[j] < minRatio) {
                minRatio = edgeAspectRatio[j];
                minIndex = j;
            }
        }
        orderedIndexEdgeAspectRatio[i] = minIndex;
        orderedEdgeAspectRatio[i] = minRatio;
        edgeAspectRatio[minIndex] = 1;
        minIndex = 0; minRatio = edgeAspectRatio[0];
    }

    for (int i = 0; i<n; i++)
        printf("ratio = %lf\n", orderedEdgeAspectRatio[i]);

     //Elimination des arretes aux plus petit ratio
     double ratio;
     unsigned int indexEdge, ip1, ip2;
     unsigned int N;
     for (int i = 0; i < n; i++) {
         indexEdge = orderedIndexEdgeAspectRatio[i];
         ratio = orderedEdgeAspectRatio[i];
         if (ratio >= epsilon) {break;}

         printf("Collapsing edge of ratio %lf", ratio);

         // indices points à "unifier" :
         if (indexEdge % 3 == 0) {
             ip1 = m_indices[indexEdge];
             ip2 = m_indices[indexEdge+1];
         } else if (indexEdge % 3 == 1) {
             ip1 = m_indices[indexEdge];
             ip2 = m_indices[indexEdge+1];
         } else {
             ip1 = m_indices[indexEdge];
             ip2 = m_indices[indexEdge-2];
         }

         // Nouveau point = milieu des 2 anciens points
         p3 = (m_positions[ip1] + m_positions[ip2]) / (float)2;
         N = m_positions.size(); // indice du nouveau point p3
         m_indices.push_back(N);
         m_positions.push_back(p3);
         // On fait pointer les anciens sommets ip1 et ip2 vers N
         this->collapseEdge(ip1, ip2, N);
     }
*/
	 return *this;
}


Mesh::Mesh(const char* filename)
{
	int j = 0;
    unsigned int tmp;
	float *n;
	FILE *file;
	int   error;
	float r;

	if((file=fopen(filename,"r"))==NULL)
	{
		std::cout << "Unable to read : " << filename << std::endl;
	}

	// create mesh
    m_positions = vector<vec3>();
    m_normals  = vector<vec3>();
    m_indices    = vector<unsigned int>();

    glm::uint nb_vertices, nb_faces;

	error = fscanf(file,"OFF\n%d %d %d\n",&(nb_vertices),&(nb_faces),&tmp);
	if(error==EOF)
	{
		std::cout << "Unable to read : " << filename << std::endl;
	}

    m_positions.resize(nb_vertices);
    m_normals.resize(nb_vertices);
    m_indices.resize(nb_faces*3);

	// reading vertices
	for(int i=0;i<nb_vertices;++i)
	{
        error = fscanf(file,"%f %f %f\n",&(m_positions[i][0]),&(m_positions[i][1]),&(m_positions[i][2]));
		if(error==EOF)
		{
			std::cout << "Unable to read vertices of : " << filename << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	// reading faces
	j = 0;
	for(int i=0;i<nb_faces;++i)
	{
        error = fscanf(file,"%d %d %d %d\n",&tmp,&(m_indices[j]),&(m_indices[j+1]),&(m_indices[j+2]));

		if(error==EOF)
		{
			std::cout << "Unable to read faces of : " << filename << std::endl;
			exit(EXIT_FAILURE);
		}

		if(tmp!=3)
		{
			printf("Error : face %d is not a triangle (%d polygonal face!)\n",i/3,tmp);
			exit(EXIT_FAILURE);
		}
		j += 3;
	}

    fclose(file);

    ColorFill(vec3(0.9));
}



void Mesh::CreateCube(Mesh& mesh)
{
    mesh.m_positions.push_back(vec3(-1, -1, -1));
    mesh.m_positions.push_back(vec3( 1, -1, -1));
    mesh.m_positions.push_back(vec3( 1,  1, -1));
    mesh.m_positions.push_back(vec3(-1,  1, -1));

    mesh.m_positions.push_back(vec3(-1, -1, -1));
    mesh.m_positions.push_back(vec3(-1,  1, -1));
    mesh.m_positions.push_back(vec3(-1,  1,  1));
    mesh.m_positions.push_back(vec3(-1, -1,  1));

    mesh.m_positions.push_back(vec3(-1, -1, -1));
    mesh.m_positions.push_back(vec3(-1, -1,  1));
    mesh.m_positions.push_back(vec3( 1, -1,  1));
    mesh.m_positions.push_back(vec3( 1, -1, -1));

    mesh.m_positions.push_back(vec3( 1,  1,  1));
    mesh.m_positions.push_back(vec3(-1,  1,  1));
    mesh.m_positions.push_back(vec3(-1, -1,  1));
    mesh.m_positions.push_back(vec3( 1, -1,  1));

    mesh.m_positions.push_back(vec3( 1,  1,  1));
    mesh.m_positions.push_back(vec3( 1, -1,  1));
    mesh.m_positions.push_back(vec3( 1, -1, -1));
    mesh.m_positions.push_back(vec3( 1,  1, -1));

    mesh.m_positions.push_back(vec3( 1,  1,  1));
    mesh.m_positions.push_back(vec3( 1,  1, -1));
    mesh.m_positions.push_back(vec3(-1,  1, -1));
    mesh.m_positions.push_back(vec3(-1,  1,  1));

    mesh.m_normals.push_back(vec3(0, 0, -1));
    mesh.m_normals.push_back(vec3(0, 0, -1));
    mesh.m_normals.push_back(vec3(0, 0, -1));
    mesh.m_normals.push_back(vec3(0, 0, -1));

    mesh.m_normals.push_back(vec3(-1, 0, 0));
    mesh.m_normals.push_back(vec3(-1, 0, 0));
    mesh.m_normals.push_back(vec3(-1, 0, 0));
    mesh.m_normals.push_back(vec3(-1, 0, 0));

    mesh.m_normals.push_back(vec3(0, -1, 0));
    mesh.m_normals.push_back(vec3(0, -1, 0));
    mesh.m_normals.push_back(vec3(0, -1, 0));
    mesh.m_normals.push_back(vec3(0, -1, 0));

    mesh.m_normals.push_back(vec3(0, 0,  1));
    mesh.m_normals.push_back(vec3(0, 0,  1));
    mesh.m_normals.push_back(vec3(0, 0,  1));
    mesh.m_normals.push_back(vec3(0, 0,  1));

    mesh.m_normals.push_back(vec3( 1, 0, 0));
    mesh.m_normals.push_back(vec3( 1, 0, 0));
    mesh.m_normals.push_back(vec3( 1, 0, 0));
    mesh.m_normals.push_back(vec3( 1, 0, 0));

    mesh.m_normals.push_back(vec3(0,  1, 0));
    mesh.m_normals.push_back(vec3(0,  1, 0));
    mesh.m_normals.push_back(vec3(0,  1, 0));
    mesh.m_normals.push_back(vec3(0,  1, 0));


    mesh.m_indices.push_back(0);
    mesh.m_indices.push_back(2);
    mesh.m_indices.push_back(1);

    mesh.m_indices.push_back(0);
    mesh.m_indices.push_back(3);
    mesh.m_indices.push_back(2);


    mesh.m_indices.push_back(4);
    mesh.m_indices.push_back(6);
    mesh.m_indices.push_back(5);

    mesh.m_indices.push_back(4);
    mesh.m_indices.push_back(7);
    mesh.m_indices.push_back(6);


    mesh.m_indices.push_back(8);
    mesh.m_indices.push_back(10);
    mesh.m_indices.push_back(9);

    mesh.m_indices.push_back(8);
    mesh.m_indices.push_back(11);
    mesh.m_indices.push_back(10);


    mesh.m_indices.push_back(12);
    mesh.m_indices.push_back(13);
    mesh.m_indices.push_back(14);

    mesh.m_indices.push_back(12);
    mesh.m_indices.push_back(14);
    mesh.m_indices.push_back(15);


    mesh.m_indices.push_back(16);
    mesh.m_indices.push_back(17);
    mesh.m_indices.push_back(18);

    mesh.m_indices.push_back(16);
    mesh.m_indices.push_back(18);
    mesh.m_indices.push_back(19);


    mesh.m_indices.push_back(20);
    mesh.m_indices.push_back(21);
    mesh.m_indices.push_back(22);

    mesh.m_indices.push_back(20);
    mesh.m_indices.push_back(22);
    mesh.m_indices.push_back(23);
}

void Mesh::CreateCube2(Mesh& mesh)
{
    for(float x = -1.0; x < 3.0; x += 2.0)
        for(float y = -1.0; y < 3.0; y += 2.0)
            for(float z = -1.0; z < 3.0; z += 2.0)
            {
                mesh.m_positions.push_back(vec3(x, y, z));
                mesh.m_normals.push_back(glm::normalize(vec3(x, y, z)));
            }

    mesh.m_indices.push_back(0);
    mesh.m_indices.push_back(1);
    mesh.m_indices.push_back(3);

    mesh.m_indices.push_back(0);
    mesh.m_indices.push_back(3);
    mesh.m_indices.push_back(2);


    mesh.m_indices.push_back(0);
    mesh.m_indices.push_back(2);
    mesh.m_indices.push_back(6);

    mesh.m_indices.push_back(0);
    mesh.m_indices.push_back(6);
    mesh.m_indices.push_back(4);


    mesh.m_indices.push_back(0);
    mesh.m_indices.push_back(4);
    mesh.m_indices.push_back(5);

    mesh.m_indices.push_back(0);
    mesh.m_indices.push_back(5);
    mesh.m_indices.push_back(1);




    mesh.m_indices.push_back(1);
    mesh.m_indices.push_back(7);
    mesh.m_indices.push_back(3);

    mesh.m_indices.push_back(3);
    mesh.m_indices.push_back(7);
    mesh.m_indices.push_back(2);


    mesh.m_indices.push_back(2);
    mesh.m_indices.push_back(7);
    mesh.m_indices.push_back(6);

    mesh.m_indices.push_back(6);
    mesh.m_indices.push_back(7);
    mesh.m_indices.push_back(4);


    mesh.m_indices.push_back(4);
    mesh.m_indices.push_back(7);
    mesh.m_indices.push_back(5);

    mesh.m_indices.push_back(5);
    mesh.m_indices.push_back(7);
    mesh.m_indices.push_back(1);
}



void Mesh::CreateSphere(Mesh& mesh, unsigned int Nu, unsigned int Nv)
{

    for(int i = 0; i < Nu; i++)
    {
        float u = float(i) / (Nu-1);
        float phi = u * M_PI * 2;

        for(int j = 0; j < Nv; j++)
        {
            float v = float(j) / (Nv-1);
            float psi = v * M_PI;

            vec3 p(cos(phi)*sin(psi), sin(phi)*sin(psi), cos(psi));
            mesh.m_positions.push_back(p);

            mesh.m_normals.push_back(p);

        }
    }

    for(int i = 0; i < Nu; i++)
    {
        for(int j = 0; j < Nv - 1; j++)
        {
            mesh.m_indices.push_back( i * Nv + j);
            mesh.m_indices.push_back((i+1)%Nu * Nv + j);
            mesh.m_indices.push_back( i * Nv + j+1);

            mesh.m_indices.push_back( i * Nv + j+1);
            mesh.m_indices.push_back((i+1)%Nu * Nv + j);
            mesh.m_indices.push_back((i+1)%Nu * Nv + j + 1);
        }
    }
}


void Mesh::RemoveDouble(float epsilon)
{
    vector<unsigned int> dbl;
    for(unsigned int i=0; i<m_positions.size(); i++)
        dbl.push_back(i);

    for(unsigned int i=0; i<m_positions.size()-1; i++)
    {
        if(dbl[i] != i)
            continue;

        for(unsigned int j=i+1; j<m_positions.size(); j++)
        {
            if(length(m_positions[i] - m_positions[j]) < epsilon)
            {
                dbl[j] = i;
            }
        }
    }


    for(unsigned int i=0; i<m_positions.size(); i++)
    {
        while(dbl[dbl[i]] != dbl[i])
        {
            dbl[i] = dbl[dbl[i]];
        }
    }

    vector<vec3> new_vertices;
    vector<int> corresp;

    for(unsigned int i=0; i<m_positions.size(); i++)
    {
        if(dbl[i] == i)
        {
            corresp.push_back(new_vertices.size());
            new_vertices.push_back(m_positions[i]);
            continue;
        }

        corresp.push_back(-1);
    }


    for(unsigned int i=0; i<m_positions.size(); i++)
    {
        while(corresp[i] == -1)
        {
            corresp[i] = corresp[dbl[i]];
        }
    }



    for(unsigned int i=0; i<m_indices.size(); i++)
    {
        m_indices[i] = corresp[m_indices[i]];
    }

    m_positions = new_vertices;

}



void Mesh::CreateIsoSurface(  Mesh& mesh
                            , const ImplicitFunction& function
                            , const float isoValue
                            , const float minX
                            , const float maxX
                            , const float minY
                            , const float maxY
                            , const float minZ
                            , const float maxZ
                            , const unsigned int resX
                            , const unsigned int resY
                            , const unsigned int resZ)
{
    for(unsigned int i=0; i < resX; i++)
    {
        float x0 = float(i  )/resX * (maxX - minX) + minX;
        float x1 = float(i+1)/resX * (maxX - minX) + minX;
        for(unsigned int j=0; j < resY; j++)
        {
            float y0 = float(j  )/resY * (maxY - minY) + minY;
            float y1 = float(j+1)/resY * (maxY - minY) + minY;
            for(unsigned int k=0; k < resZ; k++)
            {
                float z0 = float(k  )/resZ * (maxZ - minZ) + minZ;
                float z1 = float(k+1)/resZ * (maxZ - minZ) + minZ;

                vec3 p000 = vec3(x0, y0, z0);
                vec3 p001 = vec3(x1, y0, z0);
                vec3 p010 = vec3(x0, y1, z0);
                vec3 p011 = vec3(x1, y1, z0);
                vec3 p100 = vec3(x0, y0, z1);
                vec3 p101 = vec3(x1, y0, z1);
                vec3 p110 = vec3(x0, y1, z1);
                vec3 p111 = vec3(x1, y1, z1);

                vec3 p0[4] =  {p000, p001, p011, p111};
                vec3 p1[4] =  {p000, p011, p010, p111};
                vec3 p2[4] =  {p000, p010, p110, p111};
                vec3 p3[4] =  {p000, p110, p100, p111};
                vec3 p4[4] =  {p000, p100, p101, p111};
                vec3 p5[4] =  {p000, p101, p001, p111};

                ProcessTetrahedron(mesh, function, isoValue,  p0);
                ProcessTetrahedron(mesh, function, isoValue,  p1);
                ProcessTetrahedron(mesh, function, isoValue,  p2);
                ProcessTetrahedron(mesh, function, isoValue,  p3);
                ProcessTetrahedron(mesh, function, isoValue,  p4);
                ProcessTetrahedron(mesh, function, isoValue,  p5);
            }
        }
    }

}

void Mesh::ProcessTetrahedron(Mesh& mesh, const ImplicitFunction& function, const float isoValue, const vec3 p[])
{
    bool b[4] = {function.Eval(p[0]) > isoValue, function.Eval(p[1]) > isoValue, function.Eval(p[2]) > isoValue, function.Eval(p[3]) > isoValue};


    unsigned int N = mesh.m_positions.size();
    if(!b[0] && !b[1] && !b[2] && !b[3] || b[0] && b[1] && b[2] && b[3])
    {
        return;
    }

    for(unsigned int i=0; i<4; i++)
    {
        if(b[i] && !b[(i+1)%4] && !b[(i+2)%4] && !b[(i+3)%4] || !b[i] && b[(i+1)%4] && b[(i+2)%4] && b[(i+3)%4])
        {

            vec3 p0 = findRoot(function, isoValue, p[i], p[(i+1)%4]);
            vec3 p1 = findRoot(function, isoValue, p[i], p[(i+2)%4]);
            vec3 p2 = findRoot(function, isoValue, p[i], p[(i+3)%4]);


            vec3 n0 = glm::normalize(-function.EvalDev(p0));
            vec3 n1 = glm::normalize(-function.EvalDev(p1));
            vec3 n2 = glm::normalize(-function.EvalDev(p2));


            mesh.m_positions.push_back(p0);
            mesh.m_positions.push_back(p1);
            mesh.m_positions.push_back(p2);

            mesh.m_normals.push_back(n0);
            mesh.m_normals.push_back(n1);
            mesh.m_normals.push_back(n2);


            if(dot(cross(p1-p0, p2-p0), n0+n1+n2)>0)
            {
                mesh.m_indices.push_back(N);
                mesh.m_indices.push_back(N+1);
                mesh.m_indices.push_back(N+2);
            }
            else
            {
                mesh.m_indices.push_back(N);
                mesh.m_indices.push_back(N+2);
                mesh.m_indices.push_back(N+1);
            }

            return;
        }
    }

    for(unsigned int i=0; i<3; i++)
    {
        for(unsigned int j=i+1; j<4; j++)
        {
            unsigned int k = (i+1)%4;
            if(k == j)
            {
                k += 1;
                k %= 4;
            }

            unsigned int l = (k+1)%4;
            if(l == i)
            {
                l += 1;
                l %= 4;
            }
            if(l == j)
            {
                l += 1;
                l %= 4;
            }

            if(i == k || j == k || i == l || j == l || k == l)
                cerr << "indices are wrong ! " << endl;


            if(b[i] && b[j] && !b[k] && !b[l] || !b[i] && !b[j] && b[k] && b[l])
            {
                vec3 p0 = findRoot(function, isoValue, p[i], p[k]);
                vec3 p1 = findRoot(function, isoValue, p[i], p[l]);
                vec3 p2 = findRoot(function, isoValue, p[j], p[k]);
                vec3 p3 = findRoot(function, isoValue, p[j], p[l]);

                vec3 n0 = glm::normalize(-function.EvalDev(p0));
                vec3 n1 = glm::normalize(-function.EvalDev(p1));
                vec3 n2 = glm::normalize(-function.EvalDev(p2));
                vec3 n3 = glm::normalize(-function.EvalDev(p3));


                mesh.m_positions.push_back(p0);
                mesh.m_positions.push_back(p1);
                mesh.m_positions.push_back(p2);
                mesh.m_positions.push_back(p3);


                mesh.m_normals.push_back(n0);
                mesh.m_normals.push_back(n1);
                mesh.m_normals.push_back(n2);
                mesh.m_normals.push_back(n3);



                if(dot(cross(p2-p0, p3-p0), n0+n3+n2)>0)
                {
                    mesh.m_indices.push_back(N);
                    mesh.m_indices.push_back(N+2);
                    mesh.m_indices.push_back(N+3);

                    mesh.m_indices.push_back(N);
                    mesh.m_indices.push_back(N+3);
                    mesh.m_indices.push_back(N+1);
                }
                else
                {
                    mesh.m_indices.push_back(N);
                    mesh.m_indices.push_back(N+3);
                    mesh.m_indices.push_back(N+2);

                    mesh.m_indices.push_back(N);
                    mesh.m_indices.push_back(N+1);
                    mesh.m_indices.push_back(N+3);
                }


                return;
            }
        }
    }

    cerr << "no solution found in marching tetrahedron !!" << endl;
}
