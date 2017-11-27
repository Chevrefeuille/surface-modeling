#include "ImplicitFunction.h"

#include <iostream>
#include <math.h>

using namespace glm;
using namespace std;



ImplicitFunction::ImplicitFunction()
{

}

ImplicitFunction::~ImplicitFunction()
{

}

float ImplicitFunction::Eval(glm::vec3 /*p*/) const
{
    cout << "ImplicitFunction::Eval : TODO : implementation" << endl;
}

vec3 ImplicitFunction::EvalDev(glm::vec3 p) const
{
    return EvalDevFiniteDiff(p);
}

vec3 ImplicitFunction::EvalDevFiniteDiff(glm::vec3 p) const
{
    float epsilon = 1e-3;
    vec3 dx = vec3(epsilon, 0.0, 0.0);
    vec3 dy = vec3(0.0, epsilon, 0.0);
    vec3 dz = vec3(0.0, 0.0, epsilon);

    vec3 dev = vec3(Eval(p + dx) - Eval(p - dx), Eval(p + dy) - Eval(p - dy), Eval(p + dz) - Eval(p - dz));
    dev /= 2*epsilon;

    return dev;
}




SphereFunction::SphereFunction(const glm::vec3& centre, const float radius):
    m_centre(centre),
    m_radius(radius)
{

}

SphereFunction::~SphereFunction()
{

}

float SphereFunction::Eval(glm::vec3 p) const
{
    //return exp(-dot(p - m_centre, p - m_centre) / m_radius / m_radius);
    return (p[0]-m_centre[0])*(p[0]-m_centre[0])+(p[1]-m_centre[1])*(p[1]-m_centre[1])+(p[2]-m_centre[2])*(p[2]-m_centre[2])-m_radius*m_radius;
}

glm::vec3 SphereFunction::EvalDev(glm::vec3 p) const
{
//    return normalize(p - m_centre);
    return -2.0f * (p - m_centre) / m_radius * Eval(p);

//    return EvalDevFiniteDiff(p);
}




AddFunction::AddFunction(const ImplicitFunction& op0, const ImplicitFunction& op1):
    m_op0(op0),
    m_op1(op1)
{

}

AddFunction::~AddFunction()
{

}

float AddFunction::Eval(glm::vec3 p) const
{
    return m_op0.Eval(p) + m_op1.Eval(p);
}

glm::vec3 AddFunction::EvalDev(glm::vec3 p) const
{
    return m_op0.EvalDev(p) + m_op1.EvalDev(p);
}


MulFunction::MulFunction(const ImplicitFunction& op0, const ImplicitFunction& op1):
    m_op0(op0),
    m_op1(op1)
{

}

MulFunction::~MulFunction()
{

}

float MulFunction::Eval(glm::vec3 p) const
{
    return m_op0.Eval(p) * m_op1.Eval(p);
}

glm::vec3 MulFunction::EvalDev(glm::vec3 p) const
{
    return m_op0.EvalDev(p) * m_op1.Eval(p) + m_op0.Eval(p) * m_op1.EvalDev(p);
}

MinFunction::MinFunction(const ImplicitFunction& op0, const ImplicitFunction& op1):
    m_op0(op0),
    m_op1(op1)
{

}

MinFunction::~MinFunction()
{

}

float MinFunction::Eval(glm::vec3 p) const
{
    return glm::min(m_op0.Eval(p), m_op1.Eval(p));
}

glm::vec3 MinFunction::EvalDev(glm::vec3 p) const
{
    if(m_op0.Eval(p) > m_op1.Eval(p))
        return m_op1.EvalDev(p);

    return m_op0.EvalDev(p);
}


MaxFunction::MaxFunction(const ImplicitFunction& op0, const ImplicitFunction& op1):
    m_op0(op0),
    m_op1(op1)
{

}

MaxFunction::~MaxFunction()
{

}

float MaxFunction::Eval(glm::vec3 p) const
{
    return glm::max(m_op0.Eval(p), m_op1.Eval(p));
}

glm::vec3 MaxFunction::EvalDev(glm::vec3 p) const
{
    if(m_op0.Eval(p) < m_op1.Eval(p))
        return m_op1.EvalDev(p);

    return m_op0.EvalDev(p);
}



BlendFunction::BlendFunction(const ImplicitFunction& op0, const ImplicitFunction& op1, const ImplicitFunction& weight):
    m_op0(op0),
    m_op1(op1),
    m_weight(weight)
{

}

BlendFunction::~BlendFunction()
{

}

float BlendFunction::Eval(glm::vec3 p) const
{
    return m_op0.Eval(p) + m_weight.Eval(p) * (m_op1.Eval(p) - m_op0.Eval(p));
}

glm::vec3 BlendFunction::EvalDev(glm::vec3 p) const
{
    return m_op0.EvalDev(p) + m_weight.Eval(p) * (m_op1.EvalDev(p) - m_op0.EvalDev(p)) + m_weight.EvalDev(p) * (m_op1.Eval(p) - m_op0.Eval(p));
}




ScalProdFunction::ScalProdFunction(const float scalar, const ImplicitFunction& func):
    m_scalar(scalar),
    m_func(func)
{

}

ScalProdFunction::~ScalProdFunction()
{

}

float ScalProdFunction::Eval(glm::vec3 p) const
{
    return m_func.Eval(p) * m_scalar;
}

glm::vec3 ScalProdFunction::EvalDev(glm::vec3 p) const
{
    return m_func.EvalDev(p) * m_scalar;
}



AddFunction operator+(const ImplicitFunction& op0, const ImplicitFunction& op1)
{
    return AddFunction(op0, op1);
}

ScalProdFunction operator-(const ImplicitFunction& func)
{
    return ScalProdFunction(-1.0, func);
}

AddFunction operator-(const ImplicitFunction& op0, const ImplicitFunction& op1)
{
    return AddFunction(op0, -op1);
}

ScalProdFunction operator*(const float scalar, const ImplicitFunction& func)
{
    return ScalProdFunction(scalar, func);
}

ScalProdFunction operator*(const ImplicitFunction& func, const float scalar)
{
    return ScalProdFunction(scalar, func);
}


SegmentFunction::SegmentFunction(const glm::vec3& p0, const float r0, const glm::vec3& p1, const float r1):
    m_p0(p0),
    m_r0(r0),
    m_p1(p1),
    m_r1(r1)
{

}

SegmentFunction::~SegmentFunction()
{

}

float SegmentFunction::Eval(glm::vec3 p) const
{
//    float t = dot(m_p1-m_p0,p-m_p0) / dot(m_p1-m_p0, m_p1-m_p0);

//    if(t<0)
//        return exp(-dot(p - m_p0, p - m_p0) / m_r0 / m_r0);

//    if(t>=1)
//        return exp(-dot(p - m_p1, p - m_p1) / m_r1 / m_r1);

//    vec3 p_proj = m_p0+t*(m_p1-m_p0);
//    float r_proj = m_r0+t*(m_r1-m_r0);
//    return exp(-dot(p - p_proj, p - p_proj) / r_proj / r_proj);

    float t = dot(m_p1-m_p0,p-m_p0) / dot(m_p1-m_p0, m_p1-m_p0);

    vec3 p_proj = m_p0 + t * (m_p1 - m_p0);

    float n = 0; // length(p - p_proj) / length(m_p1-m_p0); //
    float cos_alpha = (m_r0 - m_r1)/length(m_p1-m_p0);

    t -= cos_alpha * n;
    p_proj = m_p0 + t * (m_p1 - m_p0);



    if(t<0)
        return exp(-dot(p - m_p0, p - m_p0) / m_r0 / m_r0);

    if(t>=1)
        return exp(-dot(p - m_p1, p - m_p1) / m_r1 / m_r1);

//    vec3 p_proj = m_p0+t*(m_p1-m_p0);

//    t = cos((t+1.0) * M_PI) * 0.5 + 0.5;
    t = t*t*(3-2*t);
    float r_proj = m_r0+t*(m_r1-m_r0);
//    p_proj = m_p0+t*(m_p1-m_p0);
    return exp(-dot(p - p_proj, p - p_proj) / r_proj / r_proj);
}

glm::vec3 SegmentFunction::EvalDev(glm::vec3 p) const
{
    return EvalDevFiniteDiff(p);

    float t = dot(normalize(m_p1-m_p0),p-m_p0) / length(m_p1-m_p0);

    if(t<0)
        return p - m_p0;
    if(t>=1)
        return p - m_p1;

    vec3 p_proj = m_p0+t*(m_p1-m_p0);
    return p-p_proj;
}



CircleFunction::CircleFunction(const glm::vec3& centre, const float radius_big, const float radius_small, const glm::vec3& normal):
    m_centre(centre),
    m_radius_big(radius_big),
    m_radius_small(radius_small),
    m_normal(normal)
{

}

CircleFunction::~CircleFunction()
{

}

float CircleFunction::Eval(glm::vec3 p) const
{
    float h = dot(p-m_centre, m_normal) / dot(m_normal, m_normal);
    vec3 p_plane_proj = p - h * m_normal;

    vec3 p_proj = m_centre + m_radius_big*normalize(p_plane_proj - m_centre);
    float d2 = dot(p - p_proj, p - p_proj);

    return exp(-d2 / m_radius_small / m_radius_small);
}

glm::vec3 CircleFunction::EvalDev(glm::vec3 p) const
{
    float h = dot(p-m_centre, m_normal) / dot(m_normal, m_normal);
    vec3 p_plane_proj = p - h * m_normal;

    vec3 p_proj = m_centre + m_radius_big*normalize(p_plane_proj - m_centre);

    return p_proj - p;
}


#define FASTFLOOR(x) ( ((x)>0) ? ((int)x) : (((int)x)-1) )

/*
 * Permutation table. This is just a random jumble of all numbers 0-255,
 * repeated twice to avoid wrapping the index at 255 for each lookup.
 * This needs to be exactly the same for all instances on all platforms,
 * so it's easiest to just keep it as static explicit data.
 * This also removes the need for any initialisation of this class.
 */
unsigned char perm[512] = {151,160,137,91,90,15,
  131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
  190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
  88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
  77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
  102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
  135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
  5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
  223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
  129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
  251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
  49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
  138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
  151,160,137,91,90,15,
  131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
  190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
  88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
  77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
  102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
  135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
  5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
  223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
  129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
  251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
  49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
  138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
};


double  grad3( int hash, double x, double y , double z ) {
    int h = hash & 15;     // Convert low 4 bits of hash code into 12 simple
    double u = h<8 ? x : y; // gradient directions, and compute dot product.
    double v = h<4 ? y : h==12||h==14 ? x : z; // Fix repeats at h = 12 to 15
    return ((h&1)? -u : u) + ((h&2)? -v : v);
}


// 3D simplex noise
double snoise3(double x, double y, double z) {

    // Simple skewing factors for the 3D case
    #define F3 0.333333333
    #define G3 0.166666667

    double n0, n1, n2, n3; // Noise contributions from the four corners

    // Skew the input space to determine which simplex cell we're in
    double s = (x+y+z)*F3; // Very nice and simple skew factor for 3D
    double xs = x+s;
    double ys = y+s;
    double zs = z+s;
    int i = FASTFLOOR(xs);
    int j = FASTFLOOR(ys);
    int k = FASTFLOOR(zs);

    double t = (double)(i+j+k)*G3;
    double X0 = i-t; // Unskew the cell origin back to (x,y,z) space
    double Y0 = j-t;
    double Z0 = k-t;
    double x0 = x-X0; // The x,y,z distances from the cell origin
    double y0 = y-Y0;
    double z0 = z-Z0;

    // For the 3D case, the simplex shape is a slightly irregular tetrahedron.
    // Determine which simplex we are in.
    int i1, j1, k1; // Offsets for second corner of simplex in (i,j,k) coords
    int i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords

    /* This code would benefit from a backport from the GLSL version! */
    if(x0>=y0) {
        if(y0>=z0)
        { i1=1; j1=0; k1=0; i2=1; j2=1; k2=0; } // X Y Z order
        else if(x0>=z0) { i1=1; j1=0; k1=0; i2=1; j2=0; k2=1; } // X Z Y order
        else { i1=0; j1=0; k1=1; i2=1; j2=0; k2=1; } // Z X Y order
    }
    else { // x0<y0
        if(y0<z0) { i1=0; j1=0; k1=1; i2=0; j2=1; k2=1; } // Z Y X order
        else if(x0<z0) { i1=0; j1=1; k1=0; i2=0; j2=1; k2=1; } // Y Z X order
        else { i1=0; j1=1; k1=0; i2=1; j2=1; k2=0; } // Y X Z order
    }

    // A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
    // a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
    // a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
    // c = 1/6.

    double x1 = x0 - i1 + G3; // Offsets for second corner in (x,y,z) coords
    double y1 = y0 - j1 + G3;
    double z1 = z0 - k1 + G3;
    double x2 = x0 - i2 + 2.0f*G3; // Offsets for third corner in (x,y,z) coords
    double y2 = y0 - j2 + 2.0f*G3;
    double z2 = z0 - k2 + 2.0f*G3;
    double x3 = x0 - 1.0f + 3.0f*G3; // Offsets for last corner in (x,y,z) coords
    double y3 = y0 - 1.0f + 3.0f*G3;
    double z3 = z0 - 1.0f + 3.0f*G3;

    // Wrap the integer indices at 256, to avoid indexing perm[] out of bounds
    int ii = i % 256;
    int jj = j % 256;
    int kk = k % 256;

    // Calculate the contribution from the four corners
    double t0 = 0.6f - x0*x0 - y0*y0 - z0*z0;
    if(t0 < 0.0f) n0 = 0.0f;
    else {
        t0 *= t0;
        n0 = t0 * t0 * grad3(perm[ii+perm[jj+perm[kk]]], x0, y0, z0);
    }

    double t1 = 0.6f - x1*x1 - y1*y1 - z1*z1;
    if(t1 < 0.0f) n1 = 0.0f;
    else {
        t1 *= t1;
        n1 = t1 * t1 * grad3(perm[ii+i1+perm[jj+j1+perm[kk+k1]]], x1, y1, z1);
    }

    double t2 = 0.6f - x2*x2 - y2*y2 - z2*z2;
    if(t2 < 0.0f) n2 = 0.0f;
    else {
        t2 *= t2;
        n2 = t2 * t2 * grad3(perm[ii+i2+perm[jj+j2+perm[kk+k2]]], x2, y2, z2);
    }

    double t3 = 0.6f - x3*x3 - y3*y3 - z3*z3;
    if(t3<0.0f) n3 = 0.0f;
    else {
        t3 *= t3;
        n3 = t3 * t3 * grad3(perm[ii+1+perm[jj+1+perm[kk+1]]], x3, y3, z3);
    }

    // Add contributions from each corner to get the final noise value.
    // The result is scaled to stay just inside [-1,1]
    return 32.0f * (n0 + n1 + n2 + n3); // TODO: The scale factor is preliminary!
}


double noise(const vec3& x,int octave,double persistency)
{
    double temp=0.0;
    double frequency=1.0;
    double current_persistency=1.0;
    for(int k=0;k<octave;k++)
    {

        temp += current_persistency*(0.5+0.5*snoise3(x.x*frequency,x.y*frequency,x.z*frequency));
        frequency *= 2.0;
        current_persistency *= persistency;
    }

    //normalize
    double epsilon=1e-5;
    if( std::abs(current_persistency-1)>epsilon )
        temp *= (1-persistency)/(1-current_persistency);

    return temp;
}






PerlinFunction::PerlinFunction(const glm::uint nb_octave, const float persistency, const float scale, const vec3 &displacement):
    m_nb_octave(nb_octave),
    m_persistency(persistency),
    m_scale(scale),
    m_displacement(displacement)
{

}

PerlinFunction::~PerlinFunction()
{

}

float PerlinFunction::Eval(glm::vec3 p) const
{
    return noise(m_scale * p + m_displacement, m_nb_octave, m_persistency);
}

glm::vec3 PerlinFunction::EvalDev(glm::vec3 p) const
{
    return EvalDevFiniteDiff(p);
}





BarthFunction::BarthFunction()
{

}

BarthFunction::~BarthFunction()
{

}

float BarthFunction::Eval(glm::vec3 p) const
{
    float phi = (1.0 + glm::sqrt(5.0))*0.5;
    float w = 1.0;
    return 4.0*(p.x*p.x*phi*phi -p.y*p.y)*(p.y*p.y*phi*phi -p.z*p.z)*(p.z*p.z*phi*phi -p.x*p.x)-(1.0+2.0*phi)*(p.x*p.x+p.y*p.y+p.z*p.z-w*w)*(p.x*p.x+p.y*p.y+p.z*p.z-w*w)*w*w;
}

glm::vec3 BarthFunction::EvalDev(glm::vec3 p) const
{
    return EvalDevFiniteDiff(p);
}


