#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include <shader.h> // Help to load shaders from files

// Include GLEW : Always include it before glfw.h et gl.h :)
#include <GL/glew.h>    // OpenGL Extension Wrangler Library : http://glew.sourceforge.net/
#include <GL/glfw.h>    // Window, keyboard, mouse : http://www.glfw.org/

#include <glm/glm.hpp>  // OpenGL Mathematics : http://glm.g-truc.net/0.9.5/index.html
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/ext.hpp>

#include "GLFW_define.h"
#include "Mesh.h"
#include "DataSet.h"
#include "ImplicitFunction.h"
#include "MeshHE.h"
#include "Object.h"


// Window size :
#define WIDTH 1000.0f
#define HEIGHT 800.0f

using namespace glm;
using namespace std;


static void cursor_position_callback(int xpos, int ypos);
static void mouse_wheel_callback(int wheel_pos);
void mouse_button_callback(int button, int action);

int prev_wheel = -1;
int prev_x = -1;
int prev_y = -1;
bool mouse_pressed = false;

float view_angle = 45.0f;

mat4 view_matrix = lookAt(vec3(1.0, 0.5, 1.0), vec3(0.0), vec3(0.0, 1.0, 0.0));
mat4 projection_matrix = perspective(view_angle, WIDTH / HEIGHT, 0.1f, 1000.0f);


void view_control(mat4& view_matrix, float dx);

int main()
{

    cout << "Starting program..." << endl;

    //==================================================
    //============= Creation de la fenetre =============
    //==================================================

    // GLFW initialization
    if( !glfwInit() )
    {
        cerr << "Failed to initialize GLFW!" << endl;
        exit(EXIT_FAILURE);
    }

    glfwOpenWindowHint(GLFW_FSAA_SAMPLES, 4); // Anti Aliasing
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MAJOR, 3); // OpenGL 3.1
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MINOR, 1);

    // Window and OpenGL conetxt creation
    if( !glfwOpenWindow(WIDTH, HEIGHT, 0,0,0,0, 32,0, GLFW_WINDOW ) )
//    if( !glfwOpenWindow(WIDTH, HEIGHT, 0,0,0,0, 32,0, GLFW_FULLSCREEN ) )
    {
        cerr << "GLFW failed to open OpenGL window!" << endl;
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    // GLFW Settings
    glfwSetWindowTitle( "TP 3A Ensimag - MMMIS - Projet" );
    glfwEnable( GLFW_STICKY_KEYS );

    // GLEW Initialization
    if (glewInit() != GLEW_OK) {
        cerr << "Failed to intialize GLEW:!" << endl;
        exit(EXIT_FAILURE);
    }



    // Soft- and Firm-ware checkings
    const GLubyte* renderer = glGetString (GL_RENDERER);
    cout << "GPU : " << renderer << endl;

    const GLubyte* version = glGetString (GL_VERSION);
    cout << "OpenGL Driver : " << version << endl;


    cout << endl;


    //==================================================
    //================= Initialization =================
    //==================================================

    cout << "Initializing..." << endl;

    glfwSetMouseButtonCallback(mouse_button_callback);
    glfwSetMousePosCallback(cursor_position_callback);
    glfwSetMouseWheelCallback(mouse_wheel_callback);



    //-------------------------------------------------
    // OpenGL Initialization

    //glClearColor(0.1, 0.1, 0.1, 1.0);       /// Dark Back ground
    glClearColor(1.0, 1.0, 1.0, 1.0);       /// Light Back ground
    glEnable(GL_DEPTH_TEST);



    //-------------------------------------------------
    // Shader program initialization

    GLuint programID = LoadShaders("../shader/vertex.glsl", "../shader/fragment.glsl");
    cout << "programID = " << programID << endl;


    //-------------------------------------------------
    // Data arrays Initialization

    // Mesh creation
    //Mesh m("../models/Horse.off");

    //BarthFunction bf;
    //float x = 1.8;
    //Mesh m; m.CreateIsoSurface(m, bf, -0.2, -x, x, -x, x, -x ,x, 100, 100, 100);

    DataSet ds("../data/bear.data");

    ds.ComputeTangentPlanes();
    ds.ComputeEMST();
    ds.AddKNeighborsEdges();
    // ds.AssignCostOnEdges();
    // ds.AssignTangentPlanesOrientation();

    //DistanceFunction f(ds);
    glm::vec3 c(0, 0, 0);
    SphereFunction fs(c, 1);

    //std::cout << ds.minX() << std::endl;

    Mesh m; m.CreateIsoSurface(m, fs, 0, ds.minX(), ds.maxX(), ds.minY(), ds.maxY(), ds.minZ(), ds.maxZ(), 20, 20, 20);
    //glm::vec3 test(0, 0, 1);
    //std::cout << "Sphere function: " << fs.Eval(test) << ", our function: " << f.Eval(test) << std::endl;

    /** Mesh creation from data set and iso function **/
    // double minX = ds.minX(); double minY = ds.minY(); double minZ = ds.minZ();
    // double maxX = ds.maxX(); double maxY = ds.maxY(); double maxZ = ds.maxZ();
    // const unsigned int resX = 10;
    // const unsigned int resY = 10;
    // const unsigned int resZ = 10;
    // minX-=2.*(maxX-minX)/resX; minY-=2.*(maxY-minY)/resY; minZ-=2.*(maxZ-minZ)/resZ;
    // maxX+=2.*(maxX-minX)/resX; maxY+=2.*(maxY-minY)/resY; maxZ+=2.*(maxZ-minZ)/resZ;

    //Mesh m(fs, minX, maxX, minY, maxY, minZ, maxZ, resX, resY, resZ);
    //printf("Mesh Created with %i point positions and %i faces\n", m.NbVertices(), m.NbFaces());

    m.Normalize();
    m.ComputeNormals();
    m.ColorFromNormals();

    // Half edge conversion
    MeshHE m_he(m);

    // Object Generation
    Object o;
    o.GenBuffers();
    o.SetMesh(&m_he);
//    o.SetMesh(&m);
    o.SetShader(programID);


    //-------------------------------------------------
    // MVP matrices initialization

    GLuint PmatrixID = glGetUniformLocation(programID, "ProjectionMatrix");
    cout << "PmatrixID = " << PmatrixID << endl;

    GLuint VmatrixID = glGetUniformLocation(programID, "ViewMatrix");
    cout << "VmatrixID = " << VmatrixID << endl;



    cout << "Initializing done." << endl;
    cout << endl;



    //==================================================
    //==================== Main Loop ===================
    //==================================================


    cout << "Starting main loop..." << endl;

    double init_time = glfwGetTime();
    double prec_time = init_time;
    double cur_time = init_time;
    double speed = 2.0;

    do{
        // Clearing Viewport
        glClear( GL_COLOR_BUFFER_BIT );
        glClear( GL_DEPTH_BUFFER_BIT );


        //==================================================
        //================== Computations ==================
        //==================================================

        prec_time = cur_time;
        cur_time = glfwGetTime() - init_time;
        float delta_time = cur_time - prec_time;

        view_control(view_matrix, speed * delta_time);


        //==================================================
        //===================== Drawing ====================
        //==================================================

        o.Draw(view_matrix, projection_matrix, VmatrixID, PmatrixID);

        glfwSwapBuffers();


        //==================================================
        //================== Stats Display =================
        //==================================================

//        cout.precision(2);
//        cout<< fixed  << "FPS: " << 1.0 / delta_time << "\t---\tElapsed time: " << setw(8) << cur_time << "\ts\r" << flush;
    }
    while( glfwGetKey( GLFW_KEY_ESC ) != GLFW_PRESS &&
           glfwGetWindowParam( GLFW_OPENED )        );

    // Closing the window
    glfwTerminate();

    cout << endl << endl << "Main loop ended." << endl;


    cout << "Program ended." << endl;


    return EXIT_SUCCESS;
}




void view_control(mat4& view_matrix, float dx)
{
    if (glfwGetKey( GLFW_KEY_LSHIFT ) == GLFW_PRESS)
    {
        dx /= 10.0;
    }

    if (glfwGetKey( GLFW_KEY_UP ) == GLFW_PRESS)
    {
        vec4 axis = vec4(1.0, 0.0, 0.0, 0.0);
        axis = inverse(view_matrix) * axis;
        view_matrix = rotate(view_matrix, dx * 180.0f, vec3(axis));
    }
    if (glfwGetKey( GLFW_KEY_DOWN ) == GLFW_PRESS)
    {
        vec4 axis = vec4(1.0, 0.0, 0.0, 0.0);
        axis = inverse(view_matrix) * axis;
        view_matrix = rotate(view_matrix, -dx * 180.0f, vec3(axis));
    }
    if (glfwGetKey( GLFW_KEY_RIGHT ) == GLFW_PRESS)
    {
        vec4 axis = vec4(0.0, 1.0, 0.0, 0.0);
        axis = inverse(view_matrix) * axis;
        view_matrix = rotate(view_matrix, dx * 180.0f, vec3(axis));
    }
    if (glfwGetKey( GLFW_KEY_LEFT ) == GLFW_PRESS)
    {
        vec4 axis = vec4(0.0, 1.0, 0.0, 0.0);
        axis = inverse(view_matrix) * axis;
        view_matrix = rotate(view_matrix, -dx * 180.0f, vec3(axis));
    }
    if (glfwGetKey( GLFW_KEY_PAGEUP ) == GLFW_PRESS)
    {
        vec4 axis = vec4(0.0, 0.0, 1.0, 0.0);
        axis = inverse(view_matrix) * axis;
        view_matrix = rotate(view_matrix, dx * 180.0f, vec3(axis));
    }
    if (glfwGetKey( GLFW_KEY_PAGEDOWN ) == GLFW_PRESS)
    {
        vec4 axis = vec4(0.0, 0.0, 1.0, 0.0);
        axis = inverse(view_matrix) * axis;
        view_matrix = rotate(view_matrix, -dx * 180.0f, vec3(axis));
    }

    if (glfwGetKey( GLFW_KEY_Z ) == GLFW_PRESS)
    {
        vec3 pos = vec3(view_matrix * vec4(0,0,0,1));
        vec4 axis = vec4(0.0, 0.0, 1.0, 0.0) * dx * length(pos) * 0.5;
        axis = inverse(view_matrix) * axis;
        view_matrix = translate(view_matrix, vec3(axis));
    }
    if (glfwGetKey( GLFW_KEY_S ) == GLFW_PRESS)
    {
        vec3 pos = vec3(view_matrix * vec4(0,0,0,1));
        vec4 axis = vec4(0.0, 0.0, 1.0, 0.0) * (-dx) * length(pos) * 0.5;
        axis = inverse(view_matrix) * axis;
        view_matrix = translate(view_matrix, vec3(axis));
    }
    if (glfwGetKey( GLFW_KEY_Q) == GLFW_PRESS)
    {
        vec4 axis = vec4(-1.0, 0.0, 0.0, 0.0) * dx;
        axis = inverse(view_matrix) * axis;
        view_matrix = translate(view_matrix, vec3(axis));
    }
    if (glfwGetKey( GLFW_KEY_D ) == GLFW_PRESS)
    {
        vec4 axis = vec4(-1.0, 0.0, 0.0, 0.0) * (-dx);
        axis = inverse(view_matrix) * axis;
        view_matrix = translate(view_matrix, vec3(axis));
    }
    if (glfwGetKey( GLFW_KEY_A ) == GLFW_PRESS)
    {
        vec4 axis = vec4(0.0, 1.0, 0.0, 0.0) * dx;
        axis = inverse(view_matrix) * axis;
        view_matrix = translate(view_matrix, vec3(axis));
    }
    if (glfwGetKey( GLFW_KEY_E ) == GLFW_PRESS)
    {
        vec4 axis = vec4(0.0, 1.0, 0.0, 0.0) * (-dx);
        axis = inverse(view_matrix) * axis;
        view_matrix = translate(view_matrix, vec3(axis));
    }
}

static void cursor_position_callback(int xpos, int ypos)
{
    float dx = float(xpos - prev_x);
    float dy = float(ypos - prev_y);
    prev_x = xpos;
    prev_y = ypos;

    if(!mouse_pressed)
        return;


    if(glfwGetKey(GLFW_KEY_LSHIFT) == GLFW_PRESS)
    {
        dx *= 0.1;
        dy *= 0.1;
    }

    if(glfwGetMouseButton(GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
    {
        vec4 axis0 = vec4(1.0, 0.0, 0.0, 0.0);
        axis0 = inverse(view_matrix) * axis0;
        view_matrix = rotate(view_matrix, dy, vec3(axis0));

        vec4 axis1 = vec4(0.0, 1.0, 0.0, 0.0);
        axis1 = inverse(view_matrix) * axis1;
        view_matrix = rotate(view_matrix, dx, vec3(axis1));
    }

    if(glfwGetMouseButton(GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
    {
        vec4 axis0 = vec4(1.0, 0.0, 0.0, 0.0) * dx * view_angle * 0.0001;
        axis0 = inverse(view_matrix) * axis0;
        view_matrix = translate(view_matrix, vec3(axis0));

        vec4 axis1 = vec4(0.0, -1.0, 0.0, 0.0) * dy * view_angle * 0.0001;
        axis1 = inverse(view_matrix) * axis1;
        view_matrix = translate(view_matrix, vec3(axis1));
    }

    if(glfwGetMouseButton(GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS)
    {
        vec4 axis2 = vec4(0.0, 0.0, 1.0, 0.0);
        axis2 = inverse(view_matrix) * axis2;
        view_matrix = translate(view_matrix, vec3(axis2)*dy * 0.01f);
        view_matrix = rotate(view_matrix, dx, vec3(axis2));
    }

}

static void mouse_wheel_callback(int wheel_pos)
{
    int delta = wheel_pos - prev_wheel;
    prev_wheel = wheel_pos;
    if(delta > 0)
    {
        if(glfwGetKey(GLFW_KEY_LSHIFT) == GLFW_PRESS)
        {
            view_angle *= 1.01;
        }
        else
        {
            view_angle *= 1.1;
        }
    }
    else
    {
        if(glfwGetKey(GLFW_KEY_LSHIFT) == GLFW_PRESS)
        {
            view_angle /= 1.01;
        }
        else
        {
            view_angle /= 1.1;
        }
    }

    view_angle = glm::min(view_angle, 180.0f);

    projection_matrix = perspective(view_angle, WIDTH / HEIGHT, 0.1f, 1000.0f);
}

void mouse_button_callback(int button, int action)
{
    mouse_pressed = !mouse_pressed;
}
