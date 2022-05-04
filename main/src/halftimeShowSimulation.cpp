/*
Author: Maria Padilla Fuentes
Description:
Makes football field with the flying UAVs
*/

#include <cmath>
#include <math.h>
#include <cstdlib>
#include <GL/glut.h>
#include <chrono>
#include <thread>
#include <array>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "iomanip"
#include <random>
#include <GL/glut.h>
#include "Bitmap.h"
#define checkImageWidth 1200
#define checkImageHeight 579
#define ESC 27

GLuint texture;

struct Image {
    unsigned long sizeX;
    unsigned long sizeY;
    char *data;
};
typedef struct Image Image;
BMP inBitmap;


GLubyte checkImage[checkImageWidth][checkImageHeight][3];
void makeCheckImage(void) {
    int i, j, c;
    for (i = 0; i < checkImageWidth; i++) {
        for (j = 0; j < checkImageHeight; j++) {
            c = ((((i & 0x8) == 0) ^ ((j & 0x8) == 0))) * 255;
            checkImage[i][j][0] = (GLubyte)c;
            checkImage[i][j][1] = (GLubyte)c;
            checkImage[i][j][2] = (GLubyte)c;
        }
    }
}

struct XYZ_Coordinates
{
    XYZ_Coordinates() = default;
    XYZ_Coordinates(double x, double y, double z)
                    :x{x}, y{y}, z{z} {}
    double x{0},y{0},z{0};
};

struct float_XYZ_Coordinates
{
    float_XYZ_Coordinates() = default;
    float_XYZ_Coordinates(float x, float y, float z)
                        :   x{x}, y{y}, z{z} {}
    float x{0},y{0},z{0};
};

struct UAVInfo
{
    UAVInfo() = default;
    UAVInfo(XYZ_Coordinates position, XYZ_Coordinates velocity, bool atSphere)
            :   position{position}, velocity{velocity}, atSphere{atSphere} {}
    XYZ_Coordinates position{0,0,0}, velocity{0,0,0};
    bool atSphere{true};
};

//----------------------------------------------------------------------
// Returns measurement in yards to meters
//
// Convert yards (used for dimensions of footall field) into meters
// (which are the dimesions for openGL functions)
//----------------------------------------------------------------------
float convertToMeters(float yardValue)
{
    return yardValue*.9144;
}

static float footBallFieldWidth_x = convertToMeters((53*3+1)/3);
static float footBallFieldLength_y = convertToMeters(100);
static const unsigned int numberOfUAVs = 15;
static const unsigned int numElements = 7;
std::array<UAVInfo, numberOfUAVs> UAVsArray;
XYZ_Coordinates origin = XYZ_Coordinates{footBallFieldWidth_x/2, footBallFieldLength_y/2, 0};
XYZ_Coordinates sphereLocation = XYZ_Coordinates{origin.x, origin.y, origin.z + 50};

bool atSphere{true};
bool decreaseSpeed{false};


float_XYZ_Coordinates eyePosition = float_XYZ_Coordinates{-90, 45, 35};
float_XYZ_Coordinates centerPosition = float_XYZ_Coordinates{footBallFieldWidth_x/2,  
                                                            footBallFieldLength_y/2, 
                                                            25};
float angleOfRotation = 0.00;

const int rcvSize = 16 * 6; 
double* rcvbuffer = new double[rcvSize];
double sendBuffer[numElements];

static double mass = 1; 
static double kSpring = 1;
static double timeStep = 0.1; 
static double maxForce = 20;
static double gravity = -10;

std::random_device rd; 
std::uniform_int_distribution<> distr(1, 3); 
std::mt19937 eng(rd()); 


//----------------------------------------------------------------------
// Display the Football field 
//
// Using bitmag image given and 'ECE_Matrix.h', image used as 
// footbal field
//----------------------------------------------------------------------
void displayFootballField()
{
    inBitmap.read("AmFBfield.bmp");
    makeCheckImage();
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); //scale linearly when image bigger than texture
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); //scale linearly when image smalled than texture
    glTexImage2D(GL_TEXTURE_2D, 0, 3, inBitmap.bmp_info_header.width, inBitmap.bmp_info_header.height, 0,
        GL_BGR_EXT, GL_UNSIGNED_BYTE, &inBitmap.data[0]);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texture);

    glPushMatrix(); //// yardValue*.9144;
        glBegin(GL_QUADS);
        glTexCoord2f(0, 0);
                            glVertex3f(0-convertToMeters(3), 
                                        footBallFieldLength_y+convertToMeters(10),
                                        0.0);
        glTexCoord2f(0, 1);
                            glVertex3f(footBallFieldWidth_x+convertToMeters(3), 
                                        footBallFieldLength_y+convertToMeters(10), 
                                        0.0);
        glTexCoord2f(1, 1);
                            glVertex3f(footBallFieldWidth_x+convertToMeters(3), 
                                        0-convertToMeters(10), 
                                        0.0);
        glTexCoord2f(1, 0);
                            glVertex3f(0-convertToMeters(3), 
                                        0-convertToMeters(10), 
                                        0.0);
        glEnd();
    glPopMatrix();
}

//----------------------------------------------------------------------
// Draw the wired sphere
//
// Sphere placed in the middle of the football field, 50 meters
// above the field 
//----------------------------------------------------------------------
void drawWireSphere()
{
    glDisable(GL_TEXTURE_2D);
    glPushMatrix();
        glColor3f(0.0, 0.0, 0.0);
        glTranslatef(sphereLocation.x, 
                        sphereLocation.y, 
                        sphereLocation.z);
        glutWireSphere(10, 12, 12);
    glPopMatrix();
    glEnable(GL_TEXTURE_2D);
}

//----------------------------------------------------------------------
// Draw UAVs in their initial position
//
// UAVs drawn in the 0,25,and 50 lines, taking into account that 
// yards have to be converted to meters
//----------------------------------------------------------------------
void startingUAVsPos()
{
    unsigned int n = 0;
    unsigned int column = 0;
    for(auto &currentUAV : UAVsArray)
    {
        currentUAV.position.y = convertToMeters(25*column);
        currentUAV.position.x = (n%3)*footBallFieldWidth_x/2;
        n++;
        if (n %3 == 0)
            column++;
    }
}

//----------------------------------------------------------------------
// drawUAVs in their updated position 
//
// UAVs drawn accoriding to the color rules and shape rules specified
// in the pdf
//----------------------------------------------------------------------
void drawUAVs()
{
    for(const auto &currentUAV : UAVsArray)
    {
        glDisable(GL_TEXTURE_2D);
        glColor3f(1.0, 1.0, 0.0); //UAV color is determined using the first letter of your last name: yellow F,P,Y
        glPushMatrix();
            glTranslatef(currentUAV.position.x, 
                        currentUAV.position.y, 
                        currentUAV.position.z);
            glScalef(1.0, 1.0, 1.0);
            glutSolidCone(1.0, 1.0, 20, 20); ////UAV shape is determined using the first letter of your first name: Cone C,M,V
        glPopMatrix();
        glEnable(GL_TEXTURE_2D);
    }
}

//----------------------------------------------------------------------
// Check for Collisions
//
// Check for collisions, and upon collisions, swap the velocites
// of the two UAVs that are colliding
//----------------------------------------------------------------------
void ifCollisionSwapVelocities()
{   
    int n = 0;
    int arrSize = UAVsArray.size();
    while ( n < arrSize )
    {
        int m = n;
        while ( m  < arrSize )
        {
            double uavsDistance = sqrt(
                                    pow(UAVsArray.at(n).position.x - UAVsArray.at(m).position.x, 2) +
                                    pow(UAVsArray.at(n).position.y - UAVsArray.at(m).position.y, 2) +
                                    pow(UAVsArray.at(n).position.z - UAVsArray.at(m).position.z, 2));
            
            if( uavsDistance < 0.51 )
            {
                XYZ_Coordinates tempVelocity = UAVsArray.at(n).velocity;
                UAVsArray.at(n).velocity = UAVsArray.at(m).velocity;
                UAVsArray.at(m).velocity = tempVelocity;
            }
            m+=1;
        }
        n+=1;
    }
}

//----------------------------------------------------------------------
// Calculate UAVs location
//
// using MPI (arguement being the current rank)
// calculates new UAVs location using conditions set in assignment pdf
//----------------------------------------------------------------------
void CalcualteUAVsLocation(unsigned int currUAVsRank)
{    
    XYZ_Coordinates force;
    double velocityVectorMagnitude = sqrt(pow(UAVsArray.at(currUAVsRank).velocity.x, 2) 
                + pow(UAVsArray.at(currUAVsRank).velocity.y, 2) 
                + pow(UAVsArray.at(currUAVsRank).velocity.z,2));
    
    
    if( velocityVectorMagnitude <= 5)
    {
        double distance = sqrt(pow(sphereLocation.x -  UAVsArray.at(currUAVsRank).position.x,2) 
                            + pow(sphereLocation.y - UAVsArray.at(currUAVsRank).position.y, 2) 
                            + pow(sphereLocation.z - UAVsArray.at(currUAVsRank).position.z,2));;
        
        double forceVectorMagnitude{0};
        
        if (UAVsArray.at(currUAVsRank).atSphere)
        {
            if(velocityVectorMagnitude >= 2)
                forceVectorMagnitude = 0;
            else
               forceVectorMagnitude =  -kSpring* (9 - distance);
            
            if (distance < 11)
                UAVsArray.at(currUAVsRank).atSphere = false;
        }
        else
        {
            forceVectorMagnitude =  -kSpring*(9 - distance);
            if(!decreaseSpeed)
                forceVectorMagnitude *= 0.8;
        }

        if (forceVectorMagnitude < -maxForce)
            forceVectorMagnitude = -maxForce;
        
        else if (forceVectorMagnitude > maxForce)
            forceVectorMagnitude = maxForce;

        if(forceVectorMagnitude > 9)
            forceVectorMagnitude = 9;
        
        else if (forceVectorMagnitude < -9)
            forceVectorMagnitude = -9;

        XYZ_Coordinates unitVec = XYZ_Coordinates{
                            (sphereLocation.x - UAVsArray.at(currUAVsRank).position.x)/(distance), 
                            (sphereLocation.y - UAVsArray.at(currUAVsRank).position.y)/(distance),  
                            (sphereLocation.z - UAVsArray.at(currUAVsRank).position.z)/(distance)};
        
        force = XYZ_Coordinates{unitVec.x*forceVectorMagnitude, 
                        unitVec.y*forceVectorMagnitude, 
                        unitVec.z*forceVectorMagnitude};

        if(!UAVsArray.at(currUAVsRank).atSphere)
        {  
            double unitx = distr(eng);
            double unity = distr(eng);
            double unitz = distr(eng);

            double counter = distr(eng);

            bool whenToBreak{false};

            XYZ_Coordinates orthoVec;

            int n = 0;
            while (n < 3)
            {
                switch( (int(n + counter)) %3 )
                {
                    case 0:
                        if( !(unitVec.x == 0) )
                            orthoVec = {(-unity*unitVec.y - unitz*unitVec.z)/unitVec.x, unity, unitz}; whenToBreak = true;
                        break;
                    case 1:
                        if( !(unitVec.y == 0) )
                            orthoVec = {(-unitx*unitVec.x - unitz*unitVec.z)/unitVec.y, unity, unitz}; whenToBreak = true;
                        break;
                    case 2:
                        if( !(unitVec.z == 0) )
                            orthoVec = {(-unitx*unitVec.x - unity*unitVec.y)/unitVec.z, unity, unitz}; whenToBreak = true;
                        break;
                }
                if(whenToBreak)
                    break;
                n++;
            }
            if(!whenToBreak)
                orthoVec = {1,1,1};
            
            double mag =  (1/(sqrt(pow(orthoVec.x, 2) 
                            + pow(orthoVec.y, 2) 
                            + pow(orthoVec.z,2))));
            
            XYZ_Coordinates orthogonalForce = XYZ_Coordinates{
                                    orthoVec.x*mag, 
                                    orthoVec.y*mag, 
                                    orthoVec.z*mag};
            
            force.x += orthogonalForce.x;
            force.y += orthogonalForce.y;
            force.z += orthogonalForce.z;
        }
    }
    else 
    {
        force = XYZ_Coordinates{UAVsArray.at(currUAVsRank).velocity.x/velocityVectorMagnitude, 
                        UAVsArray.at(currUAVsRank).velocity.y/velocityVectorMagnitude, 
                        UAVsArray.at(currUAVsRank).velocity.z/velocityVectorMagnitude};
        
        force = XYZ_Coordinates{force.x*-10, force.y*-10, force.z*-10};
    }
    
    force.z += gravity;

    double forceVectorMagnitude = sqrt(pow(force.x,2) 
                        + pow(force.y, 2) 
                        + pow(force.z, 2));

    force.z -= gravity;

     UAVsArray.at(currUAVsRank).position.x = UAVsArray.at(currUAVsRank).position.x 
                                    + UAVsArray.at(currUAVsRank).velocity.x*timeStep 
                                    + 0.5*force.x/mass*pow(timeStep,2);
    
    UAVsArray.at(currUAVsRank).position.y = UAVsArray.at(currUAVsRank).position.y 
                                    + UAVsArray.at(currUAVsRank).velocity.y*timeStep 
                                    + 0.5*force.y/mass*pow(timeStep,2);
    
    UAVsArray.at(currUAVsRank).position.z = UAVsArray.at(currUAVsRank).position.z 
                                    + UAVsArray.at(currUAVsRank).velocity.z*timeStep 
                                    + 0.5*force.z/mass*pow(timeStep,2);

    double beforeMag = sqrt(pow(UAVsArray.at(currUAVsRank).velocity.x, 2) 
                        + pow(UAVsArray.at(currUAVsRank).velocity.y, 2) 
                        + pow(UAVsArray.at(currUAVsRank).velocity.z,2));
    
    UAVsArray.at(currUAVsRank).velocity.x = UAVsArray.at(currUAVsRank).velocity.x 
                                    + force.x/mass*timeStep;
    
    UAVsArray.at(currUAVsRank).velocity.y = UAVsArray.at(currUAVsRank).velocity.y 
                                    + force.y/mass*timeStep;
    
    UAVsArray.at(currUAVsRank).velocity.z = UAVsArray.at(currUAVsRank).velocity.z 
                                    + force.z/mass*timeStep;
    
    double afterMag = sqrt(pow(UAVsArray.at(currUAVsRank).velocity.x, 2) 
                    + pow(UAVsArray.at(currUAVsRank).velocity.y, 2) 
                    + pow(UAVsArray.at(currUAVsRank).velocity.z,2));

    decreaseSpeed = afterMag < beforeMag;
}

//----------------------------------------------------------------------
// Draw the entire scene
//
// using openGL make scene appear
// (including the football field, the sphere & UAVs)
//----------------------------------------------------------------------
void drawEntireScene() 
{
    displayFootballField();
    drawWireSphere();
    drawUAVs();
}

//----------------------------------------------------------------------
// Update Recieve Array
//
// Update recieve array to hold correct UAV info 
// (including UAV's x,y,z,vx,vy, and vz)
//----------------------------------------------------------------------
void setUAVsArrayToRcvBuffer(double arr[])
{
    unsigned int n{1};
    for(auto &currentUAV : UAVsArray)
    {
        unsigned int offset = (numElements * n);
        currentUAV.position.x = arr[ 0 + (offset) ];
        currentUAV.position.y = arr[ 1 + (offset) ];
        currentUAV.position.z = arr[ 2 + (offset) ];
        currentUAV.velocity.x = arr[ 3 + (offset) ];
        currentUAV.velocity.y = arr[ 4 + (offset) ];
        currentUAV.velocity.z = arr[ 5 + (offset) ];
        currentUAV.atSphere = arr[ 6 + (offset) ];
        n++;
    }
}

//----------------------------------------------------------------------
// Update Buffer Array
//
// Update buffer to send correct UAV info 
// (including UAV's x,y,z,vx,vy, and vz)
//----------------------------------------------------------------------
void sendRankUavInfo(double arr[], unsigned int rankUAVsNum)
{
    arr[0] = UAVsArray.at(rankUAVsNum).position.x;
    arr[1] = UAVsArray.at(rankUAVsNum).position.y;
    arr[2] = UAVsArray.at(rankUAVsNum).position.z;
    arr[3] = UAVsArray.at(rankUAVsNum).velocity.x;
    arr[4] = UAVsArray.at(rankUAVsNum).velocity.y;
    arr[5] = UAVsArray.at(rankUAVsNum).velocity.z;
    arr[6] = UAVsArray.at(rankUAVsNum).atSphere;
}

//----------------------------------------------------------------------
// timerFunction  - called whenever the timer fires
//----------------------------------------------------------------------
void timerFunction(int id)
{
    glutPostRedisplay();
    glutTimerFunc(100, timerFunction, 0);
}

//----------------------------------------------------------------------
// Reshape callback
//
// Window size has been set/changed to w by h pixels. Set the camera
// perspective to 45 degree vertical field of view, a window aspect
// ratio of w/h, a near clipping plane at depth 1, and a far clipping
// plane at depth 100. The viewport is the entire window.
//
//----------------------------------------------------------------------
void changeSize(int w, int h)
{
    float ratio = ((float)w) / ((float)h); // window aspect ratio
    glMatrixMode(GL_PROJECTION); // projection matrix is active
    glLoadIdentity(); // reset the projection
    gluPerspective(45.0, ratio, 0.1, 1000.0); // perspective transformation
    //gluPerspective(60.00, ratio, 0.1, 1000.0); // perspective transformation
    glMatrixMode(GL_MODELVIEW); // return to modelview mode
    glViewport(0, 0, w, h); // set viewport (drawing area) to entire window
}

//----------------------------------------------------------------------
// User-input callbacks
//
// processNormalKeys: change view of the entireSceen using
//  keys specified below
//----------------------------------------------------------------------
void processNormalKeys(unsigned char key, int xx, int yy)
{   
    if (key == ESC || key == 'q' || key == 'Q')
        exit(0);

    if (key == 'r' || key == 'R')
        eyePosition.x +=10;

    if (key == 'l' || key == 'L')
        eyePosition.x -=10;

    if (key == 'y' || key == 'Y')
        eyePosition.y +=10;

    if (key == 'T' || key == 't')
        eyePosition.y -=10;

    if (key == 'd' || key == 'D')
        eyePosition.z -=2.50;

    if (key == 'u' || key == 'U')
        eyePosition.z +=2.50;

    if (key == 'z' || key == 'Z')
        centerPosition.z -=1.00;

    if (key == 'x' || key == 'X')
        centerPosition.z +=1.00;

    if (key == 'w' || key == 'W')
        angleOfRotation +=2.50;
    
    if (key == 'e' || key == 'E')
        angleOfRotation -=2.50;

}

//----------------------------------------------------------------------
// Draw the entire scene
//
// We first update the camera location based on its distance from the
// origin and its direction.
//----------------------------------------------------------------------
void renderScene()
{
    glClearColor(   (67.00/255.00), 
                    (105.00/255.00), 
                    (171.00/255.00), 
                    1.0); 

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Reset transformations
    glLoadIdentity();

    gluLookAt(eyePosition.x, eyePosition.y, eyePosition.z, 
                centerPosition.x, centerPosition.y, centerPosition.z, 
                0.0, 0.0, 1.0);

    glMatrixMode(GL_MODELVIEW);

    drawEntireScene();

    glutSwapBuffers(); // Make it all visible

    MPI_Allgather(sendBuffer, numElements, MPI_DOUBLE, 
                rcvbuffer, numElements, MPI_DOUBLE, MPI_COMM_WORLD);
    
    setUAVsArrayToRcvBuffer(rcvbuffer);
}

//----------------------------------------------------------------------
// mainOpenGL  - standard GLUT initializations and callbacks
//----------------------------------------------------------------------
void mainOpenGL(int argc, char**argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(400, 400);

    glutCreateWindow(argv[0]);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    //GLfloat mat_specular[] = {0.5, 0.5, 0.5, 1.0};
    GLfloat mat_specular[] = { 0.75, 0.75, 0.75, 1.0 }; 
    GLfloat mat_shininess[] = { 50.0 };

    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

    glutReshapeFunc(changeSize);
    glutDisplayFunc(renderScene);
    glutKeyboardFunc(processNormalKeys);
    glutTimerFunc(100, timerFunction, 0);
    glutMainLoop();
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Main entry point determines rank of the process and follows the 
// correct program path
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int main(int argc, char**argv)

{
    startingUAVsPos();
    int numTasks, rank;

    int rc = MPI_Init(&argc, &argv);

    if (rc != MPI_SUCCESS) 
    {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int gsize = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &gsize);


    if (rank == 0) 
    {
        mainOpenGL(argc, argv);
    }
    else
    {
        // Sleep for 5 seconds
        std::this_thread::sleep_for(std::chrono::seconds(5));
        unsigned int rankUAVsNum = rank-1;
        for (int ii = 0; ii < 600 ; ii++)
        {
            ifCollisionSwapVelocities();
            CalcualteUAVsLocation(rankUAVsNum);
            sendRankUavInfo(sendBuffer, rankUAVsNum);
            MPI_Allgather(sendBuffer, numElements, MPI_DOUBLE, 
                        rcvbuffer, numElements, MPI_DOUBLE, MPI_COMM_WORLD);
            setUAVsArrayToRcvBuffer(rcvbuffer);
        }
    }
    return 0;
}