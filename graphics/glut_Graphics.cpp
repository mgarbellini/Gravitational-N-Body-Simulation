// F. Borando - francesco.borando@studenti.unimi.it
// M. Garbellini - matteo.garbellini@studenti.unimi.it
// Università degli Studi di Milano
// Dept. of Physics
// January 2019
//
// GRAVITATIONAL N-BODY SIMULATION
// Graphics for simulation visualization

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream> 			//Managing precision from float to string std::to_string
#include <iomanip>
#include <cmath>
#include <GL/glew.h>			//OpenGL library
#include <GL/freeglut.h>		//OpenGL library
#include <glm/glm.hpp>			//Mathematical library for matrices and transformations
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "glutBody.hpp"			//Custom made class for better handling of glutBodies

//* * * G L O B A L  D E C L A R A T I O N * * * * * *

std::vector<glutBody> glut_bodies;		//All-bodies-all-iterations vector
std::ifstream input;				//Input file definition
int N;						//Number of bodies
float timestep;					//Computational timestep
int iter=0;					//Number of animation iterations (i.e. frames)
int iterations=0;				//Number of computational iterations*bodies
float rotateX = 0;				//Value of rotation around X-axis
float rotateZ = 0;				//Value of rotation around Z-axis
float scale = 1 ;				//Scaling factor (default = 1)
const float scale_value = 0.3;			//Scaling factor increment
float max=0;					//Max between x,y,z coordinates (calculated for first 3 iter)
float point_size = 1;				//Point size (variable) to allow better viewing
const float pnt_size_var = 0.2; 		//Point size variation value
int FPS;					//Value of Frame-Per-Second (default 30fps)
const float FPS_var = 15;			//Value of Frame-Per-Second variation (suggested 15)
int animating = 0;				//Bool value for animation on and off
int reverse = 0;				//Bool value for reverse playback
int is_fullscreen = 0;				//Bool value for full screen mode
int key_informations = 1;			//Bool value for key_informations on screen
int key_acknowledgment = 1;			//Bool value for acknowledgment on screen

//TEXT DECLARATIONS
//simulations parameters
std::string	text_title,text_iter, text_iter_elapsed,
text_time, text_time_elapsed, text_bodies, text_bodies_number;
//key_informations
std::string 	key_F1, key_F2, key_F3, key_F4,
key_F5, key_F6, key_F11, key_F12,
key_UP_DOWN, key_LEFT_RIGHT, key_full,
key_EXIT, key_SPACE, key_info, key_ack;

//acknowledgments
std::string 	ack_name, ack_affil, ack_dept;

//* * * End of G L O B A L  D E C L A R A T I O N * * *



//* * * F U N C T I O N S * * *

//inline function needed to adjust iter count. (Whenever a key is pressed glutPostRedisplay is called and iter
//is increased/decreases by glutDisplay function, wrongly updating iter count)
inline void adj_iter(int& iter)
{
                        if (!reverse)
                        {
                                                iter--;
                        }
                        else
                        {
                                                iter++;
                        }
}

//functions for text drawings
void drawText(const char *text, const int& length, const int& x, const int& y, const float& text_scale )
{
                        glMatrixMode(GL_PROJECTION); 						//Switch current matrix to PROJECTION
                        double matrix[16]; 							//Define double matrix (dim = 16)
                        glGetDoublev(GL_PROJECTION_MATRIX, matrix); 				//Store PROJECTION matrix value into matrix (this is used for later reset)
                        glLoadIdentity(); 							//Reset PROJECTION matrix to identity matrix
                        glOrtho(0, 1500, 0, 1500, -5, 5); 					//Define orthographic perspective (left_bottom, left_up, right_bottom, right_up, v_near, v_far)


                        glMatrixMode(GL_MODELVIEW); 						//Switch current matrix to MODELVIEW
                        glLoadIdentity(); 							//Reset MODELVIEW matrix to idetity matrix

                        glPushMatrix(); 							//Push matrix
                        glTranslatef(x, y, 0);							//Translate text position to (x,y) coordinates
                        glScalef(text_scale,text_scale,1);					//Scale text dimension
                        for(int i=0; i<length; i++)
                        {
                                                glutStrokeCharacter(GLUT_STROKE_ROMAN, (int)text[i]);		//GLUT_STROKE_ROMAN allows scaling of characters
                        }
                        glPopMatrix(); 								//Pop matrix
                        glMatrixMode(GL_PROJECTION); 						//Switch current matrix to PROJECTION
                        glLoadMatrixd(matrix); 							//Load previously saved matrix to PROJECTION matrix
                        glMatrixMode(GL_MODELVIEW); 						//Switch current matrix to MODELVIEW
}

                        //NOTE: this functions does not do anything particular, except assigning text value to globally defined std::string,
                        //convenient for developers and cleanliness of code
                        void assign_text_values()
                        {
                                                //text for simulation parameters and titles
                                                text_title 		= "GRAVITATIONAL N-BODY SIMULATION";
                                                text_bodies  		= "# of Bodies:";
                                                text_time 		= "Time elapsed: ";
                                                text_iter 		= "Computations: ";

                                                //text for key-informations
                                                key_F1 		= "F1: zoom out ";
                                                key_F2 		= "F2: reset ";
                                                key_F3 		= "F3: zoom in";
                                                key_F4 		= "F4: increase FPS";
                                                key_F5 		= "F5: reset FPS";
                                                key_F6 		= "F6: decrease FPS";
                                                key_UP_DOWN 	= "Up/Down: rotate around x-axis";
                                                key_LEFT_RIGHT 	= "Left/Right : rotate around z-axis";
                                                key_SPACE 	= "Spacebar: play/pause anim.";
                                                key_EXIT 	= "Esc: exit simulation";
                                                key_F11 	= "F11: reverse anim.";
                                                key_F12 	= "F12: restart anim. ";
                                                key_info	= "(i): hide/show informations ";
                                                key_ack 	= "(a): display acknowledgments";
                                                key_full 	= "(f): full screen";

                                                //text for acknowlegment
                                                ack_name 	= "F. BORANDO, M. GARBELLINI";
                                                ack_affil 	= "Universita' degli Studi di Milano";
                                                ack_dept 	= "Department of Physics";
                        }

                        void draw_key_informations_text(const float& text_scale)
                        {
                                                int info_x_coord = 40;
                                                int info_y_coord = 1000;
                                                int info_y_offwidth = -40;

                                                drawText(key_F1.data(), key_F1.size(), info_x_coord, info_y_coord, text_scale);
                                                drawText(key_F2.data(), key_F2.size(), info_x_coord, info_y_coord + info_y_offwidth, text_scale);
                                                drawText(key_F3.data(), key_F3.size(), info_x_coord, info_y_coord + 2*info_y_offwidth, text_scale);
                                                drawText(key_F4.data(), key_F4.size(), info_x_coord, info_y_coord + 3*info_y_offwidth, text_scale);
                                                drawText(key_F5.data(), key_F5.size(), info_x_coord, info_y_coord + 4*info_y_offwidth, text_scale);
                                                drawText(key_F6.data(), key_F6.size(), info_x_coord, info_y_coord + 5*info_y_offwidth, text_scale);
                                                drawText(key_F11.data(), key_F11.size(), info_x_coord, info_y_coord + 6*info_y_offwidth, text_scale);
                                                drawText(key_F12.data(), key_F12.size(), info_x_coord, info_y_coord + 7*info_y_offwidth, text_scale);
                                                drawText(key_UP_DOWN.data(), key_UP_DOWN.size(), info_x_coord, info_y_coord + 8*info_y_offwidth, text_scale);
                                                drawText(key_LEFT_RIGHT.data(), key_LEFT_RIGHT.size(), info_x_coord, info_y_coord + 9*info_y_offwidth, text_scale);
                                                drawText(key_SPACE.data(), key_SPACE.size(), info_x_coord, info_y_coord + 10*info_y_offwidth, text_scale);
                                                drawText(key_EXIT.data(), key_EXIT.size(), info_x_coord, info_y_coord + 11*info_y_offwidth, text_scale);
                                                drawText(key_full.data(), key_full.size(), info_x_coord, info_y_coord + 12*info_y_offwidth, text_scale);
                                                drawText(key_ack.data(), key_ack.size(), info_x_coord, info_y_coord + 13*info_y_offwidth, text_scale);
                                                drawText(key_info.data(), key_info.size(), info_x_coord, info_y_coord + 14*info_y_offwidth, text_scale);
                        }

                        void draw_simulation_text(const float& text_scale)
                        {
                                                //float to string with precision
                                                std::stringstream stream;
                                                stream << std::fixed << std::setprecision(1) << iter*timestep;
                                                text_time_elapsed = stream.str();

                                                text_iter_elapsed = std::to_string(iter);

                                                //Text position (upper left corner)
                                                const int text_x_coord = 40;
                                                const int text_y_coord = 1400;
                                                const int text_x_offwidth = 20;
                                                const int text_y_offwidth = -40;

                                                drawText(text_title.data(), text_title.size(), text_x_coord, text_y_coord,text_scale );
                                                drawText(text_bodies.data(), text_bodies.size(), text_x_coord, text_y_coord + text_y_offwidth, text_scale);
                                                drawText(text_bodies_number.data(), text_bodies_number.size(), text_x_coord + text_bodies.size()*text_x_offwidth, text_y_coord + text_y_offwidth, text_scale);
                                                drawText(text_time.data(), text_time.size(), text_x_coord , text_y_coord + 2*text_y_offwidth, text_scale);
                                                drawText(text_time_elapsed.data(), text_time_elapsed.size(), text_x_coord + text_time.size()*text_x_offwidth, text_y_coord + 2*text_y_offwidth, text_scale);
                                                drawText(text_iter.data(), text_iter.size(), text_x_coord, text_y_coord + 3* text_y_offwidth, text_scale);
                                                drawText(text_iter_elapsed.data(), text_iter_elapsed.size(), text_x_coord + text_iter.size()* text_x_offwidth, text_y_coord + 3*text_y_offwidth, text_scale);

                        }

                        void draw_acknowledgment_text(const float& text_scale)
                        {
                                                //Text position
                                                int text_x_coord = 1150; //40
                                                int text_y_coord = 1400;//200
                                                int text_x_offwidth = 30;
                                                int text_y_offwidth = -40;

                                                drawText(ack_name.data(), ack_name.size(), text_x_coord, text_y_coord, text_scale);
                                                drawText(ack_affil.data(), ack_affil.size(), text_x_coord, text_y_coord + text_y_offwidth, text_scale);
                                                drawText(ack_dept.data(), ack_dept.size(), text_x_coord, text_y_coord + 2*text_y_offwidth, text_scale);
                        }

                        //Function for fullscreen windows. Uses native glutReshapeWindow
                        void fullscreen()
                        {

                                                if(is_fullscreen)
                                                {
                                                                        glutReshapeWindow(1500, 1500);
                                                                        glutPositionWindow(0, 0);
                                                }
                                                else
                                                {
                                                                        glutFullScreen();
                                                }

                                                is_fullscreen = !is_fullscreen;

                                                adj_iter(iter);
                        }

                        //Framerate function using OpenGL native glutTimerFunc.
                        //this function is called recursively within glutTimerFunc. OpenGL
                        //definition requires that an int value is passed to idle(), although not used
                        void idle(int)
                        {
                                                if(animating)
                                                {
                                                                        glutPostRedisplay();
                                                                        glutTimerFunc(1000/FPS, idle, 0);
                                                }
                        }

                        //Function that handles playing and pausing animation and exit simulation command (Esc), restart animation, key_info and acknowledgment_info
                        void glutHandleKeyPress(const unsigned char key, const int x, const int y)
                        {
                                                if (key == ' ' )
                                                {
                                                                        if (animating)
                                                                        {
                                                                                                animating = 0;
                                                                        }
                                                                        else
                                                                        {
                                                                                                animating = 1;
                                                                                                glutTimerFunc(1000/FPS, idle, 0);
                                                                        }
                                                }

                                                if (key == 'f')
                                                {
                                                                        fullscreen();
                                                }

                                                if (key == 27)
                                                {
                                                                        exit(0);
                                                }

                                                if (key == '0') //restart animation
                                                {
                                                                        if (reverse)
                                                                        {
                                                                                                iter = iterations/N - 1;
                                                                        }
                                                                        else
                                                                        {
                                                                                                iter = 0;
                                                                        }
                                                }

                                                if (key == 'i')
                                                {
                                                                        if (key_informations)
                                                                        {
                                                                                                key_informations = 0;
                                                                        }
                                                                        else
                                                                        {
                                                                                                key_informations = 1;
                                                                        }

                                                                        adj_iter(iter);
                                                                        glutPostRedisplay();
                                                }

                                                if (key == 'a')
                                                {
                                                                        if (key_acknowledgment)
                                                                        {
                                                                                                key_acknowledgment = 0;
                                                                        }
                                                                        else
                                                                        {
                                                                                                key_acknowledgment = 1;
                                                                        }
                                                                        adj_iter(iter);
                                                                        glutPostRedisplay();
                                                }
                        }

                        //Special keys (UP, DOWN, LEFT, RIGHT, F1-6) to control coordinates rotatation,
                        //zooming in and out, speed up slow down. The definition comes directly from OpenGL documentation
                        void glutSpecialKeys(const int key, const int x, const int y)
                        {
                                                switch (key)
                                                {
                                                                        case GLUT_KEY_LEFT:
                                                                        rotateZ -= 15;
                                                                        adj_iter(iter);
                                                                        break;

                                                                        case GLUT_KEY_RIGHT:
                                                                        rotateZ += 15;
                                                                        adj_iter(iter);
                                                                        break;

                                                                        case GLUT_KEY_DOWN:
                                                                        rotateX += 15;
                                                                        adj_iter(iter);
                                                                        break;

                                                                        case GLUT_KEY_UP:
                                                                        rotateX -= 15;
                                                                        adj_iter(iter);
                                                                        break;

                                                                        case GLUT_KEY_F1:
                                                                        scale -=  scale_value;
                                                                        point_size -= pnt_size_var;
                                                                        if (scale <= 0)
                                                                        {
                                                                                                scale += scale_value;
                                                                        }
                                                                        if (point_size <= 1)
                                                                        {
                                                                                                point_size = 1;
                                                                        }

                                                                        adj_iter(iter);

                                                                        break;

                                                                        case GLUT_KEY_F2:
                                                                        scale = 1;
                                                                        point_size = 1;
                                                                        adj_iter(iter);

                                                                        break;

                                                                        case GLUT_KEY_F3:
                                                                        scale += scale_value;
                                                                        point_size += pnt_size_var;
                                                                        adj_iter(iter);

                                                                        break;

                                                                        case GLUT_KEY_F4:
                                                                        FPS -= FPS_var;
                                                                        if(FPS<=0)
                                                                        {
                                                                                                FPS = 15;
                                                                        }
                                                                        adj_iter(iter);

                                                                        break;

                                                                        case GLUT_KEY_F5:
                                                                        FPS = 30;
                                                                        adj_iter(iter);

                                                                        break;

                                                                        case GLUT_KEY_F6:
                                                                        FPS += FPS_var;
                                                                        adj_iter(iter);

                                                                        break;

                                                                        case GLUT_KEY_F11:
                                                                        if (reverse  == 0)
                                                                        {
                                                                                                reverse = 1;
                                                                                                iter = (int)iterations/N - 1;
                                                                        }
                                                                        else
                                                                        {
                                                                                                reverse = 0;
                                                                                                iter = 0;
                                                                        }

                                                                        case GLUT_KEY_F12:
                                                                        if(reverse)
                                                                        {
                                                                                                iter = iterations/N - 1;
                                                                        }
                                                                        else
                                                                        {
                                                                                                iter = 0;
                                                                        }
                                                                        break;

                                                                        default:
                                                                        std::cerr<<"Key '"<< key <<"'not assigned to any action!"<<std::endl;
                                                }

                                                glutPostRedisplay();
                        }

                        //Function finds max value for normalization (OpenGL requires coordinates in the range[0,1])
                        //to avoid 'big' normalization, max value is calculated only for first iteration
                        void find_max(const float x, const float y, const float z, float &max)
                        {

                                                if (x > max)
                                                {
                                                                        max = x;
                                                }
                                                else if ((-x) > max)
                                                {
                                                                        max = -x;
                                                }
                                                else if (y > max)
                                                {
                                                                        max = y;
                                                }
                                                else if ((-y) > max)
                                                {
                                                                        max = -y;
                                                }
                                                else if (z > max)
                                                {
                                                                        max = z;
                                                }
                                                else if ((-z) > max)
                                                {
                                                                        max = -z;
                                                }
                        }

                        //Main function for OpenGL renderization. glutDisplay is called by glutDisplayFunc(), and it
                        //must contain all the objects that need to be drawn
                        void glutDisplay(void)
                        {

                                                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);			//Clears the window for drawing
                                                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);			//Defines blending and transparency options (currently not used)
                                                glEnable(GL_BLEND);							//Enables blending and transparency
                                                glEnable(GL_POINT_SMOOTH);						//Avoids aliasing problems with many-points. Smoother animation


                                                glMatrixMode(GL_PROJECTION);						//Switch to PROJECTION mode: for 3D rotations
                                                glLoadIdentity();							//Loads Identity matrix at current position.
                                                //futures rotations will use this as a 'starting point'
                                                glRotatef(rotateX,1,0,0);						//Rotates around x-axis of rotateX degrees
                                                glRotatef(rotateZ,0,1,0);						//Rotates around y-axis of rotateY degrees


                                                //Draw text
                                                text_bodies_number = std::to_string(N);
                                                draw_simulation_text(0.25);
                                                if(key_informations)
                                                {
                                                                        draw_key_informations_text(0.20);
                                                }
                                                if(key_acknowledgment)
                                                {
                                                                        draw_acknowledgment_text(0.15);
                                                }

                                                //Automatic reset of iter counter to start over (indefinitely) with the simulation <-----LOOP CONDITION
                                                if(reverse == 0 && iter == (int)iterations/N)
                                                {
                                                                        iter = 0;
                                                }
                                                if(reverse == 1 && iter == 0)
                                                {
                                                                        iter = (int)iterations/N - 1;
                                                }

                                                if(reverse == 0)
                                                {
                                                                        for(int i=iter*N; i<N*(iter+1); i++)
                                                                        {
                                                                                                glut_bodies[i].glutScaleCoord(scale);				//Scales points according to 'scale' parameter. Used for zooming in and out
                                                                                                glut_bodies[i].glutDrawPoint(point_size);			//glutDrawPoint(float point size)
                                                                                                glut_bodies[i].glutScaleCoord((float) 1/scale);
                                                                        }

                                                                        glutSwapBuffers();						//Swap buffers to current one
                                                                        iter++;									//Incrementation of iter counter
                                                }
                                                else if (reverse==1)
                                                {
                                                                        //iter = iterations - 1;
                                                                        for(int i=iter*N; i<N*(iter+1); i++)
                                                                        {
                                                                                                //glut_bodies[i].glutDrawSphere();
                                                                                                glut_bodies[i].glutScaleCoord(scale);				//Scales points according to 'scale' parameter. Used for zooming in and out
                                                                                                glut_bodies[i].glutDrawPoint(point_size);				//glutDrawPoint(float point size, float alpha)
                                                                                                glut_bodies[i].glutScaleCoord((float) 1/scale);
                                                                        }

                                                                        glutSwapBuffers();							//Swap buffers to current one
                                                                        iter = iter - 1;
                                                }
                        }

                        //Simple routine to load vector glutbodies from file argv[1]
                        void glutLoadVector(char* argv[])
                        {
                                                std::ifstream input;
                                                float x_coord_read, y_coord_read, z_coord_read;
                                                glutBody body_to_push;

                                                //Input routine
                                                input.open(argv[1]);
                                                input>>N>>timestep;
                                                while(input.good())
                                                {
                                                                        input 	>> x_coord_read
                                                                        >> y_coord_read
                                                                        >> z_coord_read;

                                                                        if(iterations < 2*N)
                                                                        {
                                                                                                find_max(x_coord_read, y_coord_read, z_coord_read, max);
                                                                        }
                                                                        body_to_push.set_glutPosition(	x_coord_read, y_coord_read, z_coord_read);

                                                                        glut_bodies.push_back(body_to_push);				//push back to glut_bodies (all bodies) vector
                                                                        iterations++;							//iteration counter (counts body*iter)
                                                }

                                                for(int i=0; i<iterations; i++)
                                                {
                                                                        glut_bodies[i].glutNormalize(max);
                                                }

                                                input.close();
                        }
                        //* * * End of F U N C T I O N S  * * *



                        int main(int argc, char* argv[])
                        {

                                                //argv[] error outputs
                                                if(!argv[1])
                                                {
                                                                        std::cerr << "Please filename.txt as argv[1]!" << std::endl;
                                                                        std::exit (EXIT_FAILURE);
                                                }

                                                if(!argv[2])
                                                {
                                                                        std::cout<<"Default FPS selected: 30FPS"<<std::endl;
                                                                        FPS = 30;								//Default framerate 30fps
                                                }
                                                else
                                                {
                                                                        FPS = std::atoi(argv[2]);
                                                }


                                                //Load vector
                                                glutLoadVector(argv);
                                                //Assign text values
                                                assign_text_values();

                                                //OpenGL routines for setting up window
                                                glutInit(&argc, argv);							//Initialize glut
                                                glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);		//Inizialize display mode: Color:Red-Green-Blue-Alpha, Color Depth, Double-Buffering
                                                glutInitWindowSize(1500,1500);						//Define window size
                                                glutInitWindowPosition(0,0);						//Define window initial position
                                                glutCreateWindow("Gravitational N-Body Simulation");			//Inizialize window with title
                                                //glutSetColor(1.0i, 1.0f, 1.0f, 1.0f);					//Set color for drawings


                                                glViewport(0,0,1500,1500);						//Set viewport for 3D graphics
                                                glm::mat4 projection = glm::perspective(40.0f, 1.0f, 20.0f, 1000.0f); 	//Define perspective matrix. glm:perspective(FOV, W/H, nearest_point, furthest_point)

                                                glMatrixMode(GL_PROJECTION);						//Matrix mode: projection
                                                glLoadIdentity();							//Load Identity matrix for further transformations
                                                glLoadMatrixf(glm::value_ptr(projection));				//New routine (OpenGL 3.0) for perspective matrix loading

                                                glutDisplayFunc(glutDisplay);						//Rendering function
                                                glutKeyboardFunc(glutHandleKeyPress);					//Handles input from keyboard
                                                glutSpecialFunc(glutSpecialKeys);					//Handles input from 'special_keys'
                                                glutMainLoop();								//OpenGL main loop

                                                return 0;
                        }
