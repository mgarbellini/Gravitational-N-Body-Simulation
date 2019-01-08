// F. Borando - francesco.borando@studenti.unimi.it
// M. Garbellini - matteo.garbellini@studenti.unimi.it
// Università degli Studi di Milano
// Dept. of Physics
// January 2019
//
// GRAVITATIONAL N-BODY SIMULATION
// glutBody.hpp: custom class for better handling of glut bodies




#ifndef GLUT_BODY_FOR_GRAVITATIONAL_SIMULATION
#define GLUT_BODY_FOR_GRAVITATIONAL_SIMULATION

#include <iostream>
#include <cmath>

//OpenGL libraries (freeglut)
#include <GL/glew.h>
#include <GL/freeglut.h>

class glutBody
{
	public:
		
		glutBody()
		{
			m_x = 0;
			m_y = 0;
		}

		glutBody(GLfloat x, GLfloat y)
		{
			m_x = x;
			m_y = y;
		}

		~glutBody()
		{
		}

		GLfloat get_glutPosition_x(){return m_x;}
		GLfloat get_glutPosition_y(){return m_y;}
		GLfloat get_glutPosition_z(){return m_z;}

		void set_glutPosition(GLfloat x, GLfloat y, GLfloat z)
		{
			m_x = x;
			m_y = y;
			m_z = z;
		}

		void glutScaleCoord(float scale)
		{
			m_x = m_x * scale;
			m_y = m_y * scale;
			m_z = m_z * scale;
		}

		void glutDrawPoint(float point_size)
		{
			glPointSize(point_size);
			glBegin(GL_POINTS);
            	glVertex3f(m_x, m_y, m_z); 
        	glEnd();
		}

		void glutNormalize(GLfloat max)
		{
			m_x = m_x / max;
			m_y = m_y / max;
			m_z = m_z /max;
		}

		void glutClear()
		{
			m_x = 0;
			m_y = 0;
		}

	private:

		GLfloat m_x;
		GLfloat m_y;
		GLfloat m_z;

};

#endif //GLUT_BODY_FOR_GRAVITATIONAL_SIMULATION
