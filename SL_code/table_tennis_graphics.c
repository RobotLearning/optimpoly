#include "SL.h"
#include "SL_user.h"
#include "SL_man.h"

// openGL includes
#ifdef powerpc
#include <GLUT/glut.h>
#else
#include "GL/glut.h"
#endif

#include "SL_openGL.h"
#include "SL_userGraphics.h"

#include "math.h"
#include "table.h"

void display_ball_gun(void *b)
{
  double vars[14];
  double x,y,lx,ly;

  SL_quat varOrient;

  memcpy(&(vars[0]),b,14*sizeof(double));

  varOrient.q[_Q0_]=vars[7];
  varOrient.q[_Q1_]=vars[8];
  varOrient.q[_Q2_]=vars[9];
  varOrient.q[_Q3_]=vars[10];


  GLfloat  col[4] ={(float)1.0,(float)1.0,(float)0.0,(float)1.0};
  GLfloat  col2[4]={(float)0.2,(float)0.2,(float)0.2,(float)1.0};
  GLfloat black[4] = {(float)0.,(float)0.,(float)0.,(float)0.};
  GLfloat grey[4]  = {(float)0.7,(float)0.7,(float)0.7,(float)0.7};
  GLfloat red[4]   = {(float)0.,(float)0.9,(float)0.,(float)0.};

  GLfloat white[4] = {(float)1.0,(float)1.0,(float)1.0,(float)1.0};
  GLfloat forrestgreen[4] = {0, (float)0.13, (float)0.55, (float)0.13};


  /////////////////////////////////////////////////////////////////
  // Draws the Ball
  /////////////////////////////////////////////////////////////////
  glPushMatrix();
  glTranslated((GLdouble)vars[1],
	       (GLdouble)vars[2],
	       (GLdouble)vars[3]);  
  glColor4fv(col);
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, col);
  glutSolidSphere(ball_radius,8,8);
  glPopMatrix();

  glPushMatrix();
  glTranslated((GLdouble)vars[11],
	       (GLdouble)vars[12],
	       (GLdouble)vars[13]);  
  glColor4fv(col2);
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, col2);
  glutSolidSphere(ball_radius,8,8);
  glPopMatrix();

  /////////////////////////////////////////////////////////////////
  // Draws the Racket
  /////////////////////////////////////////////////////////////////  
  draw_rotated_disk(vars[4], vars[5],
		    vars[6], racket_radius, 0.012/2, black,
		    varOrient);
  draw_rotated_disk(vars[4], vars[5],
		    vars[6], racket_radius, -0.012/2, red,
		    varOrient);

  /////////////////////////////////////////////////////////////////
  // Draws the Table
  /////////////////////////////////////////////////////////////////
  x = table_center;
  y = dist_to_table - 0.5*table_length;
  lx = table_width;
  ly = table_length;
  draw_cube(x,y, floor_level-table_height-table_thickness/2.-0.001, 
	    lx, ly, table_thickness, white);
  draw_cube(x,y, floor_level-table_height-table_thickness/2., 
	    lx-0.02, ly-0.02, table_thickness, forrestgreen);
  draw_cube(x,y, floor_level-table_height - table_thickness/2. + 0.001, 
	    0.02, ly-0.02, table_thickness, white);	    	    
  draw_cube(x,y-0.005, floor_level-table_height+net_height/2+0.045/2., 
	    lx+2*net_overhang,net_thickness, net_height, white);	    

}

void add_table_tennis_graphics() {
  	addToUserGraphics("ball_gun", "Display Table, Ball and Racket", &(display_ball_gun), 14*sizeof(double));
}

