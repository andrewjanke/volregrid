/* regrid_view.c                                                             */
/*                                                                           */
/* quick widget to view paths                                                */
/*                                                                           */
/* Andrew Janke - rotor@cmr.uq.edu.au                                        */
/* Center for Magnetic Resonance                                             */
/* University of Queensland                                                  */
/*                                                                           */
/* Copyright Andrew Janke, The University of Queensland.                     */
/* Permission to use, copy, modify, and distribute this software and its     */
/* documentation for any purpose and without fee is hereby granted,          */
/* provided that the above copyright notice appear in all copies.  The       */
/* author and the University of Queensland make no representations about the */
/* suitability of this software for any purpose.  It is provided "as is"     */
/* without express or implied warranty.                                      */
/*                                                                           */
/* Wed Feb 12 13:15:50 EST 2003 - initial version                            */

#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <malloc.h>
#include <ParseArgv.h>
#include "arb_path_io.h"
#include "trackball.h"

#define TRUE   1
#define FALSE  0

/* forward declarations */
void     recalcModelView(void);

/* The number of our GLUT window */
int      winIdMain;

/* state variables */
int      W = 300, H = 300;
int      newModel = 1;
GLdouble bodyWidth = 3.0;
double   scalefactor = 1.0;
int      spinning = 0;
int      scaling = 0;
int      moving = 0;

int      beginx, beginy;

double   curquat[4];
double   lastquat[4];

/* coordinate stores */
int      n_coords = 0;
double  *x_coord;
double  *y_coord;
double  *z_coord;

double   max[3] = { -DBL_MAX, -DBL_MAX, -DBL_MAX };
double   min[3] = { DBL_MAX, DBL_MAX, DBL_MAX };

/* Argument variables */
int      verbose = FALSE;
int      points = TRUE;
int      lines = TRUE;
ArgvInfo argTable[] = {
   {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
    "Print out extra information"},
   {"-points", ARGV_CONSTANT, (char *)TRUE, (char *)&points, "Draw points"},
   {"-nopoints", ARGV_CONSTANT, (char *)FALSE, (char *)&points, "Don't draw points"},
   {"-lines", ARGV_CONSTANT, (char *)TRUE, (char *)&lines, "Draw lines"},
   {"-nolines", ARGV_CONSTANT, (char *)FALSE, (char *)&lines, "Don't draw lines"},

   {NULL, ARGV_HELP, NULL, NULL, ""},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

/* Callback function for reshaping the main window */
void mainReshape(int w, int h)
{
   glViewport(0, 0, w, h);
   W = w;
   H = h;
   }

void mainDisplay(void)
{
   int      c;

   if(newModel){
      recalcModelView();
      }
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   /* draw the bounding box */
   glLineWidth(1.0);
   glColor4f(0.0, 1.0, 0.0, 0.3);

   glBegin(GL_LINE_STRIP);
   glVertex3f(min[0], min[1], min[2]);
   glVertex3f(min[0], max[1], min[2]);
   glVertex3f(max[0], max[1], min[2]);
   glVertex3f(max[0], min[1], min[2]);
   glVertex3f(min[0], min[1], min[2]);

   glVertex3f(min[0], min[1], max[2]);
   glVertex3f(min[0], max[1], max[2]);
   glVertex3f(max[0], max[1], max[2]);
   glVertex3f(max[0], min[1], max[2]);
   glVertex3f(min[0], min[1], max[2]);
   glEnd();

   glBegin(GL_LINES);
   glVertex3f(max[0], min[1], min[2]);
   glVertex3f(max[0], min[1], max[2]);

   glVertex3f(max[0], max[1], min[2]);
   glVertex3f(max[0], max[1], max[2]);

   glVertex3f(min[0], max[1], min[2]);
   glVertex3f(min[0], max[1], max[2]);
   glEnd();

   /* Draw the points if wanted */
   if(points){
      glPointSize(2.0);
      glColor4f(0.8, 0.0, 0.9, 0.9);
      glBegin(GL_POINTS);
      for(c = 0; c < n_coords; c++){
         glVertex3f(x_coord[c], y_coord[c], z_coord[c]);
         }
      glEnd();
      }

   /* Draw the lines if wanted */
   if(lines){
      glLineWidth(0.1);
      glColor4f(0.0, 0.5, 0.9, 0.4);

      glBegin(GL_LINE_STRIP);
      for(c = 0; c < n_coords; c++){
         glVertex3f(x_coord[c], y_coord[c], z_coord[c]);
         }
      glEnd();
      }

   glFlush();
   glutSwapBuffers();
   }


void KeyPressed(unsigned char key, int x, int y)
{
   switch (key){
   case 'l':
      lines = !lines;
      break;
   case 'p':
      points = !points;
      break;

   case 'v':
      break;

   case 'd':
      break;

   case 27:                           /* escape */
      glutDestroyWindow(winIdMain);
      exit(TRUE);
      }
   glutPostRedisplay();
   }

void SpecialKey(int key, int x, int y)
{

   switch (key){
   case GLUT_KEY_RIGHT:
      break;

   case GLUT_KEY_LEFT:
      break;

   case GLUT_KEY_UP:
      break;

   case GLUT_KEY_DOWN:
      break;
      }

   glutPostRedisplay();
   }

int dump_points()
{
   int      c;

   for(c = 0; c < n_coords; c++)
      fprintf(stdout, "%03d - [%.06f][%.06f][%.06f]\n",
              c, x_coord[c], y_coord[c], z_coord[c]);
   return TRUE;
   }


void animate(void)
{
   add_quats(lastquat, curquat, curquat);
   newModel = 1;
   glutPostRedisplay();
   }

void motion(int x, int y)
{
   if(scaling){
      scalefactor = scalefactor * (1.0 + (((float)(beginy - y)) / H));
      beginx = x;
      beginy = y;
      newModel = 1;
      glutPostRedisplay();
      return;
      }
   if(moving){
      trackball(lastquat,
                (2.0 * beginx - W) / W,
                (H - 2.0 * beginy) / H, (2.0 * x - W) / W, (H - 2.0 * y) / H);
      beginx = x;
      beginy = y;
      spinning = 1;
      glutIdleFunc(animate);
      }
   }


void mouse(int button, int state, int x, int y)
{
   if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN){
      spinning = 0;
      glutIdleFunc(NULL);
      moving = 1;
      beginx = x;
      beginy = y;
      if(glutGetModifiers() & GLUT_ACTIVE_SHIFT){
         scaling = 1;
         }
      else{
         scaling = 0;
         }
      }
   if(button == GLUT_LEFT_BUTTON && state == GLUT_UP){
      moving = 0;
      }
   }

void recalcModelView(void)
{
   GLdouble m[4][4];

   glPopMatrix();
   glPushMatrix();
   build_rotmatrix(m, curquat);
   glMultMatrixd(&m[0][0]);
   if(scalefactor == 1.0){
      glDisable(GL_NORMALIZE);
      }
   else{
      glEnable(GL_NORMALIZE);
      }
   glScalef(scalefactor, scalefactor, scalefactor);
//  glTranslatef(-8, -8, -bodyWidth / 2);
   newModel = 0;
   }


void vis(int visible)
{
   if(visible == GLUT_VISIBLE){
      if(spinning)
         glutIdleFunc(animate);
      }
   else{
      if(spinning)
         glutIdleFunc(NULL);
      }
   }

int main(int argc, char **argv)
{
   int      c, i;
   Coord_list coord_buf;
   char    *coord_fn;
   int      buff_size = 2048;
   int      prev_size;

   /* Get arguments */
   if(ParseArgv(&argc, argv, argTable, 0)){
      (void)fprintf(stderr, "\nUsage: %s [options] <coords.cfg>\n", argv[0]);
      (void)fprintf(stderr, "       %s -help\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }
   coord_fn = argv[1];

   /* initialise the parser with the config file */
   fprintf(stderr, "Opening coord-file %s\n", coord_fn);
   if(!init_arb_path(coord_fn, NULL)){
      fprintf(stderr, "Failed to init arb_path, this isn't good\n");
      }

   /* suck in all the co-ordinates (brough-ha ha ha!) */
   coord_buf = get_some_arb_path_coords(buff_size);
   while (coord_buf->n_pts != 0){

      /* allocate some memory */
      prev_size = n_coords;
      n_coords += coord_buf->n_pts;
      x_coord = (double *)realloc(x_coord, n_coords * sizeof(double));
      y_coord = (double *)realloc(y_coord, n_coords * sizeof(double));
      z_coord = (double *)realloc(z_coord, n_coords * sizeof(double));

      /* add the current list */
      for(c = 0; c < coord_buf->n_pts; c++){
         x_coord[prev_size + c] = coord_buf->pts[c].coord[0];
         y_coord[prev_size + c] = coord_buf->pts[c].coord[1];
         z_coord[prev_size + c] = coord_buf->pts[c].coord[2];

         /* set min and max */
         for(i = 3; i--;){
            if(coord_buf->pts[c].coord[i] > max[i]){
               max[i] = coord_buf->pts[c].coord[i];
               }
            if(coord_buf->pts[c].coord[i] < min[i]){
               min[i] = coord_buf->pts[c].coord[i];
               }
            }

         }
      fprintf(stdout, "Got %d co-ords, realloced to %d points\n", coord_buf->n_pts,
              n_coords);

      /* get the next lot of co-ordinates */
      coord_buf = get_some_arb_path_coords(buff_size);
      }

   fprintf(stdout, "Got %d coords\n", n_coords);
   fprintf(stdout, "Min [%g:%g:%g] Max [%g:%g:%g]\n", min[0], min[1], min[2], max[0],
           max[1], max[2]);

   /* init the quat */
   trackball(curquat, 0.0, 0.0, 0.0, 0.0);

   /* GLUT stuff */
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
   glutCreateWindow("regrid viewer");
   glutDisplayFunc(mainDisplay);
   glutReshapeFunc(mainReshape);
   glutVisibilityFunc(vis);
   glutMouseFunc(mouse);
   glutMotionFunc(motion);
   glutKeyboardFunc(KeyPressed);
   glutSpecialFunc(SpecialKey);

   glMatrixMode(GL_PROJECTION);
   gluPerspective( /* field of view in degree */ 40.0,
                  /* aspect ratio */ 1.0,
                  /* Z near */ 1.0, /* Z far */ 40.0);
   glMatrixMode(GL_MODELVIEW);
   gluLookAt(0.0, 0.0, 30.0,           /* eye is at (0,0,30) */
             0.0, 0.0, 0.0,            /* center is at (0,0,0) */
             0.0, 1.0, 0.);            /* up is in positive Y direction */
   glPushMatrix();                     /* dummy push so we can pop on model
                                          recalc */

   glutMainLoop();

   return TRUE;
   }
