//---------------------------------------------------------------------------
//
// voronoi rendering
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
// (C) 2011 by Argonne National Laboratory.
// See COPYRIGHT in top-level directory.
//
//--------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "voronoi.h"
#include "ser_io.hpp"
#include <math.h>

#if defined(MAC_OSX)
#include <GLUT/glut.h> 
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h> 
#include <GL/gl.h>
#endif

using namespace std;

// 3d point or vector
struct vec3d {
  float x, y, z;
};

// 2d point or vector
struct vec2d {
  float x, y;
};

// color
struct rgba {
  float r, g, b, a;
};

// mouse button state
int press_x, press_y; 
int release_x, release_y; 

// rotate
vec2d rot = {0.0, 0.0};
float rot_rate = 0.2;

// scale
float scale = 1.0; 
float scale_rate = 0.01;
vec2d aspect; // scaling due to window aspect ratio

// near clip plane
float near;

// window size
// vec2d win_size = {1024, 512};
vec2d win_size = {1024, 1024};
// vec2d win_size = {512, 512};

// previous window size
vec2d old_win_size;

// translate
vec2d trans = {0.0, 0.0};
float trans_rate = 0.01;

// transform mode
int xform_mode = 0; 
bool block_mode = false;
#define XFORM_NONE    0 
#define XFORM_ROTATE  1
#define XFORM_SCALE   2 
#define XFORM_TRANS   3

// rendering mode
bool draw_fancy = false;

// volume filtering
float min_vol = 0.0;
float vol_step = 0.1;

vec3d sizes; // individual data sizes in each dimension
float size; // one overall number for data size, max of individual sizes

// voronoi sites
vector<vec3d> sites;
vec3d site_min;
vec3d site_max;
vec3d site_center;

// voronoi vertices, faces, cells
vector<vec3d> verts;
vector<int> num_face_verts;

// voronoi blocks
vblock_t **vblocks;
int nblocks;

// general prupose quadrics
GLUquadricObj *q;

// point sprite texture
static GLubyte sprite_intensity[5][5] = {
  {  50,    50,   50,   50,  50,  },
  {  50,   100,  100,  100,  50,  },
  {  50,   100,  255,  100,  50,  },
  {  50,   100,  100,  100,  50,  },
  {  50,    50,   50,   50,  50,  },
};
static GLubyte sprite_rgba[5][5][4];
static GLuint tex;

// function prototypes
void display();
void init_display();
void draw_cube(float *mins, float *maxs, float r, float g, float b) ;
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void key(unsigned char key, int x, int y);
void timer(int val);
void draw_sphere(rgba &color, vec3d &pos, float rad);
void draw_spheres(vector<vec3d> &sites, float rad);
void draw_sprites(vector<vec3d> &sites, float size);
void draw_axes();
void reshape(int w, int h);
void init_model();
void init_viewport(bool reset);
void headlight();
void ComputeNormal(vec3d *verts, int num_verts, vec3d &normal);
void filter_volume(float min_vol);

//--------------------------------------------------------------------------

int main(int argc, char** argv) {

  int n, m;

  if (argc < 3) {
    fprintf(stderr, "Usage: draw <filename> <swap (0 or 1)>"
	    " [min. volume (optional)]\n");
    exit(0);
  }

  int swap_bytes = atoi(argv[2]);

  if (argc > 3)
    min_vol = atof(argv[3]);

  // read the file
  SER_IO *io = new SER_IO(swap_bytes); // io object
  nblocks = io->ReadAllBlocks(argv[1], vblocks, false);

  // package rendering data
  for (int i = 0; i < nblocks; i++) { // blocks

    n = 0;
    for (int j = 0; j < vblocks[i]->num_orig_particles; j++) {

	vec3d s;
	s.x = vblocks[i]->sites[n];
	s.y = vblocks[i]->sites[n + 1];
	s.z = vblocks[i]->sites[n + 2];
	n += 3;
	sites.push_back(s);

    }

    n = 0;
    m = 0;
    for (int j = 0; j < vblocks[i]->num_complete_cells; j++) { // cells

      for (int k = 0; k < vblocks[i]->num_cell_faces[j]; k++) { // faces

	if (vblocks[i]->vols[j] >= min_vol)
	  num_face_verts.push_back(vblocks[i]->num_face_verts[n]);

	for (int l = 0; l < vblocks[i]->num_face_verts[n]; l++) { // vertices

	  int v = vblocks[i]->face_verts[m];
	  vec3d s;
	  s.x = vblocks[i]->save_verts[3 * v];
	  s.y = vblocks[i]->save_verts[3 * v + 1];
	  s.z = vblocks[i]->save_verts[3 * v + 2];
	  m++;
	  if (vblocks[i]->vols[j] >= min_vol)
	    verts.push_back(s);

	} // vertices

	n++;

      } // faces

    } // cells

  } // blocks

  // start glut
  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH); 
  glutInitWindowSize(win_size.x, win_size.y); 
  glutCreateWindow("Voronoi"); 
  glutDisplayFunc(display); 
  glutTimerFunc(10, timer, 0); 
  glutMouseFunc(mouse); 
  glutMotionFunc(motion);
  glutKeyboardFunc(key); 
  glutReshapeFunc(reshape);
  glutMainLoop(); 

}
//--------------------------------------------------------------------------
//
// rendering
//
void display() {

  static bool first = true;
  int n;

  if (first)
    init_display();
  first = false;

  // set the headlight
  headlight();

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity(); 
  gluPerspective(60.0, 1.0, near, 100.0); 

  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 
  gluLookAt(0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0); 

  // mouse interactions: pan, rotate, zoom
  glTranslatef(trans.x, trans.y, 0.0);
  glRotatef(rot.x, 0.0, 1.0, 0.0); 
  glRotatef(rot.y, 1.0, 0.0, 0.0); 
  glScalef(scale, scale, scale);

  // center the data in the window
  glTranslatef(-site_center.x, -site_center.y, -site_center.z);

  glEnable(GL_COLOR_MATERIAL);

  // draw axes
  draw_axes();

  // draw cell edges
  glDisable(GL_LIGHTING);
  glColor4f(0.5, 0.5, 0.5, 1.0);
  glLineWidth(2);
  n = 0;
  for (int i = 0; i < (int)num_face_verts.size(); i++) {
    glBegin(GL_LINE_STRIP);
    int n0 = n; // index of first vertex in this face
    for (int j = 0; j < num_face_verts[i]; j++) {
      glVertex3f(verts[n].x, verts[n].y, verts[n].z);
      n++;
    }
    glVertex3f(verts[n0].x, verts[n0].y, verts[n0].z); // repeat first vertex
    glEnd();
  }

  // draw block bounds
  if (block_mode) {
    for (int i = 0; i < nblocks; i++)
      draw_cube(vblocks[i]->mins, vblocks[i]->maxs, 1.0, 0.0, 1.0);
  }

  // draw sites
  if (draw_fancy) {
    glDisable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    GLfloat amb_mat[] = {0.6, 0.6, 0.6, 1.0};
    GLfloat spec_mat[] = {1.0, 1.0, 1.0, 1.0};
    GLfloat shine[] = {1}; // 0 - 128, 0 = shiny, 128 = dull
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec_mat);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, amb_mat);
    draw_spheres(sites, 0.04);
    //     draw_sprites(sites, 10.0);
  }
  else {
    glColor3f(1.0, 1.0, 1.0);
    glEnable(GL_POINT_SMOOTH);
    glPointSize(1.0);
    glBegin(GL_POINTS);
    for (int i = 0; i < (int)sites.size(); i++)
      glVertex3f(sites[i].x, sites[i].y, sites[i].z);
    glEnd();
    glDisable(GL_COLOR_MATERIAL);
  }

  // draw cell faces
  if (draw_fancy) {
    GLfloat front_mat[] = {0.3, 0.35, 0.5, 1.0};
    GLfloat back_mat[] = {1.0, 0.3, 0.3, 1.0}; // test if back faces ever seen
    GLfloat spec_mat[] = {0.3, 0.3, 0.3, 1.0};
    GLfloat shine[] = {128}; // 0 - 128, 0 = shiny, 128 = dull
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec_mat);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, front_mat);
    glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, back_mat);
    n = 0;
    for (int i = 0; i < (int)num_face_verts.size(); i++) {
      // flat shading, one normal per face
      vec3d normal;
      ComputeNormal(&verts[n], num_face_verts[i], normal);
      glBegin(GL_POLYGON);
      for (int j = 0; j < num_face_verts[i]; j++) {
	glNormal3f(normal.x, normal.y, normal.z);
	glVertex3f(verts[n].x, verts[n].y, verts[n].z);
	n++;
      }
      glEnd();
    }
  }

  glutSwapBuffers();

}
//--------------------------------------------------------------------------
//
// first time drawing initialization
//
void init_display() {

  // extents
  for (int i = 0; i < (int)sites.size(); i++) {
    if (i == 0) {
      site_min.x = sites[i].x;
      site_min.y = sites[i].y;
      site_min.z = sites[i].z;
      site_max.x = sites[i].x;
      site_max.y = sites[i].y;
      site_max.z = sites[i].z;
    }
    if (sites[i].x < site_min.x)
      site_min.x = sites[i].x;
    if (sites[i].y < site_min.y)
      site_min.y = sites[i].y;
    if (sites[i].z < site_min.z)
      site_min.z = sites[i].z;
    if (sites[i].x > site_max.x)
      site_max.x = sites[i].x;
    if (sites[i].y > site_max.y)
      site_max.y = sites[i].y;
    if (sites[i].z > site_max.z)
      site_max.z = sites[i].z;
  }
  site_center.x = (site_min.x + site_max.x) / 2.0;
  site_center.y = (site_min.x + site_max.y) / 2.0;
  site_center.z = (site_min.x + site_max.z) / 2.0;
  sizes.x = site_max.x - site_min.x;
  sizes.y = site_max.y - site_min.y;
  sizes.z = site_max.z - site_min.z;
  size = sizes.x;
  if (sizes.y > size)
    size = sizes.y;
  if (sizes.z > size)
    size = sizes.z;

  init_model();
  init_viewport(true);

  // background
  glClearColor(0.0, 0.0, 0.0, 1.0); 

  // gl state
//   glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_DEPTH_TEST);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_LIGHT2);
  glEnable(GL_NORMALIZE);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glShadeModel(GL_SMOOTH);

  // initialize headlight
  headlight();

  // general purpose quadrics
  q = gluNewQuadric();

  // point sprite texture
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      sprite_rgba[i][j][0] = sprite_intensity[i][j];
      sprite_rgba[i][j][1] = sprite_intensity[i][j];
      sprite_rgba[i][j][2] = sprite_intensity[i][j];
      sprite_rgba[i][j][3] = sprite_intensity[i][j];
    }
  }
  glGenTextures(1, &tex);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 5, 5, 0, GL_RGBA, 
	       GL_UNSIGNED_BYTE, sprite_rgba);

  // near clip plane
  near = 0.1;

}
//--------------------------------------------------------------------------
//
// set a headlight
//
void headlight() {

  GLfloat light_ambient[4] = {0.1, 0.1, 0.1, 1.0};  
  GLfloat light_diffuse[4] = {0.1, 0.1, 0.1, 1.0};  
  GLfloat light_specular[4] = {0.7, 0.7, 0.7, 1.0};

  glPushMatrix();

  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 
  glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glEnable(GL_LIGHT0);

  glPopMatrix();

}
//--------------------------------------------------------------------------
//
// cube, useful for block boundaries, etc.
//
void draw_cube(float *mins, float *maxs, float r, float g, float b) {

  glPushMatrix();
  glTranslatef((mins[0] + maxs[0]) / 2.0, (mins[1] + maxs[1]) / 2.0, 
	       (mins[2] + maxs[2]) / 2.0); 
  glScalef(maxs[0] - mins[0], maxs[1] - mins[1], maxs[2] - mins[2]);
  glColor3f(r, g, b); 
  glutWireCube(1.0);
  glPopMatrix();

}
//--------------------------------------------------------------------------
// 
// sphere for rendering voronoi sites (particles)
//
void draw_sphere(rgba &color, vec3d &pos, float rad) {

  glColor3f(color.r, color.g, color.b); 
  glPushMatrix();
  glTranslatef(pos.x, pos.y, pos.z);
  gluSphere(q, rad, 15, 15);
  glPopMatrix();

}
//--------------------------------------------------------------------------
// 
// all spheres for rendering voronoi sites (particles)
//
void draw_spheres(vector<vec3d> &sites, float rad) {

  for (int i = 0; i < (int)sites.size(); i++) {

    glPushMatrix();
    glTranslatef(sites[i].x, sites[i].y, sites[i].z);
    gluSphere(q, rad, 10, 10);
    glPopMatrix();

  }

}
//--------------------------------------------------------------------------
// 
// point sprite for rendering voronoi sites (particles)
//
void draw_sprites(vector<vec3d> &sites, float size) {

  glPushAttrib(GL_ALL_ATTRIB_BITS);

//   glDisable(GL_DEPTH_TEST);
//   glEnable (GL_BLEND); 

  glColor3f(1.0, 1.0, 1.0);  // color doesn't matter, will be textured over

  glPointSize(size);

  glEnable(GL_TEXTURE_2D);
  glEnable(GL_POINT_SPRITE);
  glTexEnvf(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
  glBindTexture(GL_TEXTURE_2D, tex);
  glEnable(GL_POINT_SMOOTH);
  glBegin(GL_POINTS);
  for(int i = 0; i < (int)sites.size(); i++)
    glVertex3f(sites[i].x, sites[i].y, sites[i].z);
  glEnd();
  glDisable(GL_POINT_SPRITE);

//   glDisable (GL_BLEND); 
//   glEnable(GL_DEPTH_TEST);

  glPopAttrib();

}
//--------------------------------------------------------------------------
//
// axes
//
void draw_axes() {

  glPushMatrix();

  // x
  glColor3f(1.0, 0.0, 0.0); 
  glPushMatrix();
  glRotatef(90.0, 0.0, 1.0, 0.0);
  gluCylinder(q, size * 0.007, size * 0.007, sizes.x * 1.2, 15, 1);
  glTranslatef(0.0, 0.0, sizes.x * 1.2);
  gluCylinder(q, size * 0.015, 0.0, size * .015, 20, 1);
  glPopMatrix();

  // y
  glColor3f(0.0, 1.0, 0.0); 
  glPushMatrix();
  glRotatef(-90.0, 1.0, 0.0, 0.0);
  gluCylinder(q, size * 0.007, size * 0.007, sizes.y * 1.2, 15, 1);
  glTranslatef(0.0, 0.0, sizes.y * 1.2);
  gluCylinder(q, size * 0.015, 0.0, size * .015, 20, 1);
  glPopMatrix();

  // z
  glColor3f(0.0, 0.0, 1.0); 
  glPushMatrix();
  gluCylinder(q, size * 0.007, size * 0.007, sizes.z * 1.2, 15, 1);
  glTranslatef(0.0, 0.0, sizes.z * 1.2);
  gluCylinder(q, size * 0.015, 0.0, size * .015, 20, 1);
  glPopMatrix();

  glPopMatrix();

}
//--------------------------------------------------------------------------
//
// mouse button events
//
void mouse(int button, int state, int x, int y) {

  if (state == GLUT_DOWN) {

    press_x = x;
    press_y = y; 
    if (button == GLUT_LEFT_BUTTON)
      xform_mode = XFORM_ROTATE; 
    else if (button == GLUT_RIGHT_BUTTON) 
      xform_mode = XFORM_SCALE; 
    else if (button == GLUT_MIDDLE_BUTTON) 
      xform_mode = XFORM_TRANS; 

  }
  else if (state == GLUT_UP)
    xform_mode = XFORM_NONE; 

}
//--------------------------------------------------------------------------
//
// mouse motion events
//
void motion(int x, int y) {

  if (xform_mode == XFORM_ROTATE) {

    rot.x += (x - press_x) * rot_rate; 
    if (rot.x > 180)
      rot.x -= 360; 
    else if (rot.x < -180)
      rot.x += 360; 
    press_x = x; 
	   
    rot.y += (y - press_y) * rot_rate;
    if (rot.y > 180)
      rot.y -= 360; 
    else if (rot.y <-180)
      rot.y += 360; 
    press_y = y; 

  }
  else if (xform_mode == XFORM_TRANS) {

    trans.x += (x - press_x) * trans_rate; 
    trans.y -= (y - press_y) * trans_rate;  // subtract to reverse y dir.
    press_x = x;
    press_y = y; 

  }
  else if (xform_mode == XFORM_SCALE){

    float old_scale = scale;
    scale /= (1 + (y - press_y) * scale_rate);  // divided to reverse y dir.
    if (scale < 0) 
      scale = old_scale; 
    press_y = y; 

  }

  glutPostRedisplay(); 

}
//--------------------------------------------------------------------------
//
// keyboard events
//
void key(unsigned char key, int x, int y) {

  x = x; // quiet compiler warnings
  y = y;

  switch(key) {

  case 'q':  // quit
    exit(1);
    break; 
  case 'z':  // zoom mouse motion
    xform_mode = XFORM_SCALE; 
    break; 
  case 'a':  // panning mouse motion
    xform_mode = XFORM_TRANS; 
    break; 
  case 'r': // reset rotate, pan, zoom, viewport
    init_model();
    init_viewport(true);
    break;
  case 'b': // toggle block visibility
    block_mode = !block_mode;
    break;
  case 'f': // toggle fancy rendering
    draw_fancy = !draw_fancy;
    break;
  case 'c': // increase near clip plane
    near += 0.11;
    break;
  case 'C': // decrease near clip plane
    near -= 0.1;
    break;
  case 'v': // restrict (minimum) volume range
    min_vol += vol_step;
    fprintf(stderr, "Minimum volume = %.3lf\n", min_vol);
    filter_volume(min_vol);
    break;
  case 'V': //  expand (minimum) volume range
    min_vol -= vol_step;
    if (min_vol < 0.0)
      min_vol = 0.0;
    fprintf(stderr, "Minimum volume = %.3lf\n", min_vol);
    filter_volume(min_vol);
    break;
  case 'R': // reset volume range
    min_vol = 0.0;
    fprintf(stderr, "Minimum volume = %.3lf\n", min_vol);
    filter_volume(min_vol);
    break;
  case 's': // decrease volume step size
    vol_step *= 0.1;
    if (vol_step < 0.001)
      vol_step = 0.001;
    fprintf(stderr, "Volume step size = %.3lf\n", vol_step);
    break;
  case 'S': // increase volume step size
    vol_step *= 10.0;
    fprintf(stderr, "Volume step size = %.3lf\n", vol_step);
    break;
  default:
    break;

  }
}
//--------------------------------------------------------------------------
//
// filter volume
//
void filter_volume(float min_vol) {

  num_face_verts.clear();
  verts.clear();

  // package rendering data
  for (int i = 0; i < nblocks; i++) { // blocks

    int  n = 0;
    int m = 0;

    for (int j = 0; j < vblocks[i]->num_complete_cells; j++) { // cells

      for (int k = 0; k < vblocks[i]->num_cell_faces[j]; k++) { // faces

	if (vblocks[i]->vols[j] >= min_vol)
	    num_face_verts.push_back(vblocks[i]->num_face_verts[n]);

	for (int l = 0; l < vblocks[i]->num_face_verts[n]; l++) { // vertices

	  int v = vblocks[i]->face_verts[m];
	  vec3d s;
	  s.x = vblocks[i]->save_verts[3 * v];
	  s.y = vblocks[i]->save_verts[3 * v + 1];
	  s.z = vblocks[i]->save_verts[3 * v + 2];
	  m++;
	  if (vblocks[i]->vols[j] >= min_vol)
	      verts.push_back(s);

	} // vertices

	n++;

      } // faces

    } // cells

  } // blocks

}
//--------------------------------------------------------------------------
//
// timer events
//
void timer(int val) {

  val = val; // quiet compiler warning

  glutPostRedisplay();
  glutTimerFunc(10, timer, 0); 

}
//--------------------------------------------------------------------------
//
// reshape events
//
void reshape(int w, int h) {

  // update window and viewport size and aspect ratio
  win_size.x = w;
  win_size.y = h;

  init_viewport(false);

  glutPostRedisplay();

}
//--------------------------------------------------------------------------
//
// initialize model
//
void init_model() {

  // rotate
  rot.x = rot.y = 0.0;

  // translate
  trans.x = trans.y = 0.0;

  // scale (initial scale 1.5 makes the model fill the screen better)
  scale = 1.5 / size;

}
//--------------------------------------------------------------------------
//
// initialize viewport
//
// reset: true = first time or reset to initial viewport
//        false = modify existing viewport
//
void init_viewport(bool reset) {

  if (win_size.x > win_size.y) {
    aspect.x = 1.0;
    aspect.y = win_size.x / win_size.y;
    if (reset)
      trans.y -= (win_size.x - win_size.y) / win_size.y;
  }
  else {
    aspect.x = win_size.y / win_size.x;
    aspect.y = 1.0;
    if (reset)
      trans.x -= (win_size.y - win_size.x) / win_size.x;
  }

  if (!reset) {
    trans.x += (win_size.x - old_win_size.x) / old_win_size.x;
    trans.y += (win_size.y - old_win_size.y) / old_win_size.y;
  }
  else
    near = 0.1;

  old_win_size.x = win_size.x;
  old_win_size.y = win_size.y;

  glViewport(0, 0, win_size.x * aspect.x, win_size.y * aspect.y);

}
//--------------------------------------------------------------------------
//
// compute normal of a face using Newell's method
//
void ComputeNormal(vec3d *verts, int num_verts, vec3d &normal) {

  normal.x = 0.0;
  normal.y = 0.0;
  normal.z = 0.0;

  for (int i = 0; i < num_verts; i++) {
    int cur = i;
    int next = (i + 1) % num_verts;
    normal.x += (verts[cur].y - verts[next].y) * (verts[cur].z + verts[next].z);
    normal.y += (verts[cur].z - verts[next].z) * (verts[cur].x + verts[next].x);
    normal.z += (verts[cur].x - verts[next].x) * (verts[cur].y + verts[next].y);
  }

  float mag = sqrt(normal.x * normal.x + normal.y * normal.y +
		   normal.z * normal.z);
  // normalize
  normal.x /= mag;
  normal.y /= mag;
  normal.z /= mag;

  // direction is inward, need to invert
  normal.x *= -1.0;
  normal.y *= -1.0;
  normal.z *= -1.0;

}
//--------------------------------------------------------------------------
