// 1haskhsh.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include <GL/glut.h>
#include <GL/freeglut.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <string.h>
#include <cstring>

//#####-General Global Variables-########################
using namespace std;
float EPSILON = 0.00001;
GLsizei windowsWidth = 600, windowsHeight = 500;

int numberOfVertices = 0;
int x = 0;
int counter = 0;
int numberOfPolygons = 0;

int currentP = 0;
bool colorForLinesSelected[50] = {false};
bool colorForFillSelected[50] = {false};
bool design = false;
bool PolygonMode = false;
bool ClippingMode = false;
bool polygonInterrupted = false;
bool polygonNotNormal = false;
bool polygon_triangles_draw = false;
bool clipped_polygon_triangles_draw = false;
bool showClippedPolygonsTriangles = false;
bool UserInputMode = false;
bool draw3dPolygonsMode = false;
bool motionRectangle = false;
int keyTimesPressed = -1;//etsi wste otan to pataw prwth fora na emfanizontai ta trigwna gia na einai to mod2 = 0;
int keyTimesPressed1 = -1;
//alliws tha emfanizontai me to deutero pathma!

//#######################-Camera Global Variables-#################################################################
float angleOfView = 0.0f;
float lookX = 0.0f, lookY = 0.0f, lookZ = -1.0f;
float cameraPosX = 0.0f, cameraPosZ = 5.0f;
//#################################################################################################################

//variables for extrusion length from user
vector <char> z_deph;
int extrusion_length = 0;


//structures and vectors that hold polygon information such as polygons vertices , vertices of triangles of the polygons
struct vertex {
	int x;
	int y;
};
vector <vertex> vertices, clip, clipPoints, clippedPoly;
int sizeClippedPoly;
vector <vertex> triangle_vertices ;
vector <vector<vertex> > triangle_database, clipped_polygons_triangles;
vector <vector<vertex> > database, clippedPolygons;
vertex vertexCourdinate;

vector <int> polygonSizes;//kathe stoixeio einai to megethos(poses korufes) tou kathe sxediasmenou polugwnou

vertex motionClippingRectPoints[2];

//structure for the color
struct Color {
	double red;
	double green;
	double blue;
};
Color colorForFill[50], colorForLines[50];//colors for lines and fill for polygons

//structures and vector that hiold information for the edges of polygons
struct Edge{
	vertex firstVertex;
	vertex secVertex;
};
Edge edge;
vector <Edge> edges;
vector <vector<Edge> > polygonsEdges, clippedPolygonsEdges;

//every time a new polygon is created from the user we find its edges
void findpolygonsEdges() {
	Edge temp_edge;
	vector <Edge> current_edges;
	for (int num = 0; num < database[numberOfPolygons -1].size(); num++) {
		if (num == (database[numberOfPolygons - 1].size() - 1)) {
			temp_edge.firstVertex.x = database[numberOfPolygons -1].back().x;
			temp_edge.firstVertex.y = database[numberOfPolygons -1].back().y;
			temp_edge.secVertex.x = database[numberOfPolygons -1][0].x;
			temp_edge.secVertex.y = database[numberOfPolygons -1][0].y;			
		}
		else {
			temp_edge.firstVertex.x = database[numberOfPolygons -1][num].x;
			temp_edge.firstVertex.y = database[numberOfPolygons -1][num].y;
			temp_edge.secVertex.x = database[numberOfPolygons -1][num + 1].x;
			temp_edge.secVertex.y = database[numberOfPolygons -1][num + 1].y;
		}
		current_edges.push_back(temp_edge);
	}
	/*
	for (int i = 0; i < current_edges.size(); i++) {
				cout << current_edges[i].firstVertex.x << " ";
				cout << current_edges[i].firstVertex.y << " ";
				cout << current_edges[i].secVertex.x << " ";

				cout << current_edges[i].secVertex.y << " ";
			
	}
	*/
	polygonsEdges.push_back(current_edges);
}
void findClippedPolygonsEdges(void) {
	Edge temp_edge;
	vector <Edge> current_edges;
	for (int num = 0; num < clippedPolygons.size(); num++) {
		for (int j = 0; j < polygonSizes[num]; j++){
			if (num == (polygonSizes[num] - 1)) {
				temp_edge.firstVertex.x = clippedPolygons[num].back().x;
				temp_edge.firstVertex.y = clippedPolygons[num].back().y;
				temp_edge.secVertex.x = clippedPolygons[num][0].x;
				temp_edge.secVertex.y = clippedPolygons[num][0].y;			
			}
			else {
				temp_edge.firstVertex.x = clippedPolygons[num][j].x;
				temp_edge.firstVertex.y = clippedPolygons[num][j].y;
				temp_edge.secVertex.x = clippedPolygons[num][j + 1].x;
				temp_edge.secVertex.y = clippedPolygons[num][j + 1].y;
			}
			current_edges.push_back(temp_edge);
		}
		clippedPolygonsEdges.push_back(current_edges);
		current_edges.clear();
	}
	
}
//#########################################-Functions For Triangulating Polygons and clipped Polygons-########################################
float Area(vector <vertex> contour) {
	int n = contour.size();
	float A = 0.0f;
	for (int p = n-1, q = 0; q < n; p = q++)
		A+= contour[p].x * contour[q].y - contour[q].x * contour[p].y;
	return A * 0.5f;
}

float AreaForClippingPolygons(vector <vertex> contour, int polysize) {
	int n = polysize;
	float A = 0.0f;
	for (int p = n-1, q = 0; q < n; p = q++)
		A+= contour[p].x * contour[q].y - contour[q].x * contour[p].y;
	return A * 0.5f;
}

bool InsideTriangle(int Ax, int Ay,
                      int Bx, int By,
                      int Cx, int Cy,
                      int Px, int Py)

{
  int ax, ay, bx, by, cx, cy, apx, apy, bpx, bpy, cpx, cpy;
  int cCROSSap, bCROSScp, aCROSSbp;

  ax = Cx - Bx;  ay = Cy - By;
  bx = Ax - Cx;  by = Ay - Cy;
  cx = Bx - Ax;  cy = By - Ay;
  apx= Px - Ax;  apy= Py - Ay;
  bpx= Px - Bx;  bpy= Py - By;
  cpx= Px - Cx;  cpy= Py - Cy;

  aCROSSbp = ax*bpy - ay*bpx;
  cCROSSap = cx*apy - cy*apx;
  bCROSScp = bx*cpy - by*cpx;

  return ((aCROSSbp >= 0) && (bCROSScp >= 0) && (cCROSSap >= 0));
}

bool Snip(vector <vertex> contour ,int u ,int v, int w, int n, int *V)
{
  int p;
  int Ax, Ay, Bx, By, Cx, Cy, Px, Py;

  Ax = contour[V[u]].x;
  Ay = contour[V[u]].y;

  Bx = contour[V[v]].x;
  By = contour[V[v]].y;

  Cx = contour[V[w]].x;
  Cy = contour[V[w]].y;

  if ( EPSILON > (((Bx-Ax)*(Cy-Ay)) - ((By-Ay)*(Cx-Ax))) ) return false;

  for (p=0;p<n;p++)
  {
    if( (p == u) || (p == v) || (p == w) ) continue;
    Px = contour[V[p]].x;
    Py = contour[V[p]].y;
    if (InsideTriangle(Ax,Ay,Bx,By,Cx,Cy,Px,Py)) return false;
  }

  return true;
}

vector<vertex> ProcessClippedPolygons(vector <vertex> vertices, int polysize) {
	int n = polysize;
	vertex temp;
	vector <vertex> triangles;
	//if (n < 3) exit(0);

	int *V = new int[n];
	if (0.0f < AreaForClippingPolygons(vertices, polysize))
		for(int v=0; v<n; v++) V[v] = v;
	else
		for(int v=0; v<n; v++) V[v] = (n-1)-v;
	int nv = n;
	int count = 2*nv;
	for(int m=0, v=nv-1; nv>2; )
	{
		/* if we loop, it is probably a non-simple polygon */
		if (0 >= (count--))
		{
		  //** Triangulate: ERROR - probable bad polygon!
		  //exit(0);
			polygonNotNormal = true; 
			return triangles;
		}

		/* three consecutive vertices in current polygon, <u,v,w> */
		int u = v  ; if (nv <= u) u = 0;     /* previous */
		v = u+1; if (nv <= v) v = 0;     /* new v    */
		int w = v+1; if (nv <= w) w = 0;     /* next     */

		if ( Snip(vertices,u,v,w,nv,V) )
		{
		  int a,b,c,s,t;

		  /* true names of the vertices */
		  a = V[u]; b = V[v]; c = V[w];

		  /* output Triangle */
		  triangles.push_back( vertices[a] );
		  triangles.push_back( vertices[b] );
		  triangles.push_back( vertices[c] );

		  m++;

		  /* remove v from remaining polygon */
		  for(s=v,t=v+1;t<nv;s++,t++) V[s] = V[t]; nv--;

		  /* resest error detection counter */
		  count = 2*nv;
		}
	}
	delete V;
	return triangles ;
}

vector<vertex> Process(vector <vertex> vertices) {
	int n = vertices.size();
	vertex temp;
	vector <vertex> triangles;
	//if (n < 3) exit(0);

	int *V = new int[n];
	if (0.0f < Area(vertices))
		for(int v=0; v<n; v++) V[v] = v;
	else
		for(int v=0; v<n; v++) V[v] = (n-1)-v;
	int nv = n;
	int count = 2*nv;
	for(int m=0, v=nv-1; nv>2; )
	{
		/* if we loop, it is probably a non-simple polygon */
		if (0 >= (count--))
		{
		  //** Triangulate: ERROR - probable bad polygon!
		  	
			polygonNotNormal = true; 
			return triangles;
		}

		/* three consecutive vertices in current polygon, <u,v,w> */
		int u = v  ; if (nv <= u) u = 0;     /* previous */
		v = u+1; if (nv <= v) v = 0;     /* new v    */
		int w = v+1; if (nv <= w) w = 0;     /* next     */

		if ( Snip(vertices,u,v,w,nv,V) )
		{
		  int a,b,c,s,t;

		  /* true names of the vertices */
		  a = V[u]; b = V[v]; c = V[w];

		  /* output Triangle */
		  triangles.push_back( vertices[a] );
		  triangles.push_back( vertices[b] );
		  triangles.push_back( vertices[c] );

		  m++;

		  /* remove v from remaining polygon */
		  for(s=v,t=v+1;t<nv;s++,t++) V[s] = V[t]; nv--;

		  /* resest error detection counter */
		  count = 2*nv;
		}
	}
	delete V;
	return triangles ;
}
//##########################################################################################################################################

//#################################################-Initialize , change to 3dSpace And Reshape Functions-#############################################################
void initializeColors(void) {// dinei arxiko xrwma gemismatos kai grammhs gia kaue polygwno
				// to leuko kai mauro antristoixa
	for (int i = 0; i < 50; i++) {
		colorForFill[i].red = 1.0;
		colorForFill[i].green = 1.0;
		colorForFill[i].blue = 1.0;
		colorForLines[i].red = 0.0;
		colorForLines[i].green = 0.0;
		colorForLines[i].blue = 0.0;
	}
}
void initializeArrays(void) {
	for (int i = 0; i < 50; i++) {
		colorForLinesSelected[i] = false;
		colorForFillSelected[i] = false;
	}
}
void init(void)
{
	glClearColor(1.0, 1.0, 1.0, 0.0);
	// 2d enviroment. draw just normal polygons
	glViewport(0, 0, windowsWidth, windowsHeight);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0f, windowsWidth, 0.0f, windowsHeight, 0.0f, 1.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glEnable(GL_LINE_SMOOTH);
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA_SATURATE, GL_ONE);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	initializeColors();
	
}

void reshape(GLint width, GLint height) {
	glutReshapeWindow(600, 500);
}

void changeInto3dEnviroment(void) {//3d enviroment for navigate threw extruded polygons
	glEnable(GL_DEPTH_TEST);
	glViewport(0, 0, windowsWidth, windowsHeight);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, float(windowsWidth)/float(windowsHeight), 1.0, 100.0);
	glMatrixMode(GL_MODELVIEW);
	gluLookAt(cameraPosX, 0.0f, cameraPosZ, lookX, lookY, lookZ, 0.0f, 1.0f, 0.0f);
}
//######################################################################################################################

//######################################-Functions For checking If a polygon is Self-Intersect-############################################


bool isIntersecting(vertex p1, vertex p2, vertex q1, vertex q2) {
    return (float(((q1.x-p1.x)*(p2.y-p1.y) - (q1.y-p1.y)*(p2.x-p1.x)))
            * float(((q2.x-p1.x)*(p2.y-p1.y) - (q2.y-p1.y)*(p2.x-p1.x))) < 0)
            &&
           (float(((p1.x-q1.x)*(q2.y-q1.y) - (p1.y-q1.y)*(q2.x-q1.x)))
            * float(((p2.x-q1.x)*(q2.y-q1.y) - (p2.y-q1.y)*(q2.x-q1.x))) < 0);
}



/*
bool onLine(vertex p1, vertex p2, vertex p) {
	if (p.x <= max(p1.x, p2.x) && p.x <= min(p1.x, p2.x) && (p.y <= max(p1.y, p2.y) && p.y <= min(p1.y, p2.y)))
		return true;
	return false;
}

int direction(vertex a, vertex b, vertex c) {
	int val = (b.y - a.y)*(c.x-b.x) - (b.x-a.x)*(c.y-b.y);
	
	if (val == 0)
		return 0;
	else if (val < 0)
		return 2;
	else
		return 1;	
}

bool isIntersecting(vertex a, vertex b, vertex c, vertex d) {
	int dir1 = direction(a, b, c);
	int dir2 = direction(a, b, d);
	int dir3 = direction(c, d, a);
	int dir4 = direction(c, d, b);
	
	if (dir1 != dir2 && dir3 != dir4)
		return true;
	if (dir1 == 0 && onLine (a, b, c))
		return true;
	if (dir2 == 0 && onLine (a, b, d))
		return true;
	if (dir3 == 0 && onLine (c, d, a))
		return true;
	if (dir4 == 0 && onLine (c, d, b))
		return true;
	
	return false;
}
*/
void findEdges(void) {
	int i;
	for (i = 0; i < vertices.size(); i++) {
		if (i == vertices.size() - 1) 
			break;
		else {
			edge.firstVertex = vertices[i];
			edge.secVertex = vertices[i+1];
			edges.push_back(edge);
		}
	}
}


void polygonCheckForIntersection(void) {
	findEdges();
	for (int i = 0; i < edges.size(); i++) {
		for (int j = i+1; j < edges.size(); j++) {
			if (isIntersecting(edges[i].firstVertex, edges[i].secVertex, edges[j].firstVertex, edges[j].secVertex)) {
				polygonInterrupted = true;
			}
					
		}
	}
}

void lastPolygonCheckForIntersection(void) {
	findEdges();
	Edge lastEdge;
	lastEdge.firstVertex.x = vertices[vertices.size()-1].x;
	lastEdge.firstVertex.y = vertices[vertices.size()-1].y;
	lastEdge.secVertex.x = vertices[0].x;
	lastEdge.secVertex.y = vertices[0].y;
	for (int i = 0; i < edges.size(); i++) {
		if (isIntersecting(edges[i].firstVertex, edges[i].secVertex, lastEdge.firstVertex, lastEdge.secVertex)) {
			polygonInterrupted = true;
		}				
	}
}


//#############################-Functions For Drawing User Polygons-#####################################################################
void drawPreviusPolygons(void)
{

	
	if (numberOfPolygons == 0)
		return;
	for (int i = 0; i < database.size(); i++) {
		if (database[i].size() == 3) {
			glLineWidth(4.0);
			glBegin(GL_LINE_LOOP);
			glColor3f(colorForLines[i].red, colorForLines[i].green, colorForLines[i].blue);
			glVertex2i(database[i][0].x, database[i][0].y);
			glVertex2i(database[i][1].x, database[i][1].y);
			glVertex2i(database[i][2].x, database[i][2].y);
			glEnd();
			glLineWidth(4.0);
			glColor3f(colorForFill[i].red, colorForFill[i].green, colorForFill[i].blue);
			glBegin(GL_TRIANGLES);
			glVertex2i(database[i][0].x, database[i][0].y);
			glVertex2i(database[i][1].x, database[i][1].y);
			glVertex2i(database[i][2].x, database[i][2].y);
			glEnd();

		}
		else {
			glLineWidth(4.0);
			glBegin(GL_LINE_LOOP);
			glColor3f(colorForLines[i].red, colorForLines[i].green, colorForLines[i].blue);
			for (int num = 0; num < database[i].size(); num++) {
				glVertex2i(database[i][num].x, database[i][num].y);
			}
			glEnd();
		}
	}
}



void drawCurrentPolygon(void)
{
		glPointSize(6.0);
		glBegin(GL_POINTS);
		glColor3f(colorForLines[currentP].red, colorForLines[currentP].green, colorForLines[currentP].blue);
		for (int vertexCounter = 0; vertexCounter < vertices.size(); vertexCounter++)
			glVertex2i(vertices[vertexCounter].x, vertices[vertexCounter].y);
		glEnd();
		glLineWidth(4.0);
		glBegin(GL_LINE_STRIP);
		glColor3f(colorForLines[currentP].red, colorForLines[currentP].green, colorForLines[currentP].blue);
		for (int vertexCounter = 0; vertexCounter < vertices.size(); vertexCounter++)
			glVertex2i(vertices[vertexCounter].x, vertices[vertexCounter].y);
		glEnd();
			
}

void drawPolygonsAsTriangles(void) {
	if (numberOfPolygons == 0)
		return;
	
	for (int i = 0; i < triangle_database.size(); i++) {
			for(int j = 0; j < triangle_database[i].size()/3; j++) {
				glColor3f(colorForFill[i].red, colorForFill[i].green, colorForFill[i].blue);
				glBegin(GL_TRIANGLES);
				glVertex2i(triangle_database[i][j*3+0].x, triangle_database[i][j*3+0].y);
				glVertex2i(triangle_database[i][j*3+1].x, triangle_database[i][j*3+1].y);
				glVertex2i(triangle_database[i][j*3+2].x, triangle_database[i][j*3+2].y);
				glEnd();	
			}
		}
		
}


void drawPolygonsTriangles(void) {
		
		for (int i = 0; i < triangle_database.size(); i++) {
			for(int j = 0; j < triangle_database[i].size()/3; j++) {
				glLineWidth(3.0);
				glBegin(GL_LINE_LOOP);
				glColor3f(0.0, 1.0, 0.0);
				glVertex2i(triangle_database[i][j*3+0].x, triangle_database[i][j*3+0].y);
				glVertex2i(triangle_database[i][j*3+1].x, triangle_database[i][j*3+1].y);
				glVertex2i(triangle_database[i][j*3+2].x, triangle_database[i][j*3+2].y);
				glEnd();	
			}
		}
	
}
//################################################################################################################################

//#####################-Functions for Calculating new Polygons after clipped by a a clipped rectangle-#########################

//calculate the x point of intersection between the line of the clipping area and the edge of the polygon
int findXpointOfIntersection(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4) {
	int x = (x1*y2 - y1*x2) * (x3-x4) - 
              (x1-x2) * (x3*y4 - y3*x4); 
	int y = (x1-x2) * (y3-y4) - (y1-y2) * (x3-x4); 

	return int(x/y);
}

//calculate the y point of intersection between the line of the clipping area and the edge of the polygon
int findYpointOfIntersection(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4) {
	int x = (x1*y2 - y1*x2) * (y3-y4) - 
              (y1-y2) * (x3*y4 - y3*x4); 
	int y = (x1-x2) * (y3-y4) - (y1-y2) * (x3-x4); 

	return int(x/y);
}



void clip_polygon(int x1, int y1, int x2, int y2) {//shmeia pou kathorizoun thn akmh tou orthogwnioy apokophs
	vertex temp;                               //
	vector <vertex> newPoints;
	for (int i = 0; i < polygonSizes[counter]; i++) {// gia kathe akmh tou polugwnou
		int j = (i+1) % polygonSizes[counter];
		int fp1 = clippedPolygons[counter][i].x, sp1 = clippedPolygons[counter][i].y;
		int fp2 = clippedPolygons[counter][j].x, sp2 = clippedPolygons[counter][j].y;

		int first_pos = (x2-x1) * (sp1-y1) - (y2-y1) * (fp1-x1); //position of first point
        	int sec_pos = (x2-x1) * (sp2-y1) - (y2-y1) * (fp2-x1); //position of second point

		//analoga me thn akmh tou orthogwnioy apokophs kai thn akmh tou polugwnou
		if (first_pos < 0 && sec_pos <0) {// otan kai ta 2 shmeia einai mesa 
			temp.x = fp2;
			temp.y = sp2;
			newPoints.push_back(temp); //prosthese mono to deytero shmeio
		}
		else if (first_pos >= 0  && sec_pos < 0) {// otan mono to prwo shmeio einai eksw
			//tote prosthetetai mono to shmeio tomhs kai to deytero shmeio 
			temp.x = findXpointOfIntersection(x1, y1, x2, y2, fp1, sp1, fp2, sp2);
			temp.y = findYpointOfIntersection(x1, y1, x2, y2, fp1, sp1, fp2, sp2);
			newPoints.push_back(temp);
			temp.x = fp2;
			temp.y = sp2;
			newPoints.push_back(temp);
		}
		else if (first_pos < 0  && sec_pos >= 0) {// otan mono to deytero shmeio einai eksw
			temp.x = findXpointOfIntersection(x1, y1, x2, y2, fp1, sp1, fp2, sp2);
			temp.y = findYpointOfIntersection(x1, y1, x2, y2, fp1, sp1, fp2, sp2);
			newPoints.push_back(temp);//prostithete mono to shmeio tomhs
		}
		else { // otan kai ta dyo shmeia vriskontai eksw 
			//den prostithentai korufes
		}
	}

	polygonSizes[counter] = newPoints.size();//ananenwnoume to plhthos twn koryfwn twn polygwnwn
	for (int i = 0; i < polygonSizes[counter]; i++) { //antigrafoume tis trexouses korufes pou vrhkame
							  //stis prohgoumenes korufes pou orizan to polugwno
		clippedPolygons[counter][i].x = newPoints[i].x;
		clippedPolygons[counter][i].y = newPoints[i].y;
	}
	
}


void Sutherland_Hodgman(vector<vertex> clipperPoints, int clipperSize) {
	for (int i = 0; i < clippedPolygons.size(); i++) {
			for (int j = 0; j < clipperSize; j++) { //gia kathe akmh tou orthogwnioy apokophs
				int k = (j+1) % clipperSize;    //apokopse to polygwno trexontas oles tis akmes tou polugwnou
				clip_polygon(clipperPoints[j].x, clipperPoints[j].y, clipperPoints[k].x, clipperPoints[k].y);
			}
			counter++;
	}
}

//##################################################################################################################

//########################################-Function For Drawing Clippin-Rectangle###################################
void drawClippingRectangle(void) {
	if (clip.size() > 1) {
		glBegin(GL_LINE_LOOP);
		glColor3f(0.0, 0.0, 0.0);
		glVertex2i(clip[0].x, clip[0].y);
		glVertex2i(clip[1].x, clip[0].y);
		glVertex2i(clip[1].x, clip[1].y);
		glVertex2i(clip[0].x, clip[1].y);
		glEnd();
	}
}

void drawMotionClippingRectangle(void) {
	if (!motionRectangle) {
		glBegin(GL_LINE_LOOP);
		glColor3f(0.0, 0.0, 0.0);
		glVertex2i(motionClippingRectPoints[0].x, motionClippingRectPoints[0].y);
		glVertex2i(motionClippingRectPoints[1].x, motionClippingRectPoints[0].y);
		glVertex2i(motionClippingRectPoints[1].x, motionClippingRectPoints[1].y);
		glVertex2i(motionClippingRectPoints[0].x, motionClippingRectPoints[1].y);
		glEnd();
	}
}
void drawRectanglePoints(void) {
	glPointSize(5.0);
	glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < clip.size(); i++) 
		glVertex2i(clip[i].x, clip[i].y);
	glEnd();
}
//###################################################################################################################

//########################-Functions for drawing and triangulating Clipped Polygons-#########################################################
void drawClippedPolygons(void) {
	for (int i = 0; i < clippedPolygons.size(); i++) {
		glLineWidth(4.0);
		glBegin(GL_LINE_LOOP);
		glColor3f(0.0, 1.0, 0.0);
		for (int num = 0; num < polygonSizes[i]; num++) {
			glVertex2i(clippedPolygons[i][num].x, clippedPolygons[i][num].y);
		}
		glEnd();
	}
}
void drawClippedPolygonsTriangles(void) {
	if (showClippedPolygonsTriangles) {
		for (int i = 0; i < clipped_polygons_triangles.size(); i++) {
				for(int j = 0; j < clipped_polygons_triangles[i].size()/3; j++) {
					glLineWidth(3.0);
					glColor3f(0.0, 1.0, 0.0);
					glBegin(GL_LINE_LOOP);
					glVertex2i(clipped_polygons_triangles[i][j*3+0].x, clipped_polygons_triangles[i][j*3+0].y);
					glVertex2i(clipped_polygons_triangles[i][j*3+1].x, clipped_polygons_triangles[i][j*3+1].y);
					glVertex2i(clipped_polygons_triangles[i][j*3+2].x, clipped_polygons_triangles[i][j*3+2].y);
					glEnd();	
				}
		}
	}
}

void drawClippedPolygonsAsTriangles(void) {
	for (int i = 0; i < clipped_polygons_triangles.size(); i++) {
			for(int j = 0; j < clipped_polygons_triangles[i].size()/3; j++) {
				glColor3f(colorForFill[i].red, colorForFill[i].green, colorForFill[i].blue);
				glBegin(GL_TRIANGLES);
				glVertex2i(clipped_polygons_triangles[i][j*3+0].x, clipped_polygons_triangles[i][j*3+0].y);
				glVertex2i(clipped_polygons_triangles[i][j*3+1].x, clipped_polygons_triangles[i][j*3+1].y);
				glVertex2i(clipped_polygons_triangles[i][j*3+2].x, clipped_polygons_triangles[i][j*3+2].y);
				glEnd();	
			}
		}
}

void triangulateClippedPolygons(void) {
	
	for (int i = 0; i < clippedPolygons.size(); i++) {
		if (polygonSizes[i] > 3) {
			triangle_vertices = ProcessClippedPolygons(clippedPolygons[i], polygonSizes[i]);
			clipped_polygons_triangles.push_back(triangle_vertices);
		}
	}
}
//####################################################################################################################################
//

//############################################-draw Text Area Input as a Rectangle and User Input-#################################################
void drawTextArea(void) {
	//glPushMatrix();
	glBegin(GL_LINE_LOOP);
	glColor3f(0.0, 0.0, 0.0);
	glVertex2i(250, 200);
	glVertex2i(350, 200);
	glVertex2i(350, 245);
	glVertex2i(250, 245);
	glEnd();
	//glPopMatrix();
}

void showScreenText(float x, float y, const char *screenText) {
	
	glPushMatrix();
	glColor3f(0.0f, 0.0f, 0.0f);
	glTranslatef(x, y, 0);
	glScalef(0.2, 0.2, 0.2);
	for (const char* p = screenText; *p; p++)
		glutStrokeCharacter(GLUT_STROKE_ROMAN, *p);
	glPopMatrix();
}

void showUserInput(float x, float y, vector <char> z_value) {
	
	glPushMatrix();
	glColor3f(1.0f, 0.0f, 0.0f);
	glTranslatef(x, y, 0);
	glScalef(0.2, 0.2, 0.2);
	for (int i = 0; i < z_value.size(); i++)
		glutStrokeCharacter(GLUT_STROKE_ROMAN, z_value[i]);
	glPopMatrix();
}
//#################################-Functions for Draw Polygons in 3 Dimensions-##################################################################

void drawEdgesAsQuads(void) {
	for (int i = 0; i < clippedPolygonsEdges.size(); i ++) {
		for (int j = 0; j < clippedPolygonsEdges[i].size(); j++) {	
			glBegin(GL_QUADS);
			glColor3f(colorForLines[i].red, colorForLines[i].green, colorForLines[i].blue);
			glVertex3i(clippedPolygonsEdges[i][j].firstVertex.x, clippedPolygonsEdges[i][j].firstVertex.y, 0);
			glVertex3i(clippedPolygonsEdges[i][j].secVertex.x, clippedPolygonsEdges[i][j].secVertex.y, 0);
			glVertex3i(clippedPolygonsEdges[i][j].secVertex.x, clippedPolygonsEdges[i][j].secVertex.y, -1*extrusion_length);
			glVertex3i(clippedPolygonsEdges[i][j].firstVertex.x, clippedPolygonsEdges[i][j].firstVertex.y, -1*extrusion_length);
			glEnd();
						
		}
	}
}

void drawClippedPolygonsAsTriangles3d(void) {
	for (int i = 0; i < clipped_polygons_triangles.size(); i++) {
			for(int j = 0; j < clipped_polygons_triangles[i].size()/3; j++) {
				glColor3f(colorForFill[i].red, colorForFill[i].green, colorForFill[i].blue);
				glBegin(GL_TRIANGLES);
				glVertex3i(clipped_polygons_triangles[i][j*3+0].x, clipped_polygons_triangles[i][j*3+0].y, 0);
				glVertex3i(clipped_polygons_triangles[i][j*3+1].x, clipped_polygons_triangles[i][j*3+1].y, 0);
				glVertex3i(clipped_polygons_triangles[i][j*3+2].x, clipped_polygons_triangles[i][j*3+2].y, 0);
				glEnd();	
			}
		}
}
void draw3dimPolygons(void) {
	drawEdgesAsQuads();
	drawClippedPolygonsAsTriangles3d();
}

////######################################-Actuall Function For Drawing the Scene-#####################################################
void drawScene(void) {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if (PolygonMode) {
		drawPreviusPolygons();
		drawCurrentPolygon();
		drawPolygonsAsTriangles();
		if (polygon_triangles_draw && keyTimesPressed % 2 == 0) {
			drawPolygonsTriangles();
			drawPreviusPolygons();
		}
		else
			drawPreviusPolygons();
	}

	else if (ClippingMode) {
		drawPreviusPolygons();
		drawRectanglePoints();
		drawMotionClippingRectangle();
		drawClippingRectangle();
		drawPolygonsAsTriangles();
		drawClippedPolygons();
		drawClippedPolygonsAsTriangles();
		if(clipped_polygon_triangles_draw && keyTimesPressed1 % 2 == 0) {
			drawClippedPolygonsTriangles();		
		}
		else 
			drawClippedPolygons();
		
	}
	
	else if (UserInputMode) {
		drawTextArea();
		showScreenText(170.0, 250.0, "Give extrude size pls");
		showUserInput(295.0, 210.0, z_deph);
	}
	
	else if (draw3dPolygonsMode) {
		gluLookAt(cameraPosX, 0.0f, cameraPosZ, lookX, lookY, lookZ, 0.0f, 1.0f, 0.0f);
		draw3dimPolygons();
	}
	
	glFlush();
	
}

void allocateSpace(void) {
	vertex temp;
	for (int i = 0; i < clippedPolygons.size(); i++) {
		for (int j = 0; j < 100; j++) {
			temp.x = 0;
			temp.y = 0;
			clippedPolygons[i].push_back(temp);
		}
	}
}

void findPolygonsSizes(void) {
	for (int i = 0; i < clippedPolygons.size(); i++) {
		polygonSizes.push_back(clippedPolygons[i].size());
	}
}



void processKey(unsigned char key, int x, int y) {
	if (UserInputMode) {
		if (key == 13) {
			UserInputMode = false;
			draw3dPolygonsMode = true;
			for (int i = 0; i < z_deph.size(); i++) {
				extrusion_length = extrusion_length * 10 + (z_deph[i] - '0');
			}
			cout << extrusion_length << " ";
			changeInto3dEnviroment();//change to gluPerspective mode so we can see 3d polygons...
		}
		else if (key == 127) {
			if (z_deph.size() > 0)
				z_deph.pop_back();
		}
		else if (z_deph.size() < 3) {
			if (key == 48) {
				z_deph.push_back('0');
			}
			else if (key == 49) {
				z_deph.push_back('1');
			}
			else if (key == 50) {
				z_deph.push_back('2');
			}
			else if (key == 51) {
				z_deph.push_back('3');
			}
			else if (key == 52) {
				z_deph.push_back('4');
			}
			else if (key == 53) {
				z_deph.push_back('5');
			}
			else if (key == 54) {
				z_deph.push_back('6');
			}
			else if (key == 55) {
				z_deph.push_back('7');
			}
			else if (key == 56) {
				z_deph.push_back('8');
			}
			else if (key == 57) {
				z_deph.push_back('9');
			}
		
		}
	}
	if ((key == 84 || key == 116) && PolygonMode) {
		polygon_triangles_draw = true;
		keyTimesPressed ++;
	}
	else if ((key == 84 || key == 116) && ClippingMode) {
			clipped_polygon_triangles_draw = true;
			keyTimesPressed1 ++;
	}
		
	glutPostRedisplay();
}

void specialKeys(int key, int x, int y) {
	if (key == GLUT_KEY_UP) {
		cameraPosX += lookX * 0.1f;
		cameraPosZ += lookZ * 0.1f;
	}
	else if (key == GLUT_KEY_DOWN) {
		cameraPosX -= lookX * 0.1f;
		cameraPosZ -= lookZ * 0.1f;
	}
	else if (key == GLUT_KEY_LEFT) {
		angleOfView -= 0.1f;
		lookX = sin(angleOfView);
		lookZ = -cos(angleOfView);
	}
	else if (key == GLUT_KEY_RIGHT) {
		angleOfView += 0.1f;
		lookX = sin(angleOfView);
		lookZ = -cos(angleOfView);
	}
	glutPostRedisplay();
}
//################################-Function For Mouse-################################################################
void mouse(GLint button, GLint state, GLint x, GLint y)
{

	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		//previusX = x;
		//previusY = y;
		//move =  state == GLUT_UP;
		if (state == GLUT_DOWN) {
			if (design) {
				vertexCourdinate.x = x;
				vertexCourdinate.y = windowsHeight - y;
				vertices.push_back(vertexCourdinate);
				numberOfVertices++;
				if (numberOfVertices > 3) {//an o xrhsths exei dwsei 3 oi perissoteres akmes elekse an oi uparxouses akmes temnontai. ana 2.
					polygonCheckForIntersection();
					if (polygonInterrupted) { // an uparxoun autotemnouses akmes diegrapse oti korufes exei dwsei o xrhsths.
						vertices.clear();
						numberOfVertices = 0;
						edges.clear();
						polygonInterrupted = false;
					}
					else
						edges.clear();
					
				}
				
			}
		
			else if (ClippingMode && clip.size() < 2) {//gia na mhn apothikeytei kai zwgrafistei kati allo meta thn epilogh deuterhs afou to clip size that einai 2
				vertexCourdinate.x = x;// korufhs gia to orthogwnio apokophs apo ton xrhsth
				vertexCourdinate.y = windowsHeight - y;
				clip.push_back(vertexCourdinate);
				motionClippingRectPoints[0].x = x;
				motionClippingRectPoints[0].y = windowsHeight - y;
				motionClippingRectPoints[1].x = x;
				motionClippingRectPoints[1].y = windowsHeight - y;
				if (clip.size() == 2) {//exoyn epilegei oi 2 korufes apo ton xrhsth tou orthogwnioy apokophs
							// opote arxise thn diadikasia apokophs kai trigwnopoihshs twn apokomenwn polugwnwn

					//apothikeuese tis 4 korufes pou orizoun to polugwno apokophs
					vertex x1,x2,x3,x4;
					x1.x = clip[0].x;
					x1.y = clip[0].y;
					clipPoints.push_back(x1);
					x2.x = clip[0].x;
					x2.y = clip[1].y;
					clipPoints.push_back(x2);
					x3.x = clip[1].x;
					x3.y = clip[1].y;
					clipPoints.push_back(x3);
					x4.x = clip[1].x;
					x4.y = clip[0].y;
					clipPoints.push_back(x4);
					clippedPolygons = database;
					findPolygonsSizes();
					allocateSpace();
					Sutherland_Hodgman(clipPoints, clipPoints.size());
					findClippedPolygonsEdges();
					triangulateClippedPolygons();
					showClippedPolygonsTriangles = true;
					motionRectangle = true;//stop draw a motion rectangle again after u already choose ur clipping rectangle area
					
				}

			}
		}
		else if (state == GLUT_UP) {
			
			motionClippingRectPoints[1].x = x;
			motionClippingRectPoints[1].y = windowsHeight - y;
			glutPostRedisplay();
		}
	}
	else if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
		if (PolygonMode == true && numberOfVertices > 2) {
			lastPolygonCheckForIntersection();
			if (!polygonInterrupted) {
				if (numberOfVertices > 3)//an den einai trigwno trigonopoihse to
					triangle_vertices = Process(vertices);
				if (!polygonNotNormal) {//second check if a polygon is interrupted actually
					database.push_back(vertices);
					triangle_database.push_back(triangle_vertices);
					numberOfPolygons++;//polygon is not self intersect so can be counted
					findpolygonsEdges();
						
					design = false;
					colorForLinesSelected[currentP] = true;
					colorForFillSelected[currentP] = true;
					currentP++;
				}
				else {
					
					polygonNotNormal = false;
				}
				
			}
			else {
				cout << "op" << "";
				polygonInterrupted  = false;
				//exit(0);
			}
			
			vertices.clear();
			numberOfVertices = 0;
			edges.clear();
			
		}
		
	}
	glutPostRedisplay();
}
//################################################################################################################################
void Motion(GLint x, GLint y) {
	if (ClippingMode) {
		motionClippingRectPoints[1].x = x;
		motionClippingRectPoints[1].y = windowsHeight - y;
		glutPostRedisplay();
	}
}

//#######################################################-Menus-##################################################################
void menuOption(GLint selectedOption)
{
	switch (selectedOption) {
	case 1: design = true;
		PolygonMode = true;
		colorForLinesSelected[currentP] = true;//exei epilegei automata to mauro gia tis grammes
		colorForFillSelected[currentP] = true;//exei epilegei automata to aspro gia gemisma
		break;
	case 2: ClippingMode = true;
		PolygonMode = false;
		break;
	case 3: UserInputMode = true;
		PolygonMode = false;
		ClippingMode = false;
		break;
	case 4: if (PolygonMode) {
			database.clear();
			triangle_database.clear();
			polygonsEdges.clear();
			numberOfPolygons = 0;
			initializeColors();
			initializeArrays();
			currentP = 0;
		}
		break;
	case 5: 
		database.clear();
		clippedPolygons.clear();
		triangle_database.clear();
		clipped_polygons_triangles.clear();
		exit(1);
	}
}

void subMenuOptionLines(GLint actionOption)
{
	if (!colorForLinesSelected[currentP]) {
		switch (actionOption) {
		case 1: colorForLines[currentP].red = 1.0;
			colorForLines[currentP].green = 0.0;
			colorForLines[currentP].blue = 0.0;
			break;
		case 2: colorForLines[currentP].red = 1.0;
			colorForLines[currentP].green = 1.0;
			colorForLines[currentP].blue = 1.0;
			break;
		case 3:colorForLines[currentP].red = 1.0;
			colorForLines[currentP].green = 1.0;
			colorForLines[currentP].blue = 0.0;
			break;
		case 4:colorForLines[currentP].red = 0.0;
			colorForLines[currentP].green = 0.1;
			colorForLines[currentP].blue = 0.0;
			break;
		case 5:colorForLines[currentP].red = 0.0;
			colorForLines[currentP].green = 1.0;
			colorForLines[currentP].blue = 0.0;
			break;
		case 6:colorForLines[currentP].red = 0.0;
			colorForLines[currentP].green = 0.0;
			colorForLines[currentP].blue = 1.0;
			break;
		case 7:colorForLines[currentP].red = 0.5;
			colorForLines[currentP].green = 1.0;
			colorForLines[currentP].blue = 1.0;
			break;
		case 8:colorForLines[currentP].red = 0.0;
			colorForLines[currentP].green = 0.0;
			colorForLines[currentP].blue = 0.0;
			break;
		case 9:colorForLines[currentP].red = 1.0;
			colorForLines[currentP].green = 0.0;
			colorForLines[currentP].blue = 1.0;
			break;
		case 10:colorForLines[currentP].red = 1.0;
			colorForLines[currentP].green = 0.5;
			colorForLines[currentP].blue = 0.0;
			break;
		case 11:colorForLines[currentP].red = 0.5;
			colorForLines[currentP].green = 0.5;
			colorForLines[currentP].blue = 0.5;
			break;
		case 12:colorForLines[currentP].red = 0.0;
			colorForLines[currentP].green = 0.5;
			colorForLines[currentP].blue = 0.5;
			break;
		case 13:colorForLines[currentP].red = 0.0;
			colorForLines[currentP].green = 0.5;
			colorForLines[currentP].blue = 1.0;
			break;
		case 14:colorForLines[currentP].red = 2.0;
			colorForLines[currentP].green = 0.5;
			colorForLines[currentP].blue = 1.0;
			break;
		case 15:colorForLines[currentP].red = 0.1;
			colorForLines[currentP].green = 0.1;
			colorForLines[currentP].blue = 0.1;
			break;
		case 16:colorForLines[currentP].red = 0.1;
			colorForLines[currentP].green = 0.0;
			colorForLines[currentP].blue = 0.1;
			break;
		}
		colorForLinesSelected[currentP] = true;
	}
	
}

void subMenuOptionFill(GLint actionOption)
{
	if (!colorForFillSelected[currentP]) {
		switch (actionOption) {
		case 1: colorForFill[currentP].red = 1.0;
			colorForFill[currentP].green = 0.0;
			colorForFill[currentP].blue = 0.0;
			break;
		case 2: colorForFill[currentP].red = 1.0;
			colorForFill[currentP].green = 1.0;
			colorForFill[currentP].blue = 1.0;
			break;
		case 3:colorForFill[currentP].red = 1.0;
			colorForFill[currentP].green = 1.0;
			colorForFill[currentP].blue = 0.0;
			break;
		case 4:colorForFill[currentP].red = 0.0;
			colorForFill[currentP].green = 0.1;
			colorForFill[currentP].blue = 0.0;
			break;
		case 5:colorForFill[currentP].red = 0.0;
			colorForFill[currentP].green = 1.0;
			colorForFill[currentP].blue = 0.0;
			break;
		case 6:colorForFill[currentP].red = 0.0;
			colorForFill[currentP].green = 0.0;
			colorForFill[currentP].blue = 1.0;
			break;
		case 7:colorForFill[currentP].red = 0.5;
			colorForFill[currentP].green = 1.0;
			colorForFill[currentP].blue = 1.0;
			break;
		case 8:colorForFill[currentP].red = 0.0;
			colorForFill[currentP].green = 0.0;
			colorForFill[currentP].blue = 0.0;
			break;
		case 9:colorForFill[currentP].red = 1.0;
			colorForFill[currentP].green = 0.0;
			colorForFill[currentP].blue = 1.0;
			break;
		case 10:colorForFill[currentP].red = 1.0;
			colorForFill[currentP].green = 0.5;
			colorForFill[currentP].blue = 0.0;
			break;
		case 11:colorForFill[currentP].red = 0.5;
			colorForFill[currentP].green = 0.5;
			colorForFill[currentP].blue = 0.5;
			break;
		case 12:colorForFill[currentP].red = 0.0;
			colorForFill[currentP].green = 0.5;
			colorForFill[currentP].blue = 0.5;
			break;
		case 13:colorForFill[currentP].red = 0.0;
			colorForFill[currentP].green = 0.5;
			colorForFill[currentP].blue = 1.0;
			break;
		case 14:colorForFill[currentP].red = 2.0;
			colorForFill[currentP].green = 0.5;
			colorForFill[currentP].blue = 1.0;
			break;
		case 15:colorForFill[currentP].red = 0.1;
			colorForFill[currentP].green = 0.1;
			colorForFill[currentP].blue = 0.1;
			break;
		case 16:colorForFill[currentP].red = 0.1;
			colorForFill[currentP].green = 0.0;
			colorForFill[currentP].blue = 0.1;
			break;
		}
		colorForFillSelected[currentP] = true;
	}
}

void glutMenu(void) {
	GLint subMenuColorLines, subMenuColorFill, subMenuAction;
	subMenuColorLines = glutCreateMenu(subMenuOptionLines);
		glutAddMenuEntry("Red", 1);
		glutAddMenuEntry("White", 2);
		glutAddMenuEntry("Yellow", 3);
		glutAddMenuEntry("Forest Green", 4);
		glutAddMenuEntry("Green", 5);
		glutAddMenuEntry("Blue", 6);
		glutAddMenuEntry("Cyan", 7);
		glutAddMenuEntry("Black", 8);
		glutAddMenuEntry("Purple", 9);
		glutAddMenuEntry("Orange", 10);
		glutAddMenuEntry("Violet", 11);
		glutAddMenuEntry("Blue/Green", 12);
		glutAddMenuEntry("Baby Blue", 13);
		glutAddMenuEntry("Lilac", 14);
		glutAddMenuEntry("Dark Grey", 15);
		glutAddMenuEntry("DarkPurple", 16);

	subMenuColorFill = glutCreateMenu(subMenuOptionFill);
		glutAddMenuEntry("Red", 1);
		glutAddMenuEntry("White", 2);
		glutAddMenuEntry("Yellow", 3);
		glutAddMenuEntry("Forest Green", 4);
		glutAddMenuEntry("Green", 5);
		glutAddMenuEntry("Blue", 6);
		glutAddMenuEntry("Cyan", 7);
		glutAddMenuEntry("Black", 8);
		glutAddMenuEntry("Purple", 9);
		glutAddMenuEntry("Orange", 10);
		glutAddMenuEntry("Violet", 11);
		glutAddMenuEntry("Blue/Green", 12);
		glutAddMenuEntry("Baby Blue", 13);
		glutAddMenuEntry("Lilac", 14);
		glutAddMenuEntry("Dark Grey", 15);
		glutAddMenuEntry("DarkPurple", 16);

	subMenuAction = glutCreateMenu(menuOption);
		glutAddMenuEntry("POLYGON", 1);
		glutAddMenuEntry("CLIPPING",2);
		glutAddMenuEntry("EXTRUDE",3);
		glutAddMenuEntry("CLEAR SCHENE",4);
		glutAddMenuEntry("EXIT", 5);

	glutCreateMenu(menuOption);
		glutAddSubMenu("ACTION", subMenuAction);
		glutAddSubMenu("LINE_COLOR", subMenuColorLines);
		glutAddSubMenu("FILL_COLOR", subMenuColorFill);

	glutAttachMenu(GLUT_MIDDLE_BUTTON);
}
//##################################################################################################################################

//############################################-Main Function Repeat over and over again till user exit-#############################
int main(int argc, char** argv)
{
	

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(windowsWidth, windowsHeight);
	glutCreateWindow("Draw Polygons");

	init();
	glutKeyboardFunc(processKey);
	glutSpecialFunc(specialKeys);
	glutDisplayFunc(drawScene);
	glutReshapeFunc(reshape);
	glutIdleFunc(drawScene);

	glutMouseFunc(mouse);
	glutMotionFunc(Motion);
	
	glutMenu();

	glutMainLoop();
}
//##################################################################################################################################
