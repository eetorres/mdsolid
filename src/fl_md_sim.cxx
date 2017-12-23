//========================================================================
// Fl_2D_Gl_Contour V 0.3 (Alpha version)
// Fl_2D_Gl_Contour.h - Fl_2D_Gl_Contour widget class file
// For the Fast Light Tool Kit (FLTK) - www.fltk.org
//
// Copyrigth 2002 by Edmanuel Eduardo Torres Amaris
// Grupo de Fisica y Tecnologia del Plasma
// Physics and Technology Plasma Group
// Universidad Industrial de Santander (UIS) - www.uis.edu.co
// etorresa@tux.uis.edu.co or ukeedman@yahoo.co.uk
// Bucaramanga, Colombia.
//
// Sent me any suggestion, modification or bugs. Don't hesitate to contact
// me for any question, I will be very grateful with your feedbacks.
// Get the last version at www.geocities.com/edmanuelt
//
// This widget is  free  software;  you  can  redistribute  it and/or
// modify it  under  the  terms  of  the   GNU Library General Public
// License  as  published  by  the  Free  Software Foundation; either
// version 2 of the License,  or  (at your option)  any later version.
//
// This  library  is  distributed  in the hope that it will be useful,
// but  WITHOUT ANY WARRANTY;  without  even  the  implied warranty of
// MERCHANTABILITY  or FITNESS FOR A PARTICULAR PURPOSE.   See the GNU
// Library General Public License for more details.
//
// You should have  received a copy  of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA.
//
//========================================================================
#include <fl_md_sim.h>

#define ALPHA 0.5
#define HAVE_GL 1
//#define _SHOW_LINES_ 1
//#define _SHOW_MESSAGES 1

#if HAVE_GL
fl_md_sim::fl_md_sim(int x,int y,int w,int h,const char *l) : Fl_Gl_Window(x,y,w,h,l)
#else
fl_md_sim::fl_md_sim(int x,int y,int w,int h,const char *l) : Fl_Box(x,y,w,h,l)
#endif /* HAVE_GL */
{
	vAng = 250.0;
	hAng = 18.0;
	zoom=1.0;

	size = 10.0;
	xshift=0.0;
	yshift=0.0;
	zshift=0.0;
	y_off=1.0;
	
	gph_ion = 1;

	// Foreground colors
	fcred = 0.5;
	fcgreen = 1.0;
	fcblue = 0.5;
	// Background colors
	bgred	= 0.0;
	bggreen	= 0.0;
	bgblue	= 0.0;

#if !HAVE_GL
    label("OpenGL is required for this demo to operate.");
    align(FL_ALIGN_WRAP | FL_ALIGN_INSIDE);
#endif /* !HAVE_GL */
}

////////////////////////////OPENGL////////////////////////////////

#if HAVE_GL

void fl_md_sim::drawPart(void){
    //gm_rgb _rgb;
    //dvector prt;
    lreal prt[3];
    int _prt;
    glPointSize(2.0);
    glRotatef(hAng,0,1,0); 
    glRotatef(vAng,1,0,0);
    glScalef(zoom,zoom,zoom);
    
    real ibndx = 1.0/get_box_size(0);
    real ibndy = 1.0/get_box_size(1);
    real ibndz = 1.0/get_box_size(2);
    
    drawBox();
    if(gph_ion){
	_prt = get_atom_number();
	glBegin(GL_POINTS);
	for(int l=0; l<_prt; l++){
	    //prt = mdGetPartPos(l);
	    get_particle_position(prt,l);
	    //_icc = mdcubepos(prt);
	    glColor3f(0.7,0.7,0.7);
	    //glColor3f(0.1,1.0,0.0);
	    //_rgb = cmHsv(dvertex[l].z());
	    //glColor3f(_rgb.r,_rgb.g,_rgb.b);
	    //glVertex3f(0,0,0);
	    glVertex3f(prt[0]*ibndx,prt[1]*ibndy,prt[2]*ibndz);
	}
	glEnd();
    }
    /*
    _prt = mdGetNumAtomIons();
    glBegin(GL_POINTS);
    for(int l=0; l<_prt; l++){
	if(ions[l]._sw){
	    prt = mdGetIonPos(l);
	    ////_icc = mdcubepos(prt);
	    glColor3f(1.0,0.3,0.0);
	    //glColor3f(0.1,0.1,0.1);
	    ////_rgb = cmHsv(dvertex[l].z());
	    ////glColor3f(_rgb.r,_rgb.g,_rgb.b);
	    ////glVertex3f(0,0,0);
	    glVertex3f(prt[0]*ibndx,prt[1]*ibndy,prt[2]*ibndz);
	}
    }
    glEnd();
	*/
}

void fl_md_sim::drawBox(void){
    glBegin(GL_LINES);
    //glColor3f(1.0,1.0,1.0);
    //glColor3f(0.0,0.0,0.0);
    glColor3f(1.0,0.0,0.0);
    glVertex3f(-0.5,-0.5,-0.5);
    glVertex3f( 0.5,-0.5,-0.5);
    glColor3f(0.0,1.0,0.0);
    glVertex3f(-0.5,-0.5,-0.5);
    glVertex3f(-0.5, 0.5,-0.5);
    glColor3f(0.0,0.0,1.0);
    glVertex3f(-0.5,-0.5,-0.5);
    glVertex3f(-0.5,-0.5, 0.5);
    glEnd();
}

void fl_md_sim::draw() {
    if (!valid()) {
	glClearColor(bgred,bggreen,bgblue,1.0);
        glLoadIdentity();
        glViewport(0,0,w(),h());
        glOrtho(-1,1,-1,1,-2,2);
	//glOrtho(-50,50,-50,50,-1000,1000);
	//glOrtho( 8, 9, 8, 9, -10, 20);
    }
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();
    
    
    if(is_not_data)
	drawPart();

    glPopMatrix();
    glFlush();
    redraw();
}



#endif /* HAVE_GL */

////////////////////////////HANDLE EVENTS////////////////////////////////

int fl_md_sim::handle(int event) {
    static int last_x;
    static int last_y;
    int delta_x, delta_y;
	//... position in Fl::event_x() and Fl::event_y()
	// get current mouse position and process event
    int x = Fl::event_x();
    int y = Fl::event_y();
	
	switch(event) {
	case FL_PUSH:
		//... mouse down event ...
		// save mouse position to track drag events
		last_x = x;
		last_y = y;
		return 1;
	case FL_DRAG:
		delta_x = x - last_x;
		delta_y = y - last_y;
		last_x = x;
		last_y = y;
		hAng += 0.2*delta_x;
		vAng += 0.2*delta_y;
		redraw();
		return 1;
		/*case FL_RELEASE:   
		... mouse up event ...
		return 1;
		case FL_FOCUS :
		case FL_UNFOCUS :
		... Return 1 if you want keyboard events, 0 otherwise
		return 1;
		case FL_KEYBOARD:
		... keypress, key is in Fl::event_key(), ascii in Fl::event_text()
		... Return 1 if you understand/use the keyboard event, 0 otherwise...
		return 1;
		case FL_SHORTCUT:
		... shortcut, key is in Fl::event_key(), ascii in Fl::event_text()
		... Return 1 if you understand/use the shortcut event, 0 otherwise...
		return 1;*/
  default:
    // pass other events to the base class...
    return Fl_Gl_Window::handle(event);
  }
}

//////////////////////////////UTILS///////////////////////////////

void fl_md_sim::bgColor(float r,float g,float b){
	bgred=r;
	bggreen=g;
	bgblue=b;
}

void fl_md_sim::fgColor(float r,float g,float b){
	fcred=r;
	fcgreen=g;
	fcblue=b;
}

void fl_md_sim::clear_all(void){
	//clear_cl3d();
	//clear_cn3d();
}
/*
void fl_md_sim::graph_cb(void){
#if _SHOW_MESSAGES
#endif // !HAVE_GL
}
*/
/////
