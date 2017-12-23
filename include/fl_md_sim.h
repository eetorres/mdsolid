//========================================================================
// Fl_2D_Gl_Contour V 0.3 (Alpha version)
// Fl_2D_Gl_Contour.h - Fl_2D_Gl_Contour widget class file
// For the Fast Light Tool Kit (FLTK) - www.fltk.org
//
// Copyrigth 2002 by Edmanuel Eduardo Torres Amaris
// Grupo de Fisica y Tecnologia del Plasma
// Physics and Technology Plasma Group
// Universidad Industrial de Santander (UIS) - www.uis.edu.co
// etorres@tux.uis.edu.co or ukeedman@yahoo.co.uk
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
#ifndef _fl_md_sim_H_
#define _fl_md_sim_H_


//#include <msmvtl/tmatrix.h>
//#include <msmvtl/tvmath.h>
//#include <msmvtl/tmmath.h>
//#include <Fl/config.h>

#include <cmd.h>
#include <timer.h>

#include <FL/Fl.H>

#define HAVE_GL 1

#if HAVE_GL
    #include <FL/Fl_Gl_Window.H>
    #include <FL/gl.h>
#else
    #include <FL/Fl_Box.H>
#endif /* HAVE_GL */

#include <stdlib.h>

typedef struct{
    double r;
    double g;
    double b;
} gm_rgb;

class fl_md_sim : public Fl_Gl_Window, public Moldyn {

public:

    //int is_not_data;
    double size;
    //Fl_Map_Contour<gm_real, 2> vdata;
    //gm_rgb oldrgb, rgb;
	
    fl_md_sim();
    fl_md_sim(int,int,int,int,const char* l=0);

    //void contourdata(TMatrix<gm_real>&);
    //void initData(void);
    void graph_cb(void);
    void drawPart(void);
    void drawBox(void);
//    void redraw(void){draw();};
    void bgColor(float r,float g,float b);
    void fgColor(float r,float g,float b);
    // Scale function
    void scale(float _s){_scl = _s;};
    // the rotation about the vertical (y ) axis
    void v_angle(float angle){vAng=angle;};
    float v_angle(){return vAng;};
    // the rotation about the horizontal (x ) axis
    void h_angle(float angle){hAng=angle;};
    float h_angle(){return hAng;};
    
    void szoom(float z){zoom=z;};
    float szoom(){return zoom;};
    
    void panx(float x){xshift=x;};
    void pany(float y){yshift=y;};
    
    void graph_ions(int _gi){gph_ion=_gi;};
    
#if HAVE_GL
#else
#endif /* HAVE_GL */
#if HAVE_GL
    void draw();    
#endif /* HAVE_GL */

protected:
    int handle(int);

private:

    int   _draw_mpts, _draw_data, gph_ion;                  //
    float vAng,hAng, y_off, _scl;
    float xshift,yshift,zshift;
    float zoom;
    float fcred, fcgreen, fcblue;
    float bgred, bggreen, bgblue;
    void  clear_all(void);
    ////
};

#endif

