// generated by Fast Light User Interface Designer (fluid) version 1.0305

#ifndef mdsolid_h
#define mdsolid_h
#include <FL/Fl.H>
#include <cmd.h>
extern int _stp, _n_ep; 
#include <FL/Fl_Double_Window.H>
extern Fl_Double_Window *w_implant;
#include <fl_md_sim.h>
extern fl_md_sim *viz3d;
#include <FL/Fl_Group.H>
#include <FL/Fl_Button.H>
extern Fl_Button *init;
extern Fl_Button *play;
extern Fl_Button *stop;
#include <FL/Fl_Slider.H>
extern Fl_Slider *_z;
extern Fl_Button *_x_f;
extern Fl_Button *_y_f;
extern Fl_Button *_z_f;
#include <FL/Fl_Box.H>
extern Fl_Button *eos;
extern Fl_Group *_u_gph;
#include<Cartesian.H>
extern Ca_Canvas *_u_ener;
extern Ca_X_Axis *_u_t;
extern Ca_Y_Axis *_u_e;
extern Fl_Group *_k_gph;
extern Ca_Canvas *_k_ener;
extern Ca_X_Axis *_k_t;
extern Ca_Y_Axis *_k_e;
extern Fl_Group *_t_gph;
extern Ca_Canvas *_t_ener;
extern Ca_X_Axis *_t_t;
extern Ca_Y_Axis *_t_e;
extern Fl_Group *_max_gph;
extern Ca_Canvas *_m_vel;
extern Ca_X_Axis *_m_v;
extern Ca_Y_Axis *_m_f;
void set_k_axis();
void set_u_axis();
void set_t_axis();
void set_m_axis();
#endif