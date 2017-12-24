// Edmanuel Torres 
// MD class for molecular dynamics
// 2007-2017

#ifndef _MD_H_
#define _MD_H_


//#include<const.h>
#include<catom.h>
#include<cforce.h>
#include<cread.h>

#include<time.h>
#include<timer.h>

//using namespace std;

const unsigned int _DIM=3;

#ifndef tab
#define tab '\t'
#endif

// Neighbouring cells
const int neighbor_cells[27][3] = {
{  0, 0, 0},
{  1, 0, 0},
{  1, 1, 0},
{  0, 1, 0},
{ -1, 1, 0},
{  0, 0, 1},
{  1, 0, 1},
{ -1, 0, 1}, 
{  1, 1, 1},
{  0, 1, 1},
{ -1, 1, 1},
{ -1,-1, 1}, 
{  0,-1, 1}, 
{  1,-1, 1},
{  0, 0,-1},
{  1, 0,-1},
{ -1, 0,-1}, 
{  1, 1,-1},
{  0, 1,-1},
{ -1, 1,-1},
{ -1,-1,-1}, 
{  0,-1,-1}, 
{  1,-1,-1},
{ -1, 0, 0},
{ -1,-1, 0},
{  0,-1, 0}, 
{  1,-1, 0}
};

class Moldyn {
private:

    string _quim_simb, _error_msg, _struc_type;
    unsigned int _atomic_number;

    unsigned int PARTICLE_NUMBER;
    unsigned int _basis_atom_number, _icell, _NUM_PRTKN, CELL_NUMBER;

    int _md_steps, step_counter, avg_steps;
    int sizeHistVel, countVel, doVelGrph;
	
    bool _if_neighbor;
    
    // Important distances
    lreal _cut_radius, _shell_thick, _lattice_constant;
    lreal _rc, _rc2, _rc_sh, _sh, _sh_2, _rc_sh2;
    // Real units
    lreal _atomic_mass, _temperature;
    // Reduced units
    lreal _sigma, _epsilon, _mu, _tao;
    // Optimization parameters
    lreal _max_vel, _sum_disp, _max_disp;
    // Properties
    lreal _rho;
    lreal k_energy, u_energy;
    lreal sum_k_energy, sum_u_energy, sum_vir;
    lreal _time_step, _real_time;
    // 
    lreal _vv_sum;
    lreal rangeVel, limitVel, _reduced_temp;
    lreal vol, density;
    //
    CForce eam;
    //
    //////////////////////////////////////////////
    TMatrix<lreal> _lattice_vectors, _basis_vectors;
    //
    TVector<lreal> histVel;
    TVector<lreal> _cell_frac, _box_size, _box_middle;
    TVector<lreal> _inv_delta, _positive_r, _shift;
    // old
    TVector<int>  _cell_side, _cell_list, _cell_head;
    TVector<uint> _cell_number;
    TVector<int>  _integer_r, _neighbor_cell, _md_pbc;
    //
    //
    // Control functions
    void adjust_temperature(void);				// Initial atoms velocities
    void set_next_iteration(void);
    // Configuration functions

    //void set_cell_xyz(int,int,int);				// Number of unit cells atoms size
    //void set_cell_xyz(TVector<unsigned int>&);			// Number of unit cells atoms size
    void set_cut_radio(lreal,lreal);
    //
    void set_inverse_cell(void);
    void set_cells(void);
    void set_cell_list(void);
    //
    void set_initial_position(void);				// Initial atomic positions
    void set_initial_velocities(void);				// Initial atomic velocities
    void set_maxwell_velocities(void);				// Initial velocity distribution
    //void set_intial_temperature(lreal);				// Initial temperature

    // Optimization functions
    inline void eval_cell_list(void);				// 
    inline void eval_neighbor(int);
    inline void eval_neighbor_list(void);
    inline void eval_embedded_energy(int);
    inline void eval_embedded_energy(void);
    // Force computation
    void eval_force(void);
    inline void eval_eam_force(int);
    // File utility functions
    void read_main_file(void);
    void read_md_file(void);
    void read_eam_file(string);
    // Reduced units
    void set_reduced_units(void);				// Set lattice constant
    void set_epsilon(lreal);					//
    void set_mass(lreal);					// Set atomic mass
    //
    void set_atomic_number(lreal r){ _atomic_number=r;};
    void set_atomic_mass(lreal r){ _atomic_mass=r;};
    void set_lattice_constant(lreal r){ _lattice_constant=r;};
    void set_basis_atom_number(lreal r){ _basis_atom_number=r;};
    void set_temperature(lreal r){ _temperature=r;};
    //
    void set_time_step(lreal r){ _time_step=r;};
    void set_cut_radius(lreal r){ _cut_radius=r;};
    void set_shell_width(lreal r){ _shell_thick=r;};
    void set_md_steps(int i){ _md_steps=i;};
    //
    void set_lattice_vectors(TMatrix<lreal> v){ _lattice_vectors=v;};
    void set_basis_vectors(TMatrix<lreal> v){ _basis_vectors=v;};
    void set_cells(TVector<uint> v){ _cell_number=v;};
    // read initial configuration file
    void read_md_config_file(void);
    // file utils
    //void find_next(::ifstream&);
    // Vector utility functions
    vector<lreal> VRot(vector<lreal>&,lreal,lreal,lreal);	//
    //dvector	VRand(unsigned int,lreal);			//
    void VRand(unsigned int,lreal,lreal _t[3]);			//
    //dvector	VRandNorm(unsigned int,lreal);			//
    void VRandNorm(unsigned int,lreal,lreal _t[3]);		//
    //TVector<lreal> ch2vt(string);
    //TVector<real>  derivative(TVector<real>&,real);
    // Timer
    Timer timer;
    float sum_time;
    // UTILS
    void summarize_parameters(void);
    // DEBUG UTILS
    void show_message(string);
    void eval_nearest_neighbor(int);
    void eval_nearest_neighbor(void);
	
public:
	
    Particles atoms;
    int is_not_data;
    // 
    Moldyn();
    ~Moldyn();
    
    // Intial conditions
    void set_intial_conditions(void);				// Set initial condition
    
    // Evalue physical properties
    void eval_rho(void);
    //
    void eval_properties(void);
    void eval_u_energy(void);					// Reduced potential atoms energy
    void eval_k_energy(void);					// Reduced kinetic atoms energy
    lreal get_u_energy(void);					// Reduced potential atoms energy
    lreal get_k_energy(void);					// Reduced kinetic atoms energy
    lreal eval_temperature(void);				// lreal temperature of the atoms
    void velocity_distribution(void);
    //
    // Get physical properties
    lreal get_rho(void);
    lreal get_volume(void);
    lreal get_latt_const(void);
    //
    void integration_verlet(void);
    void motion_step(void);
    
    void set_veldist_size(int i){
	sizeHistVel = i;
	histVel.resize(sizeHistVel);
    };
	
    void set_max_veldist(lreal r){ limitVel = r;}
    void set_range_veldist(lreal r){ rangeVel = r;}
    lreal get_veldist_frec(int j){ return histVel[j];}
    lreal get_maxfrec_veldist(void){ return histVel.get_max();}
    int  get_veldist_size(void){ return sizeHistVel;}
	
    // API
    //void 	mdResetAll(void);
    // return the number of epoch in the simulation
    int get_epoch_number(void){ return _md_steps;}
    lreal get_box_size(int i){ return _box_size[i];}
    //dvector	mdGetPartPos(int _i){ return atoms[_i].r;} 		// Return the position atom i in the crystal
    void get_particle_position(lreal*,int);
    // Memory consuming
    int get_max_memalloc(void){ return (atoms.max_size());} 	// This funtion mesure the maxime capacity of the machine
    int get_particle_memalloc(void){return (sizeof(CAtom));} 					// This is the memory size needed for each particle
    int get_total_memalloc(void){ return (atoms.capacity() * get_particle_memalloc());} // This is the size of memory needed for the crystal
    int get_atom_number(void){return PARTICLE_NUMBER;}
    //
    
    void eval_eos_lattice(void);				// Evalue the EOS of energy vs volume
    // Output files
    void write_atom_xyz(string,int);				// Write a file with coordinates of the particles
    void write_atom_pdb(string,int);				// Write a file with the coordinates of the all atoms
    void write_info_file(string);				// Write a file with the simulation conditions
    //
    void show_introduction(void);					// Introduction
    void show_summary(void);
    
};


#endif // _MD_H_
