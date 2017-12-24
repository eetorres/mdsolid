// Edmanuel Torres 
// Read input files
// 2007

#ifndef _CREAD_H_
#define _CREAD_H_

#include<msmvtl/tmatrix.h>
#include<msmvtl/tvmath.h>
#include<mdvect.h>

#include<string>
#include<string.h>

class CRead {
private:

    std::string _file_names, _eam_file, _md_config_file;
    int   _md_steps;
    uint _atomic_number, _basis_atom_number;
    lreal _cut_radius, _shell_thick, _lattice_constant;
    lreal _atomic_mass, _temperature, _time_step;
    //
    TVector<uint>  _cell_number, _cell_xyz;
    TVector<lreal> _shift;
    TMatrix<lreal> _lattice_vectors;
    TMatrix<lreal> _basis_vectors;
    //
    // File functions
    void read_md_name(void);
    void read_eam_name(void);
    //
    void read_atomic_number(void);
    void read_atomic_mass(void);
    void read_lattice_constant(void);
    void read_lattice_vectors(void);
    void read_basis_atom_number(void);
    void read_basis_vectors(void);
    //
    void read_temperature(void);
    void read_time_step(void);
    void read_md_steps(void);
    void read_cut_radius(void);
    void read_shell_thick(void);
    void read_cell_xyz(void);
    void read_shift(void);
    
    // UTILS
    void find_next(std::ifstream&);
    TVector<lreal> ch2vt(std::string);
    
	
public:
	
    // 
    CRead();
    ~CRead();
    
    void read_file_names(void);
    void read_md_config_file(void);
    void summarize_parameters(void);
    //
    std::string get_eam_name(void){ return _eam_file;};
    int   get_md_steps(void){return _md_steps;};
    uint  get_atomic_number(void){ return _atomic_number;};
    uint  get_basis_atom_number(void){return _basis_atom_number;};
    lreal get_atomic_mass(void){return _atomic_mass;};
    lreal get_lattice_constant(void){return _lattice_constant;};
    //
    lreal get_temperature(void){return _temperature;};
    lreal get_time_step(void){return _time_step;};
    lreal get_cut_radius(void){return _cut_radius;};
    lreal get_shell_thick(void){return _shell_thick;};
    
    TVector<uint>  get_cell_number(void){return _cell_number;};
    TVector<lreal> get_shift(void){return _shift;};
    TMatrix<lreal> get_basis_vectors(void){return _basis_vectors;};
    TMatrix<lreal> get_lattice_vectors(void){return _lattice_vectors;};
};


#endif // _CREAD_H_
