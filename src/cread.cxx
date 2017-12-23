// Programado por: Edmanuel Torres Amaris
// simulaciones de Dinamica Molecular de Solidos

#include <cread.h>

// READ CONFIGURATION FILES

#define __DEBUG_MESSAGES__

CRead::CRead(void){
    // configuration file name
    _file_names = "files.in";
    _eam_file  = "";
    _md_config_file = "";
    _lattice_vectors.resize(3,3);
    _cell_xyz.resize(3);
}


void CRead::read_file_names(void){
    read_eam_name();
    read_md_name();
}

// READ INITIAL CONDITIONS

void CRead::read_md_config_file(void){
    std::cout<<"[read atomic number]"<<std::endl;
    read_atomic_number();
    std::cout<<"[read atomic mass]"<<std::endl;
    read_atomic_mass();
    std::cout<<"[read lattice constant]"<<std::endl;
    read_lattice_constant();
    std::cout<<"[read lattice vectors]"<<std::endl;
    read_lattice_vectors();
    std::cout<<"[read atom number]"<<std::endl;
    read_basis_atom_number();
    std::cout<<"[read basis vectors]"<<std::endl;
    read_basis_vectors();
    std::cout<<"[read shitf]"<<std::endl;
    read_shift();
    //
    std::cout<<"[read temperature]"<<std::endl;
    read_temperature();
    std::cout<<"[read time step]"<<std::endl;
    read_time_step();
    std::cout<<"[read md steps]"<<std::endl;
    read_md_steps();
    std::cout<<"[read cut radius]"<<std::endl;
    read_cut_radius();
    std::cout<<"[read shell thick]"<<std::endl;
    read_shell_thick();
    std::cout<<"[read cell]"<<std::endl;
    read_cell_xyz();
}

void CRead::summarize_parameters(void){
    std::cout<<"===============SUMMARY================="<<std::endl;
    std::cout<<" Atomic number = "<<_atomic_number<<std::endl;
    std::cout<<" Atomic mass = "<<_atomic_mass<<std::endl;
    std::cout<<" Lattice constant= "<<_lattice_constant<<std::endl;
    std::cout<<" Lattice_vectors : "<<_lattice_vectors<<std::endl;
    std::cout<<" Unit cell atoms = "<<_basis_atom_number<<std::endl;
    std::cout<<" Basis atoms : "<<_basis_vectors<<std::endl;
    std::cout<<" Shift : "<<_shift<<std::endl;
    std::cout<<"===============SUMMARY================="<<std::endl;
}

// READ FILE NAMES

void CRead::read_eam_name(){
    char _chr_buff[1024], _ch, tag[64], _name[64];
    std::ifstream cfg(_file_names.c_str());
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    //puts(_buff);
	    if( !strncmp(_chr_buff,"EAMP",4)){
		sscanf(_chr_buff,"%s %s",tag,_name);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s : %s\n",tag,_name);
#endif
		_eam_file  = _name;
		break;
	    }
	}
	find_next(cfg);
    }
    cfg.close();
}

void CRead::read_md_name(){
    char _chr_buff[1024], _ch, tag[64], _name[64];
    std::ifstream cfg(_file_names.c_str());
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    //puts(_buff);
	    if( !strncmp(_chr_buff,"NAME",4)){
		sscanf(_chr_buff,"%s %s",tag,_name);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s : %s\n",tag,_name);
#endif
		_md_config_file = _name;
		break;
	    }
	}
	find_next(cfg);
    }
    cfg.close();
}

// READ INTIAL PARAMETERS

void CRead::read_atomic_number(void){
    char _chr_buff[1024], _ch, tag[64];
    unsigned int val;
    std::ifstream cfg(_md_config_file.c_str());
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    //puts(_buff);
	    if( !strncmp(_chr_buff,"ZNUM",4)){
		sscanf(_chr_buff,"%s %u",tag,&val);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s : %u\n",tag,val);
#endif
		_atomic_number = val;
		break;
		//NUM_LAYERS = val;
		//layers.resize(NUM_LAYERS);
	    }
	}
	find_next(cfg);
    }
    cfg.close();
}

void CRead::read_atomic_mass(void){
    char _chr_buff[1024], _ch, tag[64];
    lreal val;
    std::ifstream cfg(_md_config_file.c_str());
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    if( !strncmp(_chr_buff,"MNUM",4)){
		sscanf(_chr_buff,"%s %Lg",tag,&val);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s : %Lg\n",tag,val);
#endif
		_atomic_mass = val;
		break;
	    }
	}
    	find_next(cfg);
    }
    cfg.close();
}

void CRead::read_lattice_constant(void){
    char _chr_buff[1024], _ch, tag[64];
    lreal val;
    std::ifstream cfg(_md_config_file.c_str());
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    if( !strncmp(_chr_buff,"LATT",4)){
		sscanf(_chr_buff,"%s %Lg",tag,&val);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s : %Lg\n",tag,val);
#endif
		_lattice_constant = val;
		break;
	    }
	}
    	find_next(cfg);
    }
    cfg.close();
    //printf("...bye bye...");
}

void CRead::read_basis_atom_number(void){
    char _chr_buff[1024], _ch, tag[64];
    unsigned int val;
    std::ifstream cfg(_md_config_file.c_str());
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    if( !strncmp(_chr_buff,"NATM",4)){
		sscanf(_chr_buff,"%s %u",tag,&val);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s : %u\n",tag,val);
#endif
		_basis_atom_number = val;
		break;
	    }
	}
    	find_next(cfg);
    }
    cfg.close();
}

void CRead::read_temperature(void){
    char _chr_buff[1024], _ch, tag[64];
    lreal val;
    std::ifstream cfg(_md_config_file.c_str());
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    if( !strncmp(_chr_buff,"TEMP",4)){
		sscanf(_chr_buff,"%s %Lg",tag,&val);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s : %Lg\n",tag,val);
#endif
		_temperature = val;
		break;
	    }
	}
    	find_next(cfg);
    }
    cfg.close();
}

void CRead::read_md_steps(void){
    char _chr_buff[1024], _ch, tag[64];
    int val;
    std::ifstream cfg(_md_config_file.c_str());
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    if( !strncmp(_chr_buff,"NSTP",4)){
		sscanf(_chr_buff,"%s %i",tag,&val);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s : %i\n",tag,val);
#endif
		_md_steps = val;
		break;
	    }
	}
    	find_next(cfg);
    }
    cfg.close();
}

void CRead::read_time_step(void){
    char _chr_buff[1024], _ch, tag[64];
    lreal val;
    std::ifstream cfg(_md_config_file.c_str());
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    if( !strncmp(_chr_buff,"TIME",4)){
		sscanf(_chr_buff,"%s %Lg",tag,&val);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s : %Lg\n",tag,val);
#endif
		_time_step = val;
		break;
	    }
	}
    	find_next(cfg);
    }
    cfg.close();
}

void CRead::read_cut_radius(void){
    char _chr_buff[1024], _ch, tag[64];
    lreal val;
    std::ifstream cfg(_md_config_file.c_str());
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    if( !strncmp(_chr_buff,"CUTR",4)){
		sscanf(_chr_buff,"%s %Lg",tag,&val);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s : %Lg\n",tag,val);
#endif
		_cut_radius = val;
		break;
	    }
	}
    	find_next(cfg);
    }
    cfg.close();
}

void CRead::read_shell_thick(void){
    char _chr_buff[1024], _ch, tag[64];
    lreal val;
    std::ifstream cfg(_md_config_file.c_str());
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    if( !strncmp(_chr_buff,"SHLL",4)){
		sscanf(_chr_buff,"%s %Lg",tag,&val);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s : %Lg\n",tag,val);
#endif
		_shell_thick = val;
		break;
	    }
	}
    	find_next(cfg);
    }
    cfg.close();
}

void CRead::read_lattice_vectors(void){
    char _chr_buff[1024], _ch, tag[64];
    std::ifstream cfg(_md_config_file.c_str());
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    if( !strncmp(_chr_buff,"LVEC",4)){
		sscanf(_chr_buff,"%s",tag);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s : \n",tag);
#endif
		cfg.getline(_chr_buff,256);
		_lattice_vectors[0] = ch2vt((std::string)_chr_buff);
		cfg.getline(_chr_buff,256);
		_lattice_vectors[1] = ch2vt((std::string)_chr_buff);
		cfg.getline(_chr_buff,256);
		_lattice_vectors[2] = ch2vt((std::string)_chr_buff);
#ifdef __DEBUG_MESSAGES__
		std::cout<<"_lattice_vectors = "<<_lattice_vectors<<std::endl;
#endif
		break;
	    }
	}
    	find_next(cfg);
    }
    cfg.close();
}

void CRead::read_shift(void){
    char _chr_buff[1024], _ch, tag[64];
    std::ifstream cfg(_md_config_file.c_str());
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    if( !strncmp(_chr_buff,"SHFT",4)){
		sscanf(_chr_buff,"%s",tag);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s : \n",tag);
#endif
		cfg.getline(_chr_buff,256);
		_shift = ch2vt((std::string)_chr_buff);
		break;
	    }
	}
    	find_next(cfg);
    }
    cfg.close();
}

void CRead::read_basis_vectors(void){
    char _chr_buff[1024], _ch, tag[64];
    std::ifstream cfg(_md_config_file.c_str());
    _basis_vectors.resize(_basis_atom_number,3);
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    if( !strncmp(_chr_buff,"APOS",4)){
		sscanf(_chr_buff,"%s",tag);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s : \n",tag);
#endif
		for(unsigned int _i=0; _i<_basis_atom_number; _i++){
		    cfg.getline(_chr_buff,256);
		    _basis_vectors[_i] = ch2vt((std::string)_chr_buff);
	        }
		break;
	    }
	}
    	find_next(cfg);
    }
    cfg.close();
}

void CRead::read_cell_xyz(void){
    TVector<unsigned int> _t(3);
    char _chr_buff[1024], _ch, tag[64];
    std::ifstream cfg(_md_config_file.c_str());
    while( !cfg.eof() ){
	_ch = cfg.peek();
	if(_ch != '#' && _ch != '\n'){
	    cfg.getline(_chr_buff,1024);
	    if( !strncmp(_chr_buff,"CLLS",4)){
		sscanf(_chr_buff,"%s%u%u%u",tag,&_t[0],&_t[1],&_t[2]);
#ifdef __DEBUG_MESSAGES__
		printf("Tag: %s = %ux%ux%u\n",tag,_t[0],_t[1],_t[2]);
#endif
		_cell_number=_t;
		break;
	    }
	}
    	find_next(cfg);
    }
    cfg.close();
}

// UTILITY FUNCTIONS

TVector<lreal> CRead::ch2vt(std::string _ch){
    TVector<lreal> _t(3);
    //lreal _r1, _r2, _r3;
#ifdef __DEBUG_MESSAGES__
    std::cout<<"ch2vt: "<<_ch<<std::endl;
#endif
    sscanf(_ch.c_str(),"%Lf %Lf %Lf",&_t[0],&_t[1],&_t[2]);
    //_t
    return _t;
}

void CRead::find_next(std::ifstream& f){
    char _ch;// = f.peek();
    _ch = f.peek();
    //printf("[%c]\n",_ch);
    if(_ch == '#' || _ch==' ' || _ch=='\n'){
	f.ignore(1024,'\n');
	_ch = f.peek();
    }
    while(_ch == '#' || _ch==' ' || _ch=='\n'){
        f.ignore(1024,'\n');
        _ch = f.peek();
	//printf("[%c]\n",_ch);
    }
}

// END

