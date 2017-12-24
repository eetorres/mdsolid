// Edmanuel Torres
// Molecular dynamics simulation

#include <cmd.h>
#include <msmvtl/const.h>

// debugging settings
//#define __DEBUG_FUNCTIONS__
//#define __DEBUG_MESSAGES__

//#define _MD_DEBUG_FORCES_	1
//#define _MD_DEBUG_VERLET_	1
//#define _MD_DEBUG_NEIGHBOR_ 	1
//#define _MD_MSG_ 		1
//#define _MD_DPL_INF_ 		1
// Utils
#define _MD_TIMER_		1
// Configuration
//#define _USE_TIME_STEP_ 	1
#define _USE_SHELL_LIST_	1
// Configuration

// Constructor for the Molecular Dynamics simulation
Moldyn::Moldyn(){
    atoms.clear();
    //
    u_energy=0;
    k_energy=0;
    sum_k_energy = 0;
    sum_k_energy = 0;
    step_counter = 0;
    // particle density
    _rho = 0;
    // the default simulation time-step
    _time_step = 1E-3;
    //
    _cell_side.resize(3);
    _cell_frac.resize(3);
    _cell_number.resize(3);
    //mdkn.resize(3);
    _positive_r.resize(3);
    _shift.resize(3);
    _shift[0]=0.25;
    _shift[1]=0.25;
    _shift[2]=0.25;
    _integer_r.resize(3);
    _neighbor_cell.resize(3);
    _md_pbc.resize(3);
    _lattice_vectors.resize(3,3);
    //
    _box_size.resize(3);
    _box_middle.resize(3);
    //
    is_not_data = 0;
    //_show_inf 	= 0;
    //_show_time_step = 0;
    //_show_ion_k_ener= 0;
    // displacement acumulator
    
    ///////////////////////
    doVelGrph = 0;
    rangeVel = 50.0;
    limitVel = 50;
    countVel=0;
    sizeHistVel = 50;
    histVel.resize(sizeHistVel);
    _if_neighbor = true;
    avg_steps = 100;
    is_not_data = 0;
    // Maximun velocity
    _max_vel = 50.0;
    // Maximun displacement
    _max_disp = 0.1;
    _sum_disp = 0;
    sum_time = 0;
    // read the initial configuration file
    read_md_config_file();
}

void Moldyn::set_reduced_units(void){
#ifdef __DEBUG_FUNCTIONS__
  show_message("set_reduced_units");
#endif
  _epsilon = K_B;
  _sigma   = A;
  _mu      = _atomic_mass*AMU;
  _tao     = _sigma*sqrt(_mu/_epsilon);
  //
  vol = (_sigma * _sigma * _sigma);
  density = _rho * (_mu / vol);
#ifdef __DEBUG_FUNCTIONS__
  show_message("set_reduced_units");
#endif
  _reduced_temp = sqrt(2.0*_temperature);
  printf("Reduced temperature = %Le\n",_reduced_temp);
}

// Destructor - release the used memory during the simulation
Moldyn::~Moldyn(){
    atoms.clear();
}

/////////////////////// MD CORE FUNCTIONS //////////////////////////////////////

// Linked and shell cell configuration functions
void Moldyn::set_cells(void){ 
    _cell_side = iVScale(_box_size, (1.0/_rc_sh));
}

void Moldyn::set_inverse_cell(void){
    _cell_frac = lVDiv(_cell_side,_box_size);
}

void Moldyn::set_cell_list(void){
    CELL_NUMBER = (int)vVol(_cell_side);
    _cell_head.resize(CELL_NUMBER);
    _cell_list.resize(PARTICLE_NUMBER);
}

void Moldyn::eval_cell_list(void){
#ifdef __DEBUG_FUNCTIONS__
  show_message("eval_cell_list");
#endif
    unsigned int _n;
    for(_n=0; _n<CELL_NUMBER; _n++)
        _cell_head[_n] = -1;
    for(_n=0; _n<PARTICLE_NUMBER; _n++){
        _positive_r = lVAdd(atoms[_n].r,_box_middle);
        _integer_r = iVMul(_positive_r,_cell_frac);
    	_icell = iVLinear(_integer_r,_cell_side);
	_cell_list[_n] =  _cell_head[_icell];
	_cell_head[_icell] = _n;
    }
#ifdef __DEBUG_FUNCTIONS__
  show_message("eval_cell_list");
#endif
}

/// Initial conditions

// Set cut radius of particles interaction
void Moldyn::set_cut_radio(lreal rc, lreal rs){
    _rc 	= rc;							// Cut radius
    _rc2 	= _rc*_rc;						// Square cut radius
    _sh 	= rs;							// Shell thinkness
    _sh_2 	= 0.5*_sh;						// Half shell thinkness
    _rc_sh	= _rc+_sh;						// Shell cut radius
    _rc_sh2 	= _rc_sh*_rc_sh;					// Square shell cut radius
}

// The initial positons for the atoms in the system
void Moldyn::set_initial_position(void){
#ifdef __DEBUG_FUNCTIONS__
  show_message("set_initial_position");
#endif
    unsigned int z, x, y, nbase, i, _dim;
    lreal rand_shift, _xyz[3];
    CAtom part;
    atoms.clear();
    // Initial parameters
    part.n_i=0;
    //part.pot_par=0;
    //part.pot_emb=0;
    part._sw=true;
    for(i=0; i<_DIM; i++){
	part.v[i] = 0; //((double)rand()/RAND_MAX)-0.5;
	part.a[i]=0;
	part._a[i]=0;
    }
    srand((unsigned)time(0));
    for(z=0; z<_cell_number[2]; z++){
	_xyz[2]=(lreal)z;
	for(x=0; x<_cell_number[0]; x++){
	    _xyz[0]=(lreal)x;
	    for(y=0; y<_cell_number[1]; y++){
		_xyz[1]=(lreal)y;
		for(nbase=0; nbase<_basis_atom_number; nbase++){
		    for(_dim=0;_dim<_DIM;_dim++){
			rand_shift = 0.1*(((double)rand()/RAND_MAX)-0.5);
			// Lattice vectors have not been implemented, basis vector must be given with respect the lattice constant instead
			part.r[_dim]=((double)_xyz[_dim]+_shift[_dim]+_basis_vectors[nbase][_dim]+rand_shift)*_lattice_constant*_lattice_vectors[_dim][_dim];
			//part.r[_dim]=((double)_xyz[_dim]+_shift[_dim]+_basis_vectors[nbase][_dim]+rand_shift)*_lattice_constant;
			while(part.r[_dim] < -_box_middle[_dim])					// Periodic boundary conditions
			    part.r[_dim] += _box_size[_dim];
			while(part.r[_dim] > _box_middle[_dim])					// Periodic boundary conditions
			    part.r[_dim] -= _box_size[_dim];
		    }
		    atoms.push_back(part);
		}
	    }
	}
    }
    PARTICLE_NUMBER = atoms.size();
#ifdef __DEBUG_FUNCTIONS__
  show_message("set_initial_position");
#endif
}

// Working
// updated 01092011: not necesary this is done in the reduced units
//void Moldyn::set_intial_temperature(lreal _temp){
//    _reduced_temp = sqrt(2.0*_temp);
//    printf("Reduced temperature = %Le\n",_reduced_temp);
//}

// The initial velocities for the atoms in the crystal Iron (bcc) atoms
void Moldyn::set_initial_velocities(void){
    unsigned int i;
    _vv_sum = 0;
    srand((unsigned int)time(NULL));
    for(i=0; i<PARTICLE_NUMBER; i++){
	VRandNorm(3,_max_vel,atoms[i].v);
	//atoms[i].v = (0.30 * atoms[i].v);
    }
    //adjust_temperature();
    eval_properties();
    printf("E_K = %Lf\n",0.5*_vv_sum);
}

// The initial velocities for the atoms in the crystal Iron (bcc) atoms
void Moldyn::set_maxwell_velocities(void){
#ifdef __DEBUG_FUNCTIONS__
  show_message("set_maxwell_velocities");
#endif
    unsigned int i;
    lreal vSumx = 0, vSumy = 0, vSumz = 0;
    srand((unsigned int)time(NULL));
    for(i=0; i<PARTICLE_NUMBER; i++){
	//atoms[i].v=VRand(3,_reduced_temp);
	VRand(3,_reduced_temp,atoms[i].v);
	vSumx += atoms[i].v[0];
	vSumy += atoms[i].v[1];
	vSumz += atoms[i].v[2];
    }
    vSumx /= PARTICLE_NUMBER;
    vSumy /= PARTICLE_NUMBER;
    vSumz /= PARTICLE_NUMBER;
    for(i=0; i<PARTICLE_NUMBER; i++){
	atoms[i].v[0]-=vSumx;
	atoms[i].v[1]-=vSumy;
	atoms[i].v[2]-=vSumz;
    }
#ifdef __DEBUG_FUNCTIONS__
  show_message("set_maxwell_velocities");
#endif
}

// Working
void Moldyn::adjust_temperature(void){
    unsigned int _i, _j;
    lreal _vFac, vv, vvSum=0;
    lreal vSumx=0, vSumy=0, vSumz=0;
    for(_i=0; _i<PARTICLE_NUMBER; _i++){
        vv = 0;
        for(_j=0; _j<_DIM; _j++)
	    vv += atoms[_i].v[_j]*atoms[_i].v[_j];
	vvSum += vv;
    }
    if(vvSum>0){
	_vFac = _reduced_temp/sqrt(vvSum/PARTICLE_NUMBER);
	for(_i=0; _i<PARTICLE_NUMBER; _i++)
	    for(_j=0; _j<_DIM; _j++)
		atoms[_i].v[_j]*=_vFac;
	for(_i=0; _i<PARTICLE_NUMBER; _i++){
	    vSumx += atoms[_i].v[0];
	    vSumy += atoms[_i].v[1];
    	    vSumz += atoms[_i].v[2];
	}
	vSumx /= PARTICLE_NUMBER;
	vSumy /= PARTICLE_NUMBER;
	vSumz /= PARTICLE_NUMBER;
	for(_i=0; _i<PARTICLE_NUMBER; _i++){
	    atoms[_i].v[0]-=vSumx;
	    atoms[_i].v[0]-=vSumy;
	    atoms[_i].v[0]-=vSumz;
	}
    }
}

// Periodic boundary conditions and initial postion
void Moldyn::set_intial_conditions(void){
#ifdef __DEBUG_FUNCTIONS__
  show_message("set_intial_conditions");
#endif
    // set reduced units
    set_reduced_units();
    // start the initial positions for the particle in the system
    set_initial_position();
    // set the initial velocities for the particle in the system
    //set_initial_velocities();
    //1e-10*sqrt((mdGetTargetAtomMass()*1.6605e-27)/1.3806503e-23);
    // The lreal time step scale
    //_real_time = _tao; 
    // number of cells in each direction
    set_cells();
    // reciprocal of the cells in each direction
    set_inverse_cell();
    // vector allocation for the cell list vectors
    set_cell_list();
    // allow particle graph
    is_not_data = 1;
    // reset particles acelerations
    set_next_iteration();
    // neigbohr list evaluation
    eval_neighbor_list();
    // force evaluation
    eval_embedded_energy();
    eval_force();
    //embedded_energy();
    // summary
    //summarize_parameters();
    //
    eval_u_energy();
    eval_k_energy();
    //
    printf("[K0]: %0.3Lf\t[U0]: %0.3Lf eV\t[E]: %0.3Lf\n",k_energy,u_energy,k_energy+u_energy);
#ifdef __DEBUG_FUNCTIONS__
  show_message("set_intial_conditions");
#endif
}

/////////////////////// OPTIMIZATION NL FUNCTIONS //////////////////////////////

// Gives number of atoms closer than sqrt(r2max) from atom number i
void Moldyn::eval_neighbor(int i){
    lreal dx, r2, r, u_val, rho;
    int j;
    //TVector<lreal> delta(3);
#ifdef _MD_DEBUG_NEIGHBOR_
  show_message("End mdTargetTargetNeighbors");
  cout<<" ("<<i<<")"<<endl;
#endif

    atoms[i].neighbor.clear();
    _positive_r = lVAdd(atoms[i].r,_box_middle);
    _integer_r = iVMul(_positive_r,_cell_frac);
    //  Here, computing the neighbour list
    for (int _m=0; _m<14; _m++){
        _neighbor_cell = iVAdd(_integer_r,neighbor_cells[_m]);					// Neighbour cell
											// Free _DIM here with PCB
        for (unsigned int coord=0; coord<_DIM; coord++){				// For each _DIM
	    //if(_md_pbc[coord])							// Ask if PBC is used
    	    if(_neighbor_cell[coord] >= _cell_side[coord]){					// Check if  PBC is necessary
		_neighbor_cell[coord]= 0;							// Apply PBC to each  cell
    		//delta[coord]= _box_size[coord];						// Apply PBC to each particle
	    }else if(_neighbor_cell[coord] < 0){						// Check if  PBC is necessary
		_neighbor_cell[coord]= _cell_side[coord]-1;					// Apply PBC to each cell
		//delta[coord]=-_box_size[coord];						// Apply PBC to each particle
	    }//else{delta[coord]=0;}							// Leave it as it...!
	}
	_icell = iVLinear(_neighbor_cell,_cell_side);						// Look for the cell
	//if(_icell>=0 && _icell < CELL_NUMBER)						// Inside of the atoms
	j = _cell_head[_icell];								// Cell position
	//else j = -1;									// Out side of cell
	while(1){									// Over all the particles in the cell
	    if(j<0) break;								// Get out of the cell
	    if(_m!=0 || j>i){								// Avoid auto-interaction 
		r2 = 0;									// Set distance to cero
		for (unsigned int coord=0; coord<_DIM; coord++){					// Each _DIM
		    dx = atoms[j].r[coord] - atoms[i].r[coord];// + delta[coord];	// Distance along each coordinate direction
		    if (dx <= -_box_middle[coord]){
			dx += _box_size[coord];						// PBC
		    }else if(dx > _box_middle[coord]){
			dx -= _box_size[coord];						// PBC
		    }
		    r2 += (dx*dx);							// Square of the distance component
		}
		if(r2 < _rc_sh2){							// Particles inside the shell radious
		    atoms[i].neighbor.push_back(j);					// Add the particle to the neighbour list
		    //atoms[j].neighbor.push_back(i);					// Unecessary because of second Newton's law
		    if(r2 < _rc2){
			r = sqrt(r2);
			//ELECTRONIC DENSITY 
			eam.set_position(r,2);
			rho = eam.get_electron_density();
			//printf("rho_par = %Lf - ",rho);
			atoms[i].n_i += rho;
			atoms[j].n_i += rho;
			//PAIR POTENTIAL 
			eam.set_position(r,0);
			u_val = eam.get_pair_potential();
			u_energy += u_val;
		    }
		}
	    }
	    j = _cell_list[j];								// Change to the next atom in the list
	}
    }
#ifdef _MD_DEBUG_NEIGHBOR_
  show_message("End mdTargetTargetNeighbors");
#endif
}

void Moldyn::eval_neighbor_list(void){
    lreal u_val;
    eval_cell_list();
#ifdef _MD_DEBUG_NEIGHBOR_
  show_message("Beging get_neighbor");
#endif
    
    //u_energy=0;
    for(unsigned int i=0; i<PARTICLE_NUMBER; i++){
	eval_neighbor(i);
    }

    //printf("Acumulated Energy = %Lf\n",u_energy);
    for(unsigned int i=0; i<PARTICLE_NUMBER; i++){
	if(atoms[i].n_i<eam.eamp[1][1]){
	    eam.set_position(atoms[i].n_i,1);
	    atoms[i].n_i = eam.get_embedded_derivative();
	    u_val = eam.get_embedded_energy();
	    u_energy += u_val;
	    //printf(" u_val_emb = %Lf\n ",u_val);
	}else{atoms[i].n_i = 0;}
    }
#ifdef _MD_DEBUG_NEIGHBOR_
  show_message("End get_neighbor");
#endif
}

// Gives number of atoms closer than sqrt(r2max) from atom number i
void Moldyn::eval_embedded_energy(int i){
    lreal dx, dr[3], r2, r, u_val, rho;
    unsigned int coord;
#ifdef _MD_DEBUG_NEIGHBOR_
  show_message("End mdTargetTargetNeighbors");
  cout<<" ("<<i<<")"<<endl;
#endif
    for(unsigned int j=0; j<atoms[i].neighbor.size(); j++){
	r2 = 0;
	for(coord=0; coord<_DIM; coord++){
	    dx = atoms[atoms[i].neighbor[j]].r[coord] - atoms[i].r[coord];		// Distance along a coordinate direction
	    if (dx <= -_box_middle[coord]){
	        dx += _box_size[coord];							// PBC
	    }else if(dx >  _box_middle[coord]){
	        dx -= _box_size[coord];							// PBC
	    }
	    dr[coord]=dx;
	    r2 += (dx*dx);
	}
	if(r2 < _rc2){
	    r = sqrt(r2);
	    //ELECTRONIC DENSITY 
	    eam.set_position(r,2);
	    rho = eam.get_electron_density();
	    atoms[i].n_i += rho;
	    atoms[atoms[i].neighbor[j]].n_i += rho;
	    //PAIR POTENTIAL 
	    eam.set_position(r,0);
	    u_val = eam.get_pair_potential();
	    u_energy += u_val;    
	}
    }
#ifdef _MD_DEBUG_NEIGHBOR_
  show_message("End mdTargetTargetNeighbors");
#endif
}

void Moldyn::eval_embedded_energy(void){
    lreal u_val;
    //eval_cell_list();
#ifdef _MD_DEBUG_NEIGHBOR_
  show_message("Beging get_neighbor");
#endif
    //u_energy=0;
    for(unsigned int i=0; i<PARTICLE_NUMBER; i++){
	eval_embedded_energy(i);
    }
    //printf("Acumulated Energy = %Lf\n",u_energy);
    for(unsigned int i=0; i<PARTICLE_NUMBER; i++){
	//printf(" - %i ",(int)atoms[i].neighbor.size());
	if(atoms[i].n_i<eam.eamp[1][1]){
	    eam.set_position(atoms[i].n_i,1);
	    atoms[i].n_i = eam.get_embedded_derivative();
	    u_val = eam.get_embedded_energy();
	    u_energy += u_val;
	}else{atoms[i].n_i = 0;}
    }
#ifdef _MD_DEBUG_NEIGHBOR_
  show_message("End get_neighbor");
#endif
}

void Moldyn::eval_eam_force(int i){
//#ifdef __DEBUG_FUNCTIONS__
  //show_message("eval_eam_force");
//#endif
    lreal dx, dn_i, dpot_par;
    lreal dpot_emb_j, dpot_emb_i;
    lreal r2, r;
    lreal dr[3], force;
    unsigned int coord;
    sum_vir = 0;
    // This function use neighbor list in order to speed up this force computation
    for(unsigned int j=0; j<atoms[i].neighbor.size(); j++){
	r2 = 0;
	for (coord=0; coord<_DIM; coord++){
	    dx = atoms[atoms[i].neighbor[j]].r[coord] - atoms[i].r[coord];		// Distance along a coordinate direction
	    if (dx <= -_box_middle[coord]){
	        dx += _box_size[coord];							// PBC
	    }else if(dx >  _box_middle[coord]){
	        dx -= _box_size[coord];							// PBC
	    }
	    dr[coord]=dx;
	    r2 += (dx*dx);
	}
	if(r2 < _rc2){
	    r = sqrt(r2);
	    //DENSITY DERIVATIVE
	    eam.set_position(r,2);
	    dn_i  = eam.get_electron_derivative();
	    //PAIR POTENTIAL DERIVATIVE
	    eam.set_position(r,0);
	    dpot_par = eam.get_pair_derivative();
	    //
	    dpot_emb_i = atoms[i].n_i;
	    dpot_emb_j = atoms[atoms[i].neighbor[j]].n_i;
	    force = (dpot_par + (dpot_emb_i + dpot_emb_j) * dn_i)/r;
	    //printf("Force = %Lf - ",force);
	    sum_vir += r*force;
	    for(unsigned int k=0; k<_DIM; k++){
		dpot_par = dr[k]*force;
		atoms[i].a[k] += dpot_par;						// Force of particle j on i
		atoms[atoms[i].neighbor[j]].a[k] -= dpot_par;				// Second Newtons's law (isn't it?)
	    }
	    
	}
    }
//#ifdef __DEBUG_FUNCTIONS__
  //show_message("eval_eam_force");
//#endif
}

void Moldyn::eval_force(void){
#ifdef _MD_DEBUG_FORCES_
  show_message("Begin eval_force");
#endif
    for(unsigned int i=0; i<PARTICLE_NUMBER; i++){
	eval_eam_force(i);
    }
#ifdef _MD_DEBUG_FORCES_
  show_message("End eval_force");
#endif
}

////////////////////////// SYSTEM MOTION INTEGRATION ///////////////////////////

void Moldyn::integration_verlet(void){
    unsigned int i, j;
    lreal vv, maxvv;
    _vv_sum = 0;
#ifdef _MD_DEBUG_VERLET_
  show_message("Begin integration_verlet");
#endif
    // Start engines with Verlet algorithm to integrate the motion equations
    //set_next_iteration();
    //get_neighbor_list();
    //eval_force();
    // first step of motion
    for(i=0; i<PARTICLE_NUMBER; i++){
	//if(atoms[i]._sw){
        for(j=0; j<_DIM; j++){
	    atoms[i].v[j] += 0.5*atoms[i].a[j]*_time_step;
	    atoms[i].r[j] += atoms[i].v[j]*_time_step;
	    
	    if(atoms[i].r[j] <= -_box_middle[j])						// PBC
		atoms[i].r[j] += _box_size[j];
	    else if(atoms[i].r[j] >  _box_middle[j])						// PBC
		atoms[i].r[j] -= _box_size[j];
	}
	//}
    }
    set_next_iteration();
#ifdef _USE_SHELL_LIST_
    if(_sum_disp>_sh_2) _if_neighbor=true;
    if(_if_neighbor){
	eval_neighbor_list();
	_if_neighbor = false;
	//printf("\nd=%Lf > %Lf ",_sum_disp,_sh_2);
	_sum_disp = 0;
	//printf("atoms restarted\n");
    }else{
	eval_embedded_energy();
	//printf("No neighbor calculated\n");
    }
#else
    eval_neighbor_list();
#endif
    eval_force();
    //
    maxvv = 0;
    // second step of motion
    for(i=0; i<PARTICLE_NUMBER; i++){
	//if(atoms[i]._sw){
	vv = 0;
    	for(j=0; j<_DIM; j++){
	    atoms[i].v[j] += 0.5*atoms[i].a[j]*_time_step;
#ifdef _USE_SHELL_LIST_
	    vv += atoms[i].v[j]*atoms[i].v[j];
#endif
	}
	_vv_sum += vv;
#ifdef _USE_SHELL_LIST_
	maxvv = max(maxvv, vv);
#endif
	//}
    }
#ifdef _USE_SHELL_LIST_
    _sum_disp += sqrt(maxvv)*_time_step;
#endif
/*
#ifdef _MD_MSG_
	printf("\n*** Time step ***");
#endif
    _dt_ = _time_step;
#ifdef _USE_TIME_STEP_
    mdSetTimeStep(_max_dis/_mxv);
#ifdef _MD_DPL_INF_
    if(_show_time_step)
		cout<<"Time Step="<<mdGetTimeStep()<<endl;
#endif    
#endif
	//adjust_temperature(); // Corrije la energia cinetica
    //velocity_distribution();
    //eval_properties();
*/
/*
    _dt_ = _time_step;
#ifdef _MD_TARGET_
    adjust_temperature();
#endif    
    //mdResetParticles();
	set_next_iteration_accel();
#ifdef _MD_TARGET_
    eval_neighbor();
    eval_force();
#endif    
*/
#ifdef _MD__MD_DEBUG_VERLET_
  show_message("End integration_verlet");
#endif
}

// Heart of the program; the motion of all the particles is computed in each step time...
// Verlet Velocity with Neighbor List.
void Moldyn::motion_step(void){
    lreal e_pot, e_kin, temp;
#ifdef _MD_MSG_
  show_message("*** Begin Motion Step ***");
#endif

#ifdef _MD_TIMER_
    timer.start();
#endif
    u_energy=0;
    k_energy=0;
    //
    integration_verlet();
    //
    eval_u_energy();
    eval_k_energy();
    //
    step_counter++;
    if(step_counter == avg_steps){
	e_kin = sum_k_energy/lreal(avg_steps);
	e_pot = sum_u_energy/lreal(avg_steps);
	temp = eval_temperature();
        printf("<K>: %0.3Lf\t<U>: %0.3Lf eV\t<E>: %0.3Lf\t[T]: %0.3Lf K\n",e_kin,(e_pot*K_B/eV),e_kin+e_pot,temp);
        sum_k_energy = 0;
        sum_u_energy = 0;
        step_counter=0;
#ifdef _MD_TIMER_
	printf("Elapsed time: %f\n",sum_time);
	sum_time = 0;
#endif
    }
    
#ifdef _MD_TIMER_
    timer.stop();
    sum_time += timer.get_time();
#endif

#ifdef _MD_MSG_
  show_message("*** End Motion Step ***");
#endif
}

void Moldyn::eval_eos_lattice(void){
    // Initial atom positions
    set_intial_conditions();
    //SHOW SOME INFORMATION//
    //show_summary();
    //motion_step();
    //get_k_energy();
}

////////////////////////// PHYSICS PROPERTIES //////////////////////////////////

void Moldyn::eval_properties(void){
    //lreal _vx=0, _vy=0, _vz=0;
    lreal vv;
    _vv_sum = 0.0;
    for(unsigned int i=0; i<PARTICLE_NUMBER; i++){
	vv=0;
	//_vx += atoms[i].v[0];
	//_vy += atoms[i].v[1];
	//_vz += atoms[i].v[2];
	for(unsigned int j=0; j<_DIM; j++)
	    vv += (atoms[i].v[j]*atoms[i].v[j]);
	_vv_sum += vv;
    }
    //_ue = evalu_energy();
    //_ke = eval_k_energy();
    //_temp = eval_temperature();
#if _MD_DEBUG_ENER
    printf("Energy=%f\n",_ke+_ue);
    printf("Temp=%f\n",_temp);
#endif
	
}

void Moldyn::velocity_distribution(void){
#ifdef _MD_MSG_
  show_message("Begin velocity_distribution");
#endif
    //printf("Start His\n");
    lreal deltaV, _vvmax=0, vv;
    int j;
	
    if(countVel == 0)
	for(j=0; j<sizeHistVel; j++) histVel[j]=0;
    deltaV = rangeVel / sizeHistVel+1;
    //printf("deltaV = %f, rangeVel %f, sizeHistVel %i ,",deltaV,rangeVel,sizeHistVel);
    for(unsigned int i=0; i<PARTICLE_NUMBER; i++){
	vv=0;
	for(unsigned int m=0; m<_DIM; m++)
	    vv += atoms[i].v[m]*atoms[i].v[m];
	//cout<<vv<<" "<<deltaV<<endl;
	j = int(sqrt(vv)/deltaV);
	++histVel[min(j, sizeHistVel-1)];
	_vvmax = max(_vvmax,histVel[min(j, sizeHistVel-1)]);
    }
    //cout<<"max vel="<<rangeVel/deltaV<<endl;
    ++countVel;
    if(countVel==limitVel){
	doVelGrph = 1;
	countVel=0;
    }
    for(j=0; j<sizeHistVel; j++){
	histVel[j] /= _vvmax;
	//cout<<"-"<<histVel[j];
    }
#ifdef _MD_MSG_
  show_message("End velocity_distribution");
#endif
	
}

// EVALUATION FUNCTIONS

void Moldyn::eval_rho(void){
    lreal vol = get_volume();
    _rho = (lreal(PARTICLE_NUMBER)/vol);
}

// Potencial energy of the crystaline structure
void Moldyn::eval_u_energy(void){
    u_energy = u_energy/real(PARTICLE_NUMBER);
    sum_u_energy += u_energy;
}

// Kinetic energy of the crystaline structure
void Moldyn::eval_k_energy(void){
    k_energy = 0.5*_vv_sum/lreal(PARTICLE_NUMBER);
    sum_k_energy += k_energy;
}

// Kinetic energy of the crystaline structure
lreal Moldyn::eval_temperature(void){
    lreal sysTemp=get_k_energy();
    return (sysTemp/3.0);
    //return (sysTemp*1.3806503e-23)/(PARTICLE_NUMBER*1.60217733e-19);
}

// GET FUNCTIONS

lreal Moldyn::get_rho(void){
    return _rho;
}

lreal Moldyn::get_volume(void){
    return vVol(_box_size);
}

lreal Moldyn::get_latt_const(void){
    lreal side = pow(lreal(PARTICLE_NUMBER)/_rho, lreal(1.0/3.0));
    return  side / lreal(_cell_number[0]);
}

// Potencial energy of the crystaline structure
lreal Moldyn::get_u_energy(void){
    return u_energy;
}

// Kinetic energy of the crystaline structure
lreal Moldyn::get_k_energy(void){
    return k_energy;
}

////////////////////////// CONTROL FUNCTIONS ///////////////////////////////////////
// prepare atoms for the next time iteration 
void Moldyn::set_next_iteration(void){
    unsigned int i, j;
    for (i=0; i<PARTICLE_NUMBER; i++){
	atoms[i].n_i=0;
	//atoms[i].pot_par=0;
	//atoms[i].pot_emb=0;
	for(j=0; j<_DIM; j++){
	    //atoms[i].v[i] = 0;
	    atoms[i]._a[j]=atoms[i].a[j];
	    atoms[i].a[j]=0;
	}
    }
}

// READ CONFIGURATION FILES
// parameters of the simulation
void Moldyn::read_md_config_file(void){
#ifdef __DEBUG_FUNCTIONS__
  show_message("read_md_config_file");
#endif
  CRead *config = new CRead();
  config->read_file_names();
  config->read_md_config_file();
  config->summarize_parameters();

#ifdef __DEBUG_MESSAGES__
  show_message("load the EAM potencial");
#endif
  eam.read_eam_file(config->get_eam_name());
#ifdef __DEBUG_MESSAGES__
  show_message("set simulation initial parameters");
#endif
  set_atomic_number(config->get_atomic_number());
  set_atomic_mass(config->get_atomic_mass());
  set_lattice_constant(config->get_lattice_constant());
  set_basis_atom_number(config->get_basis_atom_number());
  set_temperature(config->get_temperature());
  //
  set_time_step(config->get_time_step());
  set_cut_radius(config->get_cut_radius());
  set_shell_width(config->get_shell_thick());
  set_md_steps(config->get_md_steps());
  //
  set_lattice_vectors(config->get_lattice_vectors());
  set_basis_vectors(config->get_basis_vectors());
  // the number of unit cells in each crystal axis direction
  set_cells(config->get_cell_number());
  //
  set_cut_radio(_cut_radius,_shell_thick); // solid
  //
  //set_intial_temperature(_temperature);
  // if the density is given then the lattice constant is calculated
  if(_rho > 0){
    //cout<<" _rho = "<<_rho<<endl;
    //_lattice_constant = 1.0/pow(_rho,1.0/3.0);
    _lattice_constant = get_latt_const();
    cout<<" lattice constant = "<<_lattice_constant<<endl;
  }
  for(unsigned int i=0; i<_DIM; i++){
    // Size of the box
    _box_size[i]   = (_cell_number[i]*_lattice_constant*_lattice_vectors[i][i]);
    // Half size of the box with particles
    _box_middle[i] = (_cell_number[i]*_lattice_constant*0.5*_lattice_vectors[i][i]);
  }
  // if the density is NOT give then is calculated
  if(_rho == 0){
    eval_rho();
    //cout<<" _rho = "<<_rho<<endl;
    //cout<<" lattice constant = "<<_lattice_constant<<endl;
  }
  ///////////////////////
  //_show_inf = 0;
  //_show_time_step = 0;
  //_max_disp = 0.1;
  ///////////////////////
  //set_reduced_units();
#ifdef __DEBUG_FUNCTIONS__
  show_message("read_md_config_file");
#endif
}

/////////////////////// VECTOR UTILS ///////////////////////////////////////

vector<lreal> Moldyn::VRot(vector<lreal>& _v, lreal _a, lreal _b, lreal _g){
    vector<lreal> _t(_v.size());
    lreal _s, _c;
    _s = sin(_a);
    _c = cos(_a);
    //_t[0] = _v[0];
    _t[1] = ( _c*_v[1] + _s*_v[2]);
    _t[2] = (-_s*_v[1] + _c*_v[2]);
    
    _s = sin(_b);
    _c = cos(_b);
    _t[0] = ( _c*_v[0] + _s*_t[2]);
    //_t[1] = _t[1];
    _t[2] = (-_s*_v[0] + _c*_t[2]);
	/*  _s = sin(_g);
    _c = cos(_g);
    _t[0] = ( _c*_t[0] + _s*_t[1]);
    _t[1] = (-_s*_t[0] + _c*_t[1]);
    //_t[2] = _v[2];*/
    return _t;
}

void Moldyn::VRand(unsigned int _s, lreal _scl, lreal _t[3]){
    unsigned int i;
    //dvector _t;
    //_t.resize(_s);
    for(i=0; i<_s; i++){
	_t[i] = _scl*((2.0*rand()/(RAND_MAX+1.0))-1.0);
    }
    //return _t;
}

void Moldyn::VRandNorm(unsigned int _s, lreal _scl, lreal _t[3]){
    unsigned int i;
    lreal _vnor=0;
    //dvector _t;
    //_t.resize(_s);
    VRand(_s,_scl,_t);
    for(i=0; i<_s; i++)
	_vnor+=(_t[i]*_t[i]);
    for(i=0; i<_s; i++)
	_t[i] *= _scl/sqrt(_vnor);
    //return _t;
}

/////////////////////////// INFORMATION FUNCTIONS

void Moldyn::write_info_file(string file){
#ifdef __DEBUG_FUNCTIONS__
  show_message("write_info_file");
#endif
    double tm;
    std::fstream md;
    md.open(file.c_str(), ios::out);
    md.precision(6);
    md.setf(ios::scientific | ios::fixed);
    ////////////////////////////OUTPUT FILE///////////////////////////////////
    md<<"Constante de red\t=\t"<<_lattice_constant<<"\tA"<<endl;
    md<<"Masa atomica del Fe\t=\t"<<_atomic_mass<<"\tk"<<endl;
    md<<"Estrcutura Cristalina\t=\t"<<_struc_type.c_str()<<endl;
    md<<"Paso del tiempo (h) \t=\t"<<_time_step<<endl;
    tm = _time_step * ANGSTROM * sqrt((_atomic_mass*AMU)/K_B);
    md<<"Paso lreal de tiempo \t=\t"<<tm<<"\ts"<<endl;
    md<<"\nNumero de particulas\t=\t"<<atoms.size()<<" particulas"<<endl;
    //////////////////////////////////////////////////////////////////////////
    md.close();
#ifdef __DEBUG_FUNCTIONS__
  show_message("write_info_file");
#endif
}

void Moldyn::show_introduction(void){
    cout<<"** 	Programa de simulacion de Dinamica		**"<<endl;
    cout<<"**	Molecular de alfa-Hierro.			**"<<endl;
    cout<<"**	Edmanuel Torres A.				**"<<endl;
}

void Moldyn::summarize_parameters(void){
    ////////////////////////////////////////
    cout<<"<-----------SUMMARISE----------->"<<endl;
    cout<<"Atomic number = "<<_atomic_number<<endl;
    cout<<"Atomic mass = "<<_atomic_mass<<endl;
    cout<<"Lattice constant= "<<_lattice_constant<<endl;
    cout<<"Lattice_vectors : "<<_lattice_vectors<<endl;
    cout<<"Unit cell atoms = "<<_basis_atom_number<<endl;
    cout<<"Basis atoms : "<<_basis_vectors<<endl;
    cout<<"Shift : "<<_shift<<endl;
    ///////////////////////
    printf("Reduced Energy   = %Le\n",_epsilon);
    printf("Reduced Length = %Le\n",_sigma);
    printf("Reduced Mass      = %Le\n",_mu);
    printf("Reduced Time    = %Le\n",_tao);
    printf("Real Density =  %Le (kg/m^3) = %Le (g/cm^3)\n",density,0.001*density);
    cout<<"<-----------SUMMARISE----------->"<<endl;
}

void Moldyn::show_summary(void){
#ifdef __DEBUG_FUNCTIONS__
  show_message("show_summary");
#endif
    ///////////////////////////SHOW SOME INFORMATION//////////////////////////
    cout<<endl;
    //cout<<"[Atoms]"<<endl;
    //cout<<"Simbolo ["<<_quim_simb.c_str()<<"]\t"<<"Numero Atomico ["<<_atomic_number<<"]"<<endl;
    cout<<endl; 
    cout<<"|Name\t\t\t|\tValue\t\t\t|Unit\t\t|"<<endl;
    cout<<endl;
    cout<<"Laticce Constant\t=\t"<<_lattice_constant<<"\tA"<<endl;
    cout<<"Atomic Mass \t=\t"<<_atomic_mass<<"\tuma"<<endl;
//	cout<<"Estrcutura Cristalina\t=\t"<<_struc_type.c_str()<<endl;
    cout<<"Time step (h) \t=\t"<<_time_step<<endl;
//	tm = _time_step * 1e-10 * sqrt((_atomic_mass*AMU)/K_B);
//	cout<<"Paso lreal de tiempo \t=\t"<<tm<<"\ts"<<endl;
    //
    cout<<"\nNumber of particles\t=\t"<<PARTICLE_NUMBER<<" particles"<<endl;
    cout<<"Volume \t\t\t= \t"<<get_volume()<<endl;
    cout<<"Density \t\t= \t"<<get_rho()<<endl;
    cout<<"Cut radio\t\t=\t"<<_rc<<endl;
    cout<<"Cut radio^2\t\t=\t"<<_rc<<endl;
    cout<<"Shell width\t\t=\t"<<_sh<<endl;
    cout<<"Z side \t\t=\t"<<_box_size[0]<<endl;
    cout<<"Y side \t\t=\t"<<_box_size[1]<<endl;
    cout<<"X side \t\t=\t"<<_box_size[2]<<endl;
    cout<<"Sample size \t\t=\t("<<_cell_number[0]<<"X"<<_cell_number[1]<<"X"<<_cell_number[2]<<") [unit cells]"<<endl;
    cout<<"Cells \t\t=\t("<<_cell_side[0]<<"X"<<_cell_side[1]<<"X"<<_cell_side[2]<<")"<<"="<<CELL_NUMBER<<endl;
    cout<<"Cell list size \t=\t"<<_cell_list.size()<<endl;
//	cout<<"Bytes por particula\t=\t"<<mdGetPartAloc()<<" bytes"<<endl;
//	cout<<"Capacidad Maxima\t=\t"<<mdGetMaxSize()<<" bytes"<<endl;
//	cout<<"Memoria ocupada  \t=\t"<<mdGetMemAloc()<<" bytes"<<endl<<endl;
    cout<<endl;
#ifdef __DEBUG_FUNCTIONS__
  show_message("show_summary");
#endif
}

///////////////////////////////////// FILES STUFFS //////////////////////////////////////

void Moldyn::write_atom_xyz(string file, int p){
	std::fstream es;
	es.open(file.c_str(), ios::out);
	es.precision(p);
	es.setf(ios::scientific | ios::fixed);
	es<<atoms.size()<<endl;
	es<<"alpha-Fe"<<endl;
	for(unsigned int i=0; i<atoms.size(); i++){
	    es<<"Fe\t";
	    for(unsigned int k=0; k<_DIM; k++)
		es<<atoms[i].r[k]<<tab;
	    es<<endl;
	}
	es.unsetf(ios::scientific | ios::fixed);
	es.close();
}

void Moldyn::write_atom_pdb(string file, int p){
    std::fstream es;
    es.open(file.c_str(), ios::out);
    es.precision(p);
    es.setf(ios::scientific | ios::fixed);
    for(unsigned int i=0; i<atoms.size(); i++){
	for(unsigned int k=0; k<_DIM; k++)
	    es<<atoms[i].r[k]<<tab;
	es<<_atomic_number<<tab<<_atomic_mass<<tab<<"0"<<tab<<_quim_simb<<endl;
    }
    es.unsetf(ios::scientific | ios::fixed);
    es.close();
}


/*
void Moldyn::mdErrorMesage(string msg){
	cout<<endl<<"WARNING: "<<msg.c_str()<<endl;
}
*/

//////////////////////////// DEBUG UTILS ///////////////////////////////////////////////

void Moldyn::eval_nearest_neighbor(int i){
    lreal dx, r2;
    int j;
#ifdef _MD_DEBUG_NEIGHBOR_
  show_message("End mdTargetTargetNeighbors");
  cout<<" ("<<i<<")"<<endl;
#endif
    atoms[i].neighbor.clear();
    _positive_r = lVAdd(atoms[i].r,_box_middle);
    _integer_r = iVMul(_positive_r,_cell_frac);
    //  Here, computing the neighbour list
    for (int _m=0; _m<14; _m++){
        _neighbor_cell = iVAdd(_integer_r,neighbor_cells[_m]);					// Neighbour cell
											// Free _DIM here with PCB
        for (unsigned int coord=0; coord<_DIM; coord++){				// For each _DIM
	    //if(_md_pbc[coord])							// Ask if PBC is used
    	    if(_neighbor_cell[coord] >= _cell_side[coord]){					// Check if  PBC is necessary
		_neighbor_cell[coord]= 0;							// Apply PBC to each  cell
	    }else if(_neighbor_cell[coord] < 0){						// Check if  PBC is necessary
		_neighbor_cell[coord]= _cell_side[coord]-1;					// Apply PBC to each cell
	    }
	}
	_icell = iVLinear(_neighbor_cell,_cell_side);						// Look for the cell
	//if(_icell>=0 && _icell < CELL_NUMBER)						// Inside of the atoms
	j = _cell_head[_icell];								// Cell position
	//else j = -1;									// Out side of cell
	while(1){									// Over all the particles in the cell
	    if(j<0) break;								// Get out of the cell
	    if(_m!=0 || j>i){								// Avoid auto-interaction 
		r2 = 0;									// Set distance to cero
		for (unsigned int coord=0; coord<_DIM; coord++){					// Each _DIM
		    dx = atoms[j].r[coord] - atoms[i].r[coord];				// Distance along each coordinate direction
		    if (dx <= -_box_middle[coord]){
			dx += _box_size[coord];						// PBC
		    }else if(dx >  _box_middle[coord]){
			dx -= _box_size[coord];						// PBC
		    }
		    r2 += (dx*dx);							// Square of the distance component
		}
		if(r2 < _rc2){								// Particles inside the shell radious
 		    atoms[i].neighbor.push_back(j);					// Add the particle to the neighbour list
		    atoms[j].neighbor.push_back(i);					// Unecessary because of second Newton's law
		}
	    }
	    j = _cell_list[j];								// Change to the next atom in the list
	}
    }
#ifdef _MD_DEBUG_NEIGHBOR_
  show_message("End mdTargetTargetNeighbors");
#endif
}

void Moldyn::eval_nearest_neighbor(void){
    eval_cell_list();
#ifdef _MD_DEBUG_NEIGHBOR_
  show_message("Beging nearest_neighbor");
#endif
    for(unsigned int i=0; i<PARTICLE_NUMBER; i++){
	eval_nearest_neighbor(i);
    }
    for(unsigned int i=0; i<PARTICLE_NUMBER; i++){
	cout<<"NNN: ("<<i<<") ="<<atoms[i].neighbor.size()<<endl;
    }
#ifdef _MD_DEBUG_NEIGHBOR_
  show_message("End nearest_neighbor");
#endif
}

void Moldyn::show_message(string msg){
  cout<<" [MESSAGE] "<<msg<<endl;
}

// Return the position atom i in the crystal target
void Moldyn::get_particle_position(lreal p[3], int i){ 
    for(unsigned int j=0; j<_DIM; j++)
	p[j]=atoms[i].r[j];
}


// END

