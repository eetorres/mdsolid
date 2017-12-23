
#ifndef _CFORCE_H_
#define _CFORCE_H_

#include<vector>
#include<msmvtl/tmatrix.h>
#include<msmvtl/tvmath.h>
#include<mdvect.h>

#include<string>
#include<string.h>
#include<iostream>


// T. Hammerschmidt, A. Kersch and P. Vogl

const lreal rho_ri[4] = { 5.09113442, 4.38171405, 2.7, 2.5};
const lreal rho_ai[4] = { 0.5476138534, -0.5512656007, 40.0, -25.0};

const lreal pair_ri[8] = {5.09113442, 5.00767320, 4.67382832, 3.96440795, 3.33844880, 2.950800645, 2.9, 2.7};
const lreal pair_ai[8] = {-0.7857149938, 1.110966254, -0.2994497740, -0.1430612718, 1.025367788, 0.4942930559, 6.0, -3.41};

const lreal F_rho[6] = { 0.00, 10.00, 20.00, 30.00, 40.00, 48.00};
const lreal F[6] = {0.00, -2.45, -4.85, -5.45, -6.30, -6.9282032302755};
const lreal Fpp[6] = { -2.324535081465755E-2, -1.349070162931513E-2, 3.371745570260303E-2, -1.337912118109709E-2, 4.799029021785401E-3, -1.578195999192991E-5};

// T. Hammerschmidt, A. Kersch and P. Vogl

class CForce{

private:    
    unsigned int _tabledat;
    int rk;
    lreal wt;
    
    TVector<lreal> eam_par, eam_emb, eam_den;
    TVector<lreal> deam_par, deam_emb, deam_den;
    TVector<lreal> _inv_delta;

public:
    TMatrix<lreal> eamp, eamp2;
    
    void read_eam_file(std::string);
    void read_eam2();
    TVector<lreal> derivative(TVector<lreal>&, lreal);
    
    void set_number_data(unsigned int nd){ _tabledat=nd;};
    
    void set_position(lreal _r, unsigned int _i){
	lreal ri;
	ri = (_r-eamp[0][_i])*_inv_delta[_i];
	rk = (int)floorl(ri);
	if (rk < 0) rk = 0;                         			// unlikely but just to protect
	wt = ri - (lreal)rk;                   				// fractional part, in [0,1]
    };
    
    //EAM methods
    
    lreal get_electron_density(void){
	//ELECTRONIC DENSITY 
	return (wt* eam_den[rk+1] + (1.0-wt)* eam_den[rk]); 		// do linear interpolation.
    }
    
    lreal get_embedded_energy(void){
	//EMBEDDED ENERGY
	return (wt* eam_emb[rk+1] + (1.0-wt)* eam_emb[rk]);  		// do linear interpolation.
    }
    
    lreal get_pair_potential(void){
	//PAIR POTENTIAL
	return (wt * eam_par[rk+1] + (1.0-wt)* eam_par[rk]);  		// do linear interpolation.
    }
    
    lreal get_embedded_derivative(void){
	//EMBEDDED ENERGY DERIVATIVE
	return (wt* deam_emb[rk+1] + (1.0-wt)* deam_emb[rk]);  		// do linear interpolation.
    }
    
    lreal get_electron_derivative(void){
	//ELECTRONIC DENSITY 
	return (wt* deam_den[rk+1] + (1.0-wt)* deam_den[rk]); 		// do linear interpolation.
    }
    
    lreal get_pair_derivative(void){
	//PAIR POTENTIAL
	return (wt * deam_par[rk+1] + (1.0-wt)* deam_par[rk]);  	// do linear interpolation.
    }
    
    // T. Hammerschmidt, A. Kersch and P. Vogl
    
    lreal step_func(lreal _x){
	return _x>0?1:0;
    }
        
    lreal F0(lreal r){
	int i = ceil(r/10);
	int j = i-1;
        
	//  `i´ and `j´ now bracket the input value of x.
	lreal h = F_rho[j] - F_rho[i];
	//  The xa's must be distinct
	//if (h == 0.0) {fprintf(stderr, "Bad `xa[]´ input to splint()\n"); exit(1);}
	//  Cubic spline polynomial is now evaluated
	lreal A = (F_rho[j] - r) / h;
	lreal B = (r - F_rho[i]) / h;
	return( A*F[i] + B*F[j] + ((A*A*A-A)*Fpp[i] + (B*B*B-B)*Fpp[j]) * (h*h) / 6.0 );
    }

    lreal rho(lreal _r){
	lreal _rho=0, _x;
	for(unsigned int i=0; i<4; i++){
	    _x = (rho_ri[i]-_r);
	    _rho += rho_ai[i]*(_x*_x*_x)*step_func(_x);
	}
	return _rho;
    }

    lreal pair_pot(lreal _r){
	lreal _pair=0, _x=0;
	for(unsigned int i=0; i<8; i++){
	    _x = (pair_ri[i]-_r);
	    _pair += pair_ai[i]*(_x*_x*_x)*step_func(_x);
	}
	return _pair;
    }    
    
};

#endif

