//========================================================================
// Mini Scientific Matrix and Vector Template Library - MSMVTL
// MSMVTL is a simple  matrix and vector template library
//
// FILE: const.h
//       Physical constants
//
// Copyrigth 2007-2012 by Edmanuel Torres, eetorres@gmail.com
//
// Don't hesitate to contact me for any question, suggestion,
// modification or bug. Get the latest version at:
// http://msmvtl.sourceforge.net
//
// This library is free  software;  you  can  redistribute  it and/or
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


#ifndef _CONST_H_
#define _CONST_H_ 1

//#include<math.h>

/*
Internal quantity       Multiply with to obtain SI 
-----------------       --------------------------
       E                1.6022e-19 J
       r                1e-10 m
       F                1.6022e-9 N
       m                1.6606e-27*m(u) kg
       t                10.1805e-15*sqrt(m(u)) s
       v                9822.66/sqrt(m(u)) m/s
       a                9.6485e17/m(u) m/s^2
*/
// Matematical Constants
#ifndef PI
const double PI         = 3.141592653589793238;                 // Pi constant
#endif
const double C_PI       = 3.141592653589793238;			// Pi constant
const double C_2PI      = 2.0*C_PI;
const double C_PI2      = C_PI/2.0;
//
// Convertions constants
const double DEG_RAD    = C_PI/180.0;                             // Convert from degrees to radians
const double RAD_DEG    = 180.0/C_PI;                             // Convert from radians to degrees
//#endif

/// Some useful International System units

const double AMU        = 1.6605402e-27;			// Atomic Mass unity (Kg)
const double AMU_g      = 1.6605402e-24;			// Atomic Mass unity (g)
const double ERG        = 1e-7;					// Ergio unity
const double eV         = 1.60217733e-19;			// Electron Volt
const double PS         = 1e-12;				// (s)
const double A          = 1e-10;				// in (m)
const double A_m        = 1e-10;				// in (m)
const double A_cm       = 1e-8;					// in (cm)
const double ANGSTROM   = 1e-10;				// Angstrom Unity (m)

// Important Physical Constants
const double m_p        = 1.6726231e-27;			// Proton Mass (Kg)
const double m_e        = 9.10938188e-31;			// Proton Mass (Kg)
const double E          = 1.602176462e-19;			// Electron charge (C)
const double K_B        = 1.3806503e-23;			// Bolzman Constant
const double K_B_J_K    = 1.3806503e-23;			// Bolzman Constant (SI)
const double K_B_erg_K  = 1.3806503e-24;			// Bolzman Constant (cgs)
const double A_0        = 0.5291772083e-10;			// Bohr's Radious (A)
const double _A_0       = 0.5291772083;				// Bohr's Radious (dimensionless)
const double E_0        = 8.854187817e-12;			// Epsilon Constant
const double N_A        = 6.0221e23;				// Avogadro's Number
const double h_0        = 1.054571596e-34;			// Reduced Plank's constant (J·s)

// Some Useful Ion Implating Constants 
const double ao         = 0.24005 ;				// 
const double ab         = 0.5291772083e-10;			//
const double Vo         = 2.1879e6;				// Fermi's velocity (2.1879e+008 cm/s)
const double C_i        = 0.5;					// C Imaplant Constant
//const double alpha    = pow(4.0/(9.0*PI),1.0/3.0);		// alpha Imaplant Constant
const double p1         = 1.0/3.0; 				//
const double p2         = 2.0/3.0;				//
//const double backdens = (4.0/2.8664*2.8664*1e-20);		//
const double k          = h_0/m_e;
const double K_0        = (5.0 * eV) / K_B;			// Umbral kinetics ion' energy
const double N_P        = (14.007*AMU)/m_p;			//
/////////////////////////////////////////////////////////



// Units conversion
const double HARTREE_TO_KCALMOL =  627.509;  //1 Hartree 627.509 kcal mol-1 

//

#endif

// END
