//========================================================================
// Mini Scientific Matrix and Vector Template Library - MSMVTL
// MSMVTL is a simple  matrix and vector template library
//
// FILE: linalg.h
//       linear algebra
//
// Copyrigth 2002-2012 by Edmanuel Torres, eetorres@gmail.com
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

#ifndef _LNALG_H_
#define _LNALG_H_

#include <msmvtl/tvector.h>
#include <msmvtl/tmatrix.h>
#include <msmvtl/laext.h>

//////////////////////////////////////////////////////////////////////////////

template<class Tv, class Tr> Tv vrot_y(const Tv& v, const Tr r){
  // cos  0  sin //
  // 0    1  0   //
  //-sin  0  cos //
  Tv _v1;
  TMatrix<Tr> _R(3,3);
  _R[0][0]= cos(r);
  _R[0][2]= sin(r);
  _R[1][1]= 1;
  _R[2][0]=-sin(r);
  _R[2][2]= cos(r);
  _v1 = _R * v;
  return _v1;
}

template<class Tv, class Tr> Tv vrot_z(Tv& v, Tr r){
  // cos -sin 0 //
  // sin  cos 0 //
  // 0    0   1 //
  Tv _v1;
  TMatrix<Tr> _R(3,3);
  _R[0][0]= cos(r);
  _R[0][1]=-sin(r);
  _R[1][0]= sin(r);
  _R[1][1]= cos(r);
  _R[2][2]= 1;
  _v1 = _R * v;
  return _v1;
}

template<class Tv> Tv v_xyz_to_rpt(const Tv& v){
  Tv _v1(3);
  _v1[0]=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); // r
  if(_v1[0]!=0){
    _v1[1]=acos(v[2]/_v1[0]);                 // phi
    _v1[2]=atan2(v[1],v[0]);                  // theta
  }else{
    _v1[1]=_v1[2]=0;
  }
  return _v1;
}

template<class Tv> Tv v_rpt_to_xyz(const Tv& v){
  Tv _v1(3);
  _v1[0] = v[0]*cos(v[2])*sin(v[1]); // x
  _v1[1] = v[0]*sin(v[2])*sin(v[1]); // y
  _v1[2] = v[0]*cos(v[1]);           // z
  return _v1;
}

template<class _T> lreal vleng(_T _v){
  _T _sum = 0;
  for(unsigned int i=0; i<_v.size(); i++)
    _sum += _v[i]*_v[i];
  return sqrt(_sum);
}

template<class _T> real vvol(_T _v){
  real prd = 1;
  for( unsigned int i=0; i<_v.size(); i++)
    prd *= _v[i];
  return prd;
}

///////////////////////////////////////////////////////////////////////
     //////////// AREA 51 /////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

#endif

///
