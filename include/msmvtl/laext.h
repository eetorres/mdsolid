//========================================================================
// Mini Scientific Matrix and Vector Template Library - MSMVTL
// MSMVTL is a simple  matrix and vector template library
//
// FILE: laext.h
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

#ifndef _LNALG_EXT_H_
#define _LNALG_EXT_H_

///////////////////////////////////////////////////////////////////////
     //////////// AREA 51 /////////////////////////////////////
///////////////////////////////////////////////////////////////////////
template<class _Ta> bool bVPos(_Ta _iv){
  return bool( _iv[0]>=0 && _iv[1]>=0 && _iv[2]>=0);
}

template<class _T> TVector<lreal> lVCopy(_T _v){
  TVector<lreal> _t(_v.size());
  for( unsigned int i=0; i<_v.size(); i++)
    _t[i] = lreal(_v[i]);
  return _t;
}

template<class _T> TVector<real> fVCopy(_T _v){
  TVector<real> _t(_v.size());
  for( unsigned int i=0; i<_v.size(); i++)
    _t[i] = lreal(_v[i]);
  return _t;
}

template<class _T> real vVol(_T _v){
    real prd = 1;
    for( unsigned int i=0; i<_v.size(); i++)
        prd *= _v[i];
    return prd;
}

template<class _Ta, class _Tb> bool bVComp(_Ta _iv, _Tb  _m){
  return bool( _iv[0]<_m[0] && _iv[1]<_m[1] && _iv[2]<_m[2]);
}

template<class _Ta, class _Tb> int iVLinear(_Ta _iv, _Tb  _m){
  return int(((_iv[2]*_m[1]+_iv[1])*_m[0])+_iv[0]);
}

template<class _T, class _Tu> TVector<lreal> lVMov(_T _v, _Tu  _add){
  TVector<lreal> _t(_v.size());
  for( unsigned int i=0; i<_v.size(); i++)
    _t[i] = lreal(_v[i] + _add);
  return _t;
}

template<class _T, class _Tu> TVector<real> fVMov(_T _v, _Tu  _add){
  TVector<real> _t(_v.size());
  for( unsigned int i=0; i<_v.size(); i++)
    _t[i] = real(_v[i] + _add);
  return _t;
}

template<class _T1, class _T2> TVector<int> iVAdd(_T1 _v1, _T2  _v2){
  TVector<int> _t(_v1.size());
  for( unsigned int i=0; i<_t.size(); i++)
    _t[i] = int(_v1[i] + _v2[i]);
  return _t;
}

template<class _T1, class _T2> TVector<real> fVAdd(_T1 _v1, _T2  _v2){
  TVector<real> _t(_v1.size());
  for( unsigned int i=0; i<_t.size(); i++)
    _t[i] = real(_v1[i]+_v2[i]);
  return _t;
}

template<class _T1, class _T2> TVector<lreal> lVAdd(_T1 _v1, _T2  _v2){
  TVector<lreal> _t(_v2.size());
  for( unsigned int i=0; i<_t.size(); i++)
    _t[i] = lreal(_v1[i]+_v2[i]);
  return _t;
}

template<class _T, class _Tu> TVector<int> iVScale(_T _v, _Tu  _scale){
  TVector<int> _t(_v.size());
  for( unsigned int i=0; i<_v.size(); i++)
    _t[i] = int(_v[i] * _scale);
  return _t;
}

template<class _T, class _Tu> TVector<lreal> lVScale(_T _v, _Tu  _scale){
  TVector<lreal> _t(_v.size());
  for( unsigned int i=0; i<_v.size(); i++)
    _t[i] = lreal(_v[i] * _scale);
  return _t;
}

template<class _T, class _Tu> TVector<real> fVScale(_T _v, _Tu  _scale){
  TVector<real> _t(_v.size());
  for( unsigned int i=0; i<_v.size(); i++)
    _t[i] = real(_v[i] * _scale);
  return _t;
}

template<class _T1, class _T2> TVector<lreal> lVVScale(_T1 _v, _T2 _vs){
  TVector<lreal> _t(_v.size());
  for( unsigned int i=0; i<_v.size(); i++)
    _t[i] = lreal(_v[i] * _vs[i]);
  return _t;
}

template<class _T1, class _T2> TVector<real> fVVScale(_T1 _v, _T2 _vs){
  TVector<real> _t(_v.size());
  for( unsigned int i=0; i<_v.size(); i++)
    _t[i] = real(_v[i] * _vs[i]);
  return _t;
}

template<class _Ta, class _Tb> TVector<int> iVMul(_Ta _v1, _Tb _v2){
  TVector<int> _t(_v1.size());
  for( unsigned int i=0; i<_v1.size(); i++)
    _t[i] = int(_v1[i] * _v2[i]);
  return _t;
}

template<class _Ta, class _Tb> TVector<real> fVDiv(_Ta _va, _Tb  _vb){
  TVector<real> _t(_va.size());
  for( unsigned int i=0; i<_t.size(); i++)
    _t[i] = (real)_va[i] /(real)_vb[i];
  return _t;
}

template<class _Ta, class _Tb> TVector<lreal> lVDiv(_Ta _va, _Tb  _vb){
  TVector<lreal> _t(_va.size());
  for( unsigned int i=0; i<_t.size(); i++)
    _t[i] = (lreal)_va[i] /(lreal)_vb[i];
  return _t;
}

//////////////////////////////////////////////////////////////////////////////

#endif

///
