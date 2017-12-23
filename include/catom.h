
#ifndef _CATOMS_H_
#define _CATOMS_H_

//#include<vector>
#include<msmvtl/tmatrix.h>
#include<msmvtl/tvmath.h>
#include<mdvect.h>

#include<string>
#include<string.h>
#include<iostream>

using namespace std;

//
class CAtom{
public:
    vector<int> neighbor;
    //lreal mass;
    lreal r[3];
    lreal v[3];
    lreal a[3];
    lreal _a[3];
    lreal n_i;
    bool _sw;
    CAtom(){};
    ~CAtom(){};
};

typedef vector< CAtom > Particles;

template<class T> class CSystem{

public:
    Particles atoms;
    lreal k_energy, u_energy;
    lreal sum_k_energy, sum_u_energy, sum_vir;

};

#endif

