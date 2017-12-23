// Programado por: Edmanuel E. Torres Amaris
// simulaciones de Dinamica Molecular de Solidos (Iron BCC)


#include <cforce.h>

#define __FORCE_DEBUG_MESSAGES__

// load the EAM potencial for the crystaline structure
void CForce::read_eam_file(std::string eamf){
#ifdef __FORCE_DEBUG_MESSAGES__
    std::cout<<"!!!read_eam_file!!!"<<std::endl;
    std::cout<<" file: "<<eamf<<std::endl;
#endif
    char ch;
    int ndat=0;
    unsigned int j; // n;
    lreal min, max, pot, h[3];
    //vector<lreal> dt;
    std::ifstream es;

    es.open(eamf.c_str());
    
    if(!es.is_open()){
	std::string _error_msg = "The file " + eamf + " cannot be opened";
	std::cerr<<_error_msg <<std::endl;
	exit(1);
    }else{
      std::string _error_msg = "The file " + eamf + " was opened";
      std::cout<<_error_msg <<std::endl;
    }
    eamp.resize(2,3);
    eamp2.resize(2,3);
    //dt.resize(_tabledat);
    std::cout<<"reading...";
    while(!es.eof()){
        std::cout<<" ...";
	do{
	    es.get(ch);
	    if(ch=='~'){
	      es>>_tabledat;
	      std::cout<<" data: "<<_tabledat<<std::endl;
            }//else es.ignore(1024,NLN);
	}while(ch!='~' && !es.eof());
	
	switch(ndat){
	    case 0:
		eam_par.resize(_tabledat);
		break;
	    case 1:
		eam_emb.resize(_tabledat);
		break;
	    case 2:
		eam_den.resize(_tabledat);
		break;
	    default:
		break;
	}
	if(ndat==3) break;
	//if(ch=='~'){
	    es>>min;
	    std::cout<<" min: "<<min<<std::endl;
	    es>>max;
	    std::cout<<" max: "<<max<<std::endl;
	    eamp[0][ndat]=min;
	    eamp[1][ndat]=max;
	    eamp2[0][ndat]=min*min;
	    eamp2[1][ndat]=max*max;
	    h[ndat]=(max-min)/(_tabledat-1.0);
	    _inv_delta.push_back(1/h[ndat]);
	    for(j=0; j<_tabledat; j++){
		es>>pot;
		switch(ndat){
		    case 0:
			eam_par[j]=pot;
			break;
		    case 1:
			eam_emb[j]=pot;
			break;
		    case 2:
			eam_den[j]=pot;
			break;
		    default:
			break;
		}
	    }
	    ndat++;
	    
	//}
    }
    es.close();
    //
    //eam_par = eam_par;
    //eam_emb = eam_emb;
    //std::cout<<eamp;
    std::cout<<"Potential limits: "<<eamp;
    //std::cout<<eam_par.size();
    //eam_par.save_file("eam_par.dat");
    //std::cout<<eam_den.size();
    //eam_den.save_file("eam_den.dat");
    //std::cout<<eam_emb.size();
    //eam_emb.save_file("eam_emb.dat");
    //
    deam_par = derivative(eam_par,h[0]);
    //deam_par.write_file("deam_par.plot");
    //std::cout<<deam_par;
    deam_emb = derivative(eam_emb,h[1]);
    //deam_emb.write_file("deam_emb.plot");
    //std::cout<<deam_den;
    deam_den = derivative(eam_den,h[2]);
    //std::cout<<deam_emb;
    //std::cout<<eam_den;
    //std::cout<<eam_emb;
    //std::cout.flush();
    /////////////////////////////////////////////////////
    //read_eam2();
#ifdef __FORCE_DEBUG_MESSAGES__
    std::cout<<"!!!read_eam_file!!!"<<std::endl;
#endif
}

TVector<lreal> CForce::derivative(TVector<lreal>& _v, lreal _h){
    unsigned int i;
    TVector<lreal> _dv(_v.size());
    for(i=0; i<2; i++){
	_dv[i]=(2.0*_v[i+3]-9.0*_v[i+2]+18.0*_v[i+1]-11.0*_v[i])/(6.0*_h);
    }
    for(i=2; i<(_v.size()-2); i++){
	_dv[i]=(_v[i-2]-8.0*_v[i-1]+8.0*_v[i+1]-_v[i+2])/(12.0*_h);
    }
    for(i=_v.size()-2; i<_v.size(); i++){
	_dv[i]=(11.0*_v[i]-18.0*_v[i-1]+9.0*_v[i-2]-2.0*_v[i-3])/(6.0*_h);
    }
    return _dv;
}

// END
