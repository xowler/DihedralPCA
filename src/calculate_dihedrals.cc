// Include support for more than 99999 atoms
// Exclude waters
// Support for ACE NMe NH2
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <cmath>
#include "ccxtc.h"

using namespace std;

map<string,string> opts;
vector<int> bb;
const char *usage = ""
"   Usage: \n"
"       do struct.gro file.xtc prefix \n"
"\n"
"   Arguments: \n"
"          struct.gro: GRO file, used to identify the dihedrals\n"
"          file.xtc:   XTC file to analyze\n"
"          prefix:     Sets the prefix for the three output files\n"
"                      dihedrals.xvg:   Dihedrals file\n" 
"                      pcaed.xvg:       PCAed data\n" 
"                      components.json: PCA information\n";

bool check_file(string fname, string ext=""){
    ifstream f(fname.c_str());
    if( !f.is_open() ) {
        cout << "ERR: Not existing file <" << fname << ">" << endl;
        return false;
    }
    if(ext.size()==0) return true;

    if(!(fname.substr(fname.size()-3,3)==ext)){
        cout << "ERR: Wrong filetype (" << ext << ") <" << fname << ">" <<  endl;
        return false;
    }
    return true;
}

int parse_opts(int argc, char* argv[]){
    if(argc < 4) return 1;
    // GRO
    opts["gro"] = argv[1];
    opts["xtc"] = argv[2];
    opts["prefix"] = argv[3];

    return !(check_file(opts["gro"],string("gro"))) || !((check_file(opts["xtc"],string("xtc"))));
    return 0;
}

void selectbb(vector< pair<string,string> > atoms){
    map<string,int> res;
    string last("");
    int c=0, nres=0;
    cout << "Residues identified (Please check)" << endl;
    atoms.push_back(make_pair("fake","fake"));
    for(vector<pair<string,string> >::iterator i=atoms.begin(); i!=atoms.end(); ++i){
        if( (last != i->first) && (res.size()>0) ){
            if(bb.size()>0){
                bb.push_back(res["N"]);
                bb.push_back(res["CA"]);
            }
            bb.push_back(res["C"]);
            cout << "\t" << last;
            if(++nres==5){
                nres = 0;
                cout << endl;
            }
            res.clear();
        }
        res[i->second] = c++;
        last = i->first;
    }
    cout << endl;
}

void readgro(string fname){
    ifstream f(fname.c_str());
    string temp, res, name;
    int natoms;

    vector< pair<string,string> > atoms;
    
    cout << "Reading GRO file: " <<  fname << endl;
    getline(f,temp);
    f >> natoms;
    for(int i=0; i<natoms; i++){
        getline(f,temp);
        f >> res >> name;
        atoms.push_back(make_pair(res,name));
    }

    selectbb(atoms);
}

inline float dot(float *x, float *y){
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

inline void vminus(float *x, float *y, float *z){
    z[0] = x[0] - y[0];
    z[1] = x[1] - y[1];
    z[2] = x[2] - y[2];
}

inline void vmultscalar(float a, float *x, float *v){
    v[0] = x[0]*a;
    v[1] = x[1]*a;
    v[2] = x[2]*a;
}

inline void cross(float *x, float *y, float *z){
    z[0] = x[1]*y[2]-x[2]*y[1];
    z[1] = x[2]*y[0]-x[0]*y[2];
    z[2] = x[0]*y[1]-x[1]*y[0];
}

void pvec(float *x){
    cout << x[0] << " " << x[1] << " " << x[2] << endl;
}

float dihedral(float *c1, float *c2, float *c3, float *c4){
    float v1[3],v2[3],v3[3],v1p[3],v3p[3];
    vminus(c1,c2,v1);
    vminus(c3,c2,v2);
    vminus(c4,c3,v3);
    float v2n = 1.0/dot(v2,v2);
    vmultscalar(v2n*dot(v1,v2),v2,v1p);
    vminus(v1,v1p,v1p);
    vmultscalar(v2n*dot(v3,v2),v2,v3p);
    vminus(v3,v3p,v3p);

    float c = dot(v1p,v3p)/sqrt(dot(v1p,v1p)*dot(v3p,v3p));
    if(c>1) c=1;
    if(c<-1) c=-1;

    float angle = acos(c);
    // cout << endl << endl;
    // pvec(v1p);
    // pvec(v3p);
    // pvec(v2);
    cross(v1p,v3p,v1);
    // pvec(v1);
    if(dot(v1,v2)<0) angle = -angle;
    return angle*180/3.1415926;
}

void dihedrals(float (*x)[3], ofstream *o){
    int i0,i1,i2,i3,i4;
    for(vector<int>::iterator i=bb.begin(); (i+4)!=bb.end(); i+=3){
        i0 = *(i);
        i1 = *(i+1);
        i2 = *(i+2);
        i3 = *(i+3);
        i4 = *(i+4);
        *o << dihedral(x[i0],x[i1],x[i2],x[i3]) << " " 
           << dihedral(x[i1],x[i2],x[i3],x[i4]) << " ";
    }
    *o << "\n";
}

void readxtc(string fname, string prefix){
    cout << "Reading XTC file: " << fname ;
    ccxtc::xtc xtc(fname.c_str());
    cout << " - " << xtc.natoms << " atoms." << endl; 

    ofstream o((prefix + string("_dihedrals.xvg")).c_str());
    // o << ""; // ToDo: Add labels
    int t=0;
    while(xtc.next()){
        dihedrals(xtc.x, &o);
        if((++t%10000)==0) cout << t;
        if((t%1000)==0){
            cout << ".";
            cout.flush();
        }
    }  
}

int main(int argc, char* argv[]){
    if( parse_opts(argc, argv)){
        cout << usage;
        return 1;
    }
    readgro(opts["gro"]);
    readxtc(opts["xtc"],opts["prefix"]);

    return 0;
}
