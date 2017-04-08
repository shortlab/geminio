/**************************************
//Based on Table 1 in T. Jourdan's paper: Efficient simulation of kinetics of radiation induced defects: A cluster dynamics approach
*************************************/

#include "MaterialParameters.h"
#include<stdio.h>
#include "MooseError.h"

#define INF 100
#define SCALE 1 //change unit from um
#define PI 3.14159265359
#define Vatom 1.182e-11 //iron atom volume um^3
#define Boltz_const 8.6173315e-5 //boltzmann constant eV/K


double energy(int s,std::string species, std::string Etype) {//unit:eV
    double E=0.0;
    if ((species == "V") && (Etype == "migration")){
///* from literature
        switch(s){
            case 1:
            {
                E = 0.83;
                break;
            }
            case 2:
            {
                E = 0.62;
                break;
            }
            case 3:
            {
                E = 0.35;
                break;
            }
            case 4:
            {
                E = 0.48;
                break;
            }
            default:
                E = INF;
        }
    }
    else if ((species == "I") && (Etype == "migration")){
        switch(s){
            case 1:
            {
                E = 0.34;
                break;
            }
            case 2:
            {
                E = 0.42;
                break;
            }
            case 3:
            {
                E = 0.43;
                break;
            }
            default:
                E = INF;
        }
    }
    else if ((species == "V") && (Etype == "binding")){
        switch(s){
            case 1:
            {
                E = INF;
                break;
            }

            case 2:
            {
                E = 0.30;
                break;
            }
            case 3:
            {
                E = 0.37;
                break;
            }
            case 4:
            {
                E = 0.62;
                break;
            }
            case 5:
            {
                E = 0.73;
                break;
            }
            default:
            {
                E = 2.2 + (0.3-2.2)/(pow(2.0,2.0/3)-1) * (pow(s,2.0/3)-pow(s-1,2.0/3));//capillary law
            }
        }
    }
    else if ((species == "I") && (Etype == "binding")) {
        switch(s){
            case 1:
            {
                E = INF;
                break;
            }
            case 2:
            {
                E = 0.83;
                break;
            }
            case 3:
            {
                E = 0.92;
                break;
            }
            case 4:
            {
                E = 1.64;
                break;
            }
            default:
            {
                //E = 3.8 - 5.06*(pow(s,2.0/3)-pow(s-1,2.0/3));//capillary law
                E = 3.64 - 4.78378*(pow(s,2.0/3)-pow(s-1,2.0/3));//capillary law
            }
        }
    }
    else
        mooseError("Energy not defined for " + Etype + " " + species);
    return E;
}

double D_prefactor(int s, std::string species) {
    double D0 = 8.2e5;//um^2/s
    return D0;
}

//size S1 and S2
double absorb(int S1, int S2, std::string C1, std::string C2,double T, int tag1, int tag2) {
    if(tag1==0 && tag2==0) return 0.0;//tag1, tag2 denotes the mobility of C1 and C2; 1: mobile, 0: immobile
    double r_vi = 0.65e-3;//recombination radius in um
    double r1 = pow(S1*Vatom*3/4/PI,1.0/3); //cluster effective radius
    double r2 = pow(S2*Vatom*3/4/PI,1.0/3); //cluster effective radius
    double D_s1 = D_prefactor(S1,C1)*exp(-energy(S1,C1,"migration")/Boltz_const/T);
    double D_s2 = D_prefactor(S2,C2)*exp(-energy(S2,C2,"migration")/Boltz_const/T);
    return 4*PI*(D_s1*tag1+D_s2*tag2)*(r1+r2+r_vi);
}

double diff(int S1, std::string C1,double T)  {
	return D_prefactor(S1,C1)*exp(-energy(S1,C1,"migration")/Boltz_const/T);
}//in um^2/s

double emit(int S1, int S2, double T, std::string C1, std::string C2, int tag1, int tag2) {
    //for now only consider self species emmision, S1 emits S2, S1==1
    double emit_c = 0.0;
    if (S1 > S2 && S2==1)
        emit_c = absorb(S1,S2,C1,C1,T,tag1,tag2)/(Vatom* pow(SCALE,3)) *exp(-energy(S1,C1,"binding")/Boltz_const/T);//unit:/s
    return emit_c;
}
