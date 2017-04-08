/*************************************************/
/*           DO NOT MODIFY THIS HEADER           */
/*                                               */
/*                     BISON                     */
/*                                               */
/*    (c) 2015 Battelle Energy Alliance, LLC     */
/*            ALL RIGHTS RESERVED                */
/*                                               */
/*   Prepared by Battelle Energy Alliance, LLC   */
/*     Under Contract No. DE-AC07-05ID14517      */
/*     With the U. S. Department of Energy       */
/*                                               */
/*     See COPYRIGHT for full restrictions       */
/*************************************************/
#include "MooseMesh.h"
#include "BCCIronProperty.h"

#define INF 100
#define SCALE 1 //change unit from um
#define PI 3.14159265359
#define Vatom 1.181e-11 //iron atom volume um^3
#define Boltz_const 8.6173315e-5 //boltzmann constant eV/K

template<>
InputParameters validParams<BCCIronProperty>()
{
  InputParameters params = validParams<GMaterialConstants>();
  params.addRequiredParam<Real>("dislocation","dislocation density in #/um^2");
  params.addParam<Real>("i_disl_bias",1.1,"dislocation bias for intersitials");
  params.addParam<Real>("v_disl_bias",1.0,"dislocation bias for vacancies");
  params.addClassDescription( "Calculate specific material properties");
  return params;
}

BCCIronProperty::BCCIronProperty(const InputParameters & parameters)
: GMaterialConstants(parameters),
 _rho_d(getParam<Real>("dislocation")),
 _i_bias(getParam<Real>("i_disl_bias")),
 _v_bias(getParam<Real>("v_disl_bias"))
{
printf("BCCIronProperty constructed\n");
}

void BCCIronProperty::initialize()
{
}

void BCCIronProperty::execute()
{
}

void BCCIronProperty::finalize()
{
}

double BCCIronProperty::energy(int s,std::string species, std::string Etype) const{//unit:eV
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
//*/
// for test
//if(s==1) E = 0.9; else E = INF;
    }
    else if ((species == "I") && (Etype == "migration")){
///* data from literature
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
//*/
//for test
//if(s<=20) E = 0.1+(s-1)/19.0*0.7; else E = 1.1;//0.1+(s-1)/19.0*0.7; 

    }
    else if ((species == "V") && (Etype == "binding")){
        if (s==1) {
           printf("called: error\n");
            //std::cout << "Error: size should be larger than 1 to have a binding energy." << std::endl;
            return 0.0;
        }
        switch(s){
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
            default:
            {
                E = 2.2 - 3.2346 * (pow(s,2.0/3)-pow(s-1,2.0/3));
            }
        }
    }
    else if ((species == "I") && (Etype == "binding")) {
        if (s==1) {
            //std::cout << "Error: size should be larger than 1 to have a binding energy." << std::endl;
            printf("called: error\n");
            return 0.0;
        }
        switch(s){
            case 2:
            {
                E = 0.80;
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
                E = 3.8 - 5.06*(pow(s,2.0/3)-pow(s-1,2.0/3));
            }
        }
    }
    else E = INF;
    return E;
}

double BCCIronProperty::D_prefactor(int s, std::string species) const{
    double D0 = 0.0;
    if (species == "V"){
///* from literature
        switch(s){
            case 1:
            {
                D0 = 7.9e-7;//m^2/s
                break;
            }
            case 2:
            {
                D0 = 3.5e-8;
                break;
            }
            default:
                D0 = 0.0;
        }
//*/
//for test
//if(s==1) D0 = 2.0e-7; else D0 = 0.0;
    }
    else if (species == "I") {
///* from literature
        switch (s){
            case 1:
            {
                D0 = 1.3e-8;
                break;
            }
            case 2:
            {
                D0 = 351.6e-8;//this number is a little weird
                break;
            }
            case 3:
            {
                D0 = 12.1e-8;
                break;
            }
            case 4:
            {
                D0 = 12.3e-8;
                break;
            }
            default:
                D0 = 9.0e-7*pow(s,-0.6);

        }
//*/
// for test
//if(s<=20) D0 = 2.0e-7/s; else D0 = 2.0e-7*pow(s,-0.7);

    }
    else D0 = 0.0;
    return D0*1.0e12*SCALE*SCALE;//change m^2/s to um^2/s
}
//size S1 and S2
double BCCIronProperty::absorb(int S1, int S2, std::string C1, std::string C2,double T, int tag1, int tag2) const{
    if(tag1==0 && tag2==0) return 0.0;//tag1, tag2 denotes the mobility of C1 and C2; 1: mobile, 0: immobile
    double r_vi = 0.65e-3*SCALE;//um
    double r1 = pow(S1*Vatom*3/4/PI,1.0/3); //cluster effective radius
    double r2 = pow(S2*Vatom*3/4/PI,1.0/3); //cluster effective radius
    double D_s1 = D_prefactor(S1,C1)*exp(-energy(S1,C1,"migration")/Boltz_const/T);
    double D_s2 = D_prefactor(S2,C2)*exp(-energy(S2,C2,"migration")/Boltz_const/T);
    return 4*PI*(D_s1*tag1+D_s2*tag2)*(r1+r2+r_vi);
}

double BCCIronProperty::diff(int S1, std::string C1,double T) const {
	return D_prefactor(S1,C1)*exp(-energy(S1,C1,"migration")/Boltz_const/T);
}//in um^2/s

double BCCIronProperty::emit(int S1, int S2, double T, std::string C1, std::string C2, int tag1, int tag2) const{
    //for now only consider self species emmision, S1 emits S2, S1==1
    double emit_c = 0.0;
    if (S1 > S2 && S2==1)
        emit_c = absorb(S1,S2,C1,C1,T,tag1,tag2)/(Vatom* pow(SCALE,3)) *exp(-energy(S1,C1,"binding")/Boltz_const/T);//unit:/s
    return emit_c;
}

double BCCIronProperty::disl_ksq(int S1, std::string C1, double T, int tag) const {
   double bias = (! C1.compare("V"))? _v_bias : _i_bias;
   return tag * diff(S1,C1,T) * _rho_d * bias;
}

double BCCIronProperty::Ebinding(double large, const char* type, double small) const
{//binding energy of small cluster, size 1, i.e. point defecs
    int s = (int) large;
    double E = 0.0;
    if(!strcmp(type,"V")){//vacancy type
        switch(s){
            case 1:
            {
                E = 99999;//a huge value => impossible for emission
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
            default:
            {
                E = 2.2 - 3.2346 * (pow(s,2.0/3)-pow(s-1,2.0/3));
            }
        }
    }
    else if (!strcmp(type,"I")) {
        switch(s){
            case 1:
            {
                E = 99999;
                break;
            }
            case 2:
            {
                E = 0.80;
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
                E = 3.8 - 5.06*(pow(s,2.0/3)-pow(s-1,2.0/3));
            }
        }
    }
    else 
       mooseError("Wrong argument to retrieve binding energy");

    return E;
}
