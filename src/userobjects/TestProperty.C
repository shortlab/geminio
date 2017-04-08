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
#include "TestProperty.h"

#define INF 100
#define SCALE 1 //change unit from um
#define PI 3.14159265359
#define Vatom 1.182e-11 //iron atom volume um^3
#define Boltz_const 8.6173315e-5 //boltzmann constant eV/K

template<>
InputParameters validParams<TestProperty>()
{
  InputParameters params = validParams<GMaterialConstants>();
  params.addRequiredParam<Real>("dislocation","dislocation density in #/um^2");
  params.addParam<Real>("i_disl_bias",1.1,"dislocation bias for intersitials");
  params.addParam<Real>("v_disl_bias",1.0,"dislocation bias for vacancies");
  params.addClassDescription( "Calculate specific material properties");
  return params;
}

TestProperty::TestProperty(const InputParameters & parameters)
: GMaterialConstants(parameters),
 _rho_d(getParam<Real>("dislocation")),
 _i_bias(getParam<Real>("i_disl_bias")),
 _v_bias(getParam<Real>("v_disl_bias"))
{
printf("TestProperty constructed\n");
}

void TestProperty::initialize()
{
}

void TestProperty::execute()
{
}

void TestProperty::finalize()
{
}

double TestProperty::energy(int s,std::string species, std::string Etype) const{//unit:eV
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

double TestProperty::D_prefactor(int s, std::string species) const{
    double D0 = 8.2e5;//um^2/s
    return D0;
}

//size S1 and S2
double TestProperty::absorb(int S1, int S2, std::string C1, std::string C2,double T, int tag1, int tag2) const{
    if(tag1==0 && tag2==0) return 0.0;//tag1, tag2 denotes the mobility of C1 and C2; 1: mobile, 0: immobile
    double r_vi = 0.65e-3;//recombination radius in um
    double r1 = pow(S1*Vatom*3/4/PI,1.0/3); //cluster effective radius
    double r2 = pow(S2*Vatom*3/4/PI,1.0/3); //cluster effective radius
    double D_s1 = D_prefactor(S1,C1)*exp(-energy(S1,C1,"migration")/Boltz_const/T);
    double D_s2 = D_prefactor(S2,C2)*exp(-energy(S2,C2,"migration")/Boltz_const/T);
    return 4*PI*(D_s1*tag1+D_s2*tag2)*(r1+r2+r_vi);
}

double TestProperty::diff(int S1, std::string C1,double T) const {
	return D_prefactor(S1,C1)*exp(-energy(S1,C1,"migration")/Boltz_const/T);
}//in um^2/s

double TestProperty::emit(int S1, int S2, double T, std::string C1, std::string C2, int tag1, int tag2) const{
    //for now only consider self species emmision, S1 emits S2, S1==1
    double emit_c = 0.0;
    if (S1 > S2 && S2==1)
        emit_c = absorb(S1,S2,C1,C1,T,tag1,tag2)/(Vatom* pow(SCALE,3)) *exp(-energy(S1,C1,"binding")/Boltz_const/T);//unit:/s
    return emit_c;
}

//dislocation sink rate k^2*Cj*Dj, return k^2*Dj in this function, where k^2= z*rho_d, P230/839 Was book
double TestProperty::disl_ksq(int S1, std::string C1, double T, int tag) const {
   double bias = (! C1.compare("V"))? _v_bias : _i_bias;
   return tag * diff(S1,C1,T) * _rho_d * bias;
}

double TestProperty::Ebinding(double large, const char* type, double small) const
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
                //E = 3.8 - 5.06*(pow(s,2.0/3)-pow(s-1,2.0/3));
                E = 3.64 - 4.78378*(pow(s,2.0/3)-pow(s-1,2.0/3));//capillary law
            }
        }
    }
    else 
       mooseError("Wrong argument to retrieve binding energy");

    return E;
}
