//********************pure iron under neutron irradiation*************//
//*****calculate: Fe irradiation,  Neutron-induced swelling and embrittlement of pure iron and pure nickel irradiated in the BN-350 and BOR-60 fast reactors
//*****parameters: Efficient simulation of kinetics of radiation induced defects: A cluster dynamics approach


#include "MooseMesh.h"
#include "GIron.h"

#define INF 100
#define SCALE 1 //change unit from um
#define PI 3.14159265359
#define Vatom 1.165e-11 //iron atom volume um^3
#define Burgers 2.4734e-4 //burgers vector (um) (sqrt(3)/2*a0)
#define Boltz_const 8.6173315e-5 //boltzmann constant eV/K
#define R0 3.3e-4 //um capture radius term

/*** reference: Influence of the picosecond defect distribution on damage accumulation in irradiated α-Fe ***/
template<>
InputParameters validParams<GIron>()
{
  InputParameters params = validParams<GMaterialConstants>();
  params.addClassDescription( "Calculate specific material properties");
  return params;
}

GIron::GIron(const InputParameters & parameters)
: GMaterialConstants(parameters)
{
//printf("GIron constructed\n");
  atomic_vol = Vatom;
  Ei_formation = 3.77; //interstitial formation energy eV
  Ev_formation = 2.07; //vacancy formation energy eV
  Eib2 = 0.8; // binding energy for interstitial cluster size 2
  Evb2 = 0.3; // binding energy for vacancy cluster size 2
  Ei_binding_factor = (Eib2-Ei_formation)/ (std::pow(2.0,2.0/3)-1);
  Ev_binding_factor = (Evb2-Ev_formation)/(std::pow(2.0,2.0/3)-1);
}

void GIron::initialize()
{
}

void GIron::execute()
{
}

void GIron::finalize()
{
}

double GIron::energy(int s,std::string species, std::string Etype) const{//unit:eV
    double E=0.0;
    if ((species == "V") && (Etype == "migration")){
        switch(s){
            case 1:
                E = 0.67;
                break;
            case 2:
                E = 0.62;
                break;
            case 3:
                E = 0.35;
                break;
            case 4:
                E = 0.48;
                break;
            default:
                E = INF;
        }
    }
    else if ((species == "I") && (Etype == "migration")){
        switch(s){
            case 1:
                E = 0.34;
                break;
            case 2:
                E = 0.42;
                break;
            case 3:
                E = 0.43;
                break;
            case 4:
                E = 0.43;
                break;
            default:
                E = INF;
        }
    }
    else if ((species == "V") && (Etype == "binding")){
        switch(s){
            case 1:
                E = INF;
                break;
            case 2:
                E = Evb2;
                break;
            case 3:
                E = 0.37;
                break;
            case 4:
                E = 0.62;
                break;
            default:
                E = Ev_formation + Ev_binding_factor * (std::pow(s*1.0,2.0/3)-std::pow(s-1.0,2.0/3));//capillary law
        }

    }
    else if ((species == "I") && (Etype == "binding")) {
        switch(s){
            case 1:
                E = INF;
                break;
            case 2:
                E = Eib2;
                break;
            case 3:
                E = 0.92;
                break;
            case 4:
                E = 1.64;
                break;
            default:
                E = Ei_formation + Ei_binding_factor * (std::pow(s*1.0,2.0/3)-std::pow(s-1.0,2.0/3)) ;
        }
    }
    else
        mooseError("Energy not defined for " + Etype + " " + species);
    return E;
}

double GIron::D_prefactor(int s, std::string species) const{
    double D0 = 8.2e5;//um^2/s
    return D0;
}

//size S1 and S2
double GIron::absorb(int S1, int S2, std::string C1, std::string C2,double T, int tag1, int tag2) const{
    if(tag1==0 && tag2==0) return 0.0;//tag1, tag2 denotes the mobility of C1 and C2; 1: mobile, 0: immobile
    double r0 = 3.3e-4;//recombination radius in um
    double bias1 = (! C1.compare("V"))? _v_bias : _i_bias;
    //double r1 = bias1 * (std::pow(S1*Vatom*3/4/PI,1.0/3)+r0); //cluster effective radius
    double r1 = (! C1.compare("V"))? (bias1*(std::pow(S1*Vatom*3/4/PI,1.0/3)+r0)):(bias1*(std::pow(S1*Vatom/Burgers/PI,1.0/2)+r0)); //cluster effective radius
    double bias2 = (! C2.compare("V"))? _v_bias : _i_bias;
    //double r2 = bias2 * (std::pow(S2*Vatom*3/4/PI,1.0/3)+r0); //cluster effective radius
    double r2 = (! C2.compare("V"))? (bias2*(std::pow(S2*Vatom*3/4/PI,1.0/3)+r0)):(bias2*(std::pow(S2*Vatom/Burgers/PI,1.0/2)+r0)); //cluster effective radius
    double D_s1 = D_prefactor(S1,C1)*exp(-energy(S1,C1,"migration")/Boltz_const/T);
    double D_s2 = D_prefactor(S2,C2)*exp(-energy(S2,C2,"migration")/Boltz_const/T);
    return 4*PI*(D_s1*tag1+D_s2*tag2)*(r1+r2);

}

//vv reaction; flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile
double GIron::absorbVV(int S1, int S2, int flag, double T) const{
    double result = 0.0;
    double r1 = _v_bias*(std::pow(S1*Vatom*3/4/PI,1.0/3)+R0); //cluster effective radius
    double r2 = _v_bias*(std::pow(S2*Vatom*3/4/PI,1.0/3)+R0); //cluster effective radius
    switch(flag){
        case 1:
        {
          double D_s1 = D_prefactor(S1,"V")*exp(-energy(S1,"V","migration")/Boltz_const/T);
          result = 4.0*PI*D_s1*(r1+r2);
          break;
        }
        case 2:
        {
          double D_s2 = D_prefactor(S2,"V")*exp(-energy(S2,"V","migration")/Boltz_const/T);
          result = 4.0*PI*D_s2*(r1+r2);
          break;
        }
        case 3:
        {
          double D_s1 = D_prefactor(S1,"V")*exp(-energy(S1,"V","migration")/Boltz_const/T);
          double D_s2 = D_prefactor(S2,"V")*exp(-energy(S2,"V","migration")/Boltz_const/T);
          result = 4*PI*(D_s1+D_s2)*(r1+r2);
          break;
        }
    }
    return result;
}

//vi reaction; flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile
double GIron::absorbVI(int S1, int S2, int flag, double T) const{
    double result = 0.0;
    double r1 = _v_bias*(std::pow(S1*Vatom*3/4/PI,1.0/3)+R0);
    double r2 = _i_bias*(std::pow(S2*Vatom/Burgers/PI,1.0/2)+R0); //cluster effective radius
    switch(flag){
        case 1:
        {
          double D_s1 = D_prefactor(S1,"V")*exp(-energy(S1,"V","migration")/Boltz_const/T);
          result = 4.0*PI*D_s1*(r1+r2);
          break;
        }
        case 2:
        {
          double D_s2 = D_prefactor(S2,"I")*exp(-energy(S2,"I","migration")/Boltz_const/T);
          result = 4.0*PI*D_s2*(r1+r2);
          break;
        }
        case 3:
        {
          double D_s1 = D_prefactor(S1,"V")*exp(-energy(S1,"V","migration")/Boltz_const/T);
          double D_s2 = D_prefactor(S2,"I")*exp(-energy(S2,"I","migration")/Boltz_const/T);
          result = 4*PI*(D_s1+D_s2)*(r1+r2);
          break;
        }
    }
    return result;
}

//ii reaction; flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile
double GIron::absorbII(int S1, int S2, int flag, double T) const{
    double result = 0.0;
    double r1 = _i_bias*(std::pow(S1*Vatom/Burgers/PI,1.0/2)+R0); //cluster effective radius
    double r2 = _i_bias*(std::pow(S2*Vatom/Burgers/PI,1.0/2)+R0); //cluster effective radius
    switch(flag){
        case 1:
        {
          double D_s1 = D_prefactor(S1,"I")*exp(-energy(S1,"I","migration")/Boltz_const/T);
          result = 4.0*PI*D_s1*(r1+r2);
          break;
        }
        case 2:
        {
          double D_s2 = D_prefactor(S2,"I")*exp(-energy(S2,"I","migration")/Boltz_const/T);
          result = 4.0*PI*D_s2*(r1+r2);
          break;
        }
        case 3:
        {
          double D_s1 = D_prefactor(S1,"I")*exp(-energy(S1,"I","migration")/Boltz_const/T);
          double D_s2 = D_prefactor(S2,"I")*exp(-energy(S2,"I","migration")/Boltz_const/T);
          result = 4*PI*(D_s1+D_s2)*(r1+r2);
          break;
        }
    }
    return result;
}

double GIron::diff(int S1, std::string C1,double T) const {
	return D_prefactor(S1,C1)*std::exp(-energy(S1,C1,"migration")/Boltz_const/T);
}//in um^2/s

double GIron::emit(int S1, int S2, double T, std::string C1, std::string C2, int tag1, int tag2) const{
    //for now only consider self species emmision, S1 emits S2, S1==1
    if (C1 == "I") return 0.0;//intersitial cluster doesnt' emit.
    double emit_c = 0.0;
    if (S1 > S2 && S2==1)
        emit_c = absorb(S1,S2,C1,C1,T,tag1,tag2)/(Vatom) *exp(-energy(S1,C1,"binding")/Boltz_const/T);//unit:/s only emit point defect of the same species 
    return emit_c;
}

//dislocation sink rate k^2*Cj*Dj, return k^2*Dj in this function, where k^2= z*rho_d, P230/839 Was book
double GIron::disl_ksq(int S1, std::string C1, double T, int tag) const {
   double bias = (! C1.compare("V"))? _v_bias : _i_bias;
   return tag * diff(S1,C1,T) * _rho_d * bias;
}
