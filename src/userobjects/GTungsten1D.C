

/****** This implements 1D movement of SIA clusters (n>=1) ******/
//reference: Cluster Dynamics Models of Irradiation Damage Accumulation in Ferritic Iron Part II: Effects of Reaction Dimensionality
//NOTE: normally n=1 would be 3D. However, in tungsten case, rotation barrier is very high (0.4eV). Probability of rotation is very low at 30K.

#include "MooseMesh.h"
#include "GTungsten1D.h"

#define INF 100
#define SCALE 1 //change unit from um
#define PI 3.14159265359
#define Vatom 1.5825e-11 //tungsten atom volume um^3
#define Burgers 2.7366e-4 //burgers vector (um) (sqrt(3)/2*a0)
#define Rvi 0.65e-3 //vacancy - intersitial reaction distance (um)

#define Boltz_const 8.6173315e-5 //boltzmann constant eV/K

template<>
InputParameters validParams<GTungsten1D>()
{
  InputParameters params = validParams<GMaterialConstants>();
  params.addClassDescription( "Calculate specific material properties");
  return params;
}

GTungsten1D::GTungsten1D(const InputParameters & parameters)
: GMaterialConstants(parameters)
{
  atomic_vol = Vatom;
  Ei_formation = 9.96; //interstitial formation energy eV
  Ev_formation = 3.23; //vacancy formation energy eV
  Eib2 = 2.12; // binding energy for interstitial cluster size 2
  Evb2 = -0.1; // binding energy for vacancy cluster size 2
  Ei_binding_factor = (Eib2-Ei_formation)/ (std::pow(2.0,2.0/3)-1);
  Ev_binding_factor = (Evb2-Ev_formation)/(pow(2.0,2.0/3)-1);
}

void GTungsten1D::initialize()
{
}

void GTungsten1D::execute()
{
}

void GTungsten1D::finalize()
{
}

double GTungsten1D::energy(int s,std::string species, std::string Etype) const{//unit:eV
//ref: Microstructural evolution of irradiated tungsten: Ab initio parameterisation of an OKMC model
    double E=0.0;
    if ((species == "V") && (Etype == "migration")){
        E = 1.66;//eV
    }
    else if ((species == "I") && (Etype == "migration")){
        E = 0.013; //.33;//eV
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
                E = 0.04;
                break;
            case 4:
                E = 0.64;
                break;
            case 5:
                E = 0.72;
                break;
            case 6:
                E = 0.89;
                break;
            case 7:
                E = 0.72;
                break;
            default:
                E = Ev_formation + Ev_binding_factor * (std::pow(s*1.0,2.0/3)-std::pow(s-1.0,2.0/3));//capillary law
            /*
            {
                double Ef=1.77,gamma=1.0/1.6022e-7;//ev/um^2
                double r = pow(s*Vatom*3/4/PI,1.0/3); //cluster effective radius
                return Ef-2*gamma*Vatom/r;
            }
            */
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
                E = 3.02;
                break;
            case 4:
                E = 3.6;
                break;
            case 5:
                E = 3.98;
                break;
            case 6:
                E = 4.27;
                break;
            case 7:
                E = 5.39;
                break;
            default://capillary law
                E = Ei_formation + Ei_binding_factor * (std::pow(s*1.0,2.0/3)-std::pow(s-1.0,2.0/3)) ;
        }
    }
    else
        mooseError("Energy not defined for " + Etype + " " + species);
    return E;
}

double GTungsten1D::D_prefactor(int n, std::string species) const{
//ref: Microstructural evolution of irradiated tungsten: Ab initio parameterisation of an OKMC model
    double D0 = 0.0;//um^2/s
    if (species == "V")
       D0 = 6.0096*std::pow(10,8.0-3.0*n);//um^2/s
    else{
       D0 = 1.0016e5*std::pow(1.0*n,-0.5);
    }
    return D0;
}

//vv reaction; flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile
double GTungsten1D::absorbVV(int S1, int S2, int flag, double T) const{
    double result = 0.0;
    switch(flag){
        case 1:
        {
          double w = pow(48.0*PI*PI/Vatom/Vatom*S2,1.0/3); 
          double D_s1 = D_prefactor(S1,"V")*exp(-energy(S1,"V","migration")/Boltz_const/T);
          result = w*Vatom*D_s1; 
          break;
        }
        case 2:
        {
          double w = pow(48.0*PI*PI/Vatom/Vatom*S1,1.0/3); 
          double D_s2 = D_prefactor(S2,"V")*exp(-energy(S2,"V","migration")/Boltz_const/T);
          result = w*Vatom*D_s2; 
          break;
        }
        case 3:
        {
          int S = (S1>S2)? S1: S2;
          double w = pow(48.0*PI*PI/Vatom/Vatom*S,1.0/3); 
          double D_s1 = D_prefactor(S1,"V")*exp(-energy(S1,"V","migration")/Boltz_const/T);
          double D_s2 = D_prefactor(S2,"V")*exp(-energy(S2,"V","migration")/Boltz_const/T);
          result = w*Vatom*(D_s1+D_s2);//add Vatom for unit concern,  P5/19 in ref
          break;
        }
    }
    return result;
}

//vi reaction; flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile; S1:V ; S2:I
double GTungsten1D::absorbVI(int S1, int S2, int flag, double T) const{
    double result = 0.0;
    if (flag == 1){//mobile 3D vacancy 
        double w = pow(4*PI/Vatom/Burgers*S2,1.0/2); 
        double D_s1 = D_prefactor(S1,"V")*exp(-energy(S1,"V","migration")/Boltz_const/T);
        result = _v_bias*w*Vatom*D_s1; //k_i,j
    }
    else if (flag = 2 || flag ==3){//mobile 1D SIA (dominate); need add additional coefficient to be reaction coefficient k_i,j
        double r1 = pow(3*Vatom*S1/4/PI,1.0/3);
        double r2 = pow(Vatom*S2/PI/Burgers,1.0/2);
        result = PI*pow(r1+r2+Rvi,2);
    }
    return result;
}

//ii reaction; flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile
double GTungsten1D::absorbII(int S1, int S2, int flag, double T) const{
    double result = 0.0;
    double r1 = pow(Vatom*S1/PI/Burgers,1.0/2);
    double r2 = pow(Vatom*S2/PI/Burgers,1.0/2);
    if (r1+r2<Rvi)
       result = PI*pow(r1+r2+Rvi,2);
    else
       result = PI*pow(r1+r2+Rvi,2)-PI*pow(r1+r2-Rvi,2);
    return result;
}

double GTungsten1D::diff(int S1, std::string C1,double T) const {
    return D_prefactor(S1,C1)*exp(-energy(S1,C1,"migration")/Boltz_const/T);
}//in um^2/s

double GTungsten1D::emit(int S1, int S2, double T, std::string C1, std::string C2, int tag1, int tag2) const{
    //for now only consider self species emmision, S1 emits S2, S1==1
    if (C1 == "I") return 0.0;//intersitial cluster doesnt' emit.
    double emit_c = 0.0;
    if (S1 > S2 && S2==1)
        emit_c = absorbVV(S1-1,S2,T,tag1+2*tag2)/(Vatom) *exp(-energy(S1,C1,"binding")/Boltz_const/T);
    //unit:/s only emit point defect of the same species 
    return emit_c;
}

double GTungsten1D::disl_ksq(int S1, std::string C1, double T, int tag) const {
    double bias = (! C1.compare("V"))? _v_bias : _i_bias;
    return tag * diff(S1,C1,T) * _rho_d * bias;
}
