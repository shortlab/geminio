//********************pure iron under neutron irradiation*************//
//*****calculate: Fe irradiation,  Neutron-induced swelling and embrittlement of pure iron and pure nickel irradiated in the BN-350 and BOR-60 fast reactors
//*****parameters: Efficient simulation of kinetics of radiation induced defects: A cluster dynamics approach

/****** This implements 1D movement of SIA clusters (n>=1) ******/
//reference: Cluster Dynamics Models of Irradiation Damage Accumulation in Ferritic Iron Part II: Effects of Reaction Dimensionality
//NOTE: normally n=1 would be 3D. This should be included if rotation barrier of crowdion is low.

#include "MooseMesh.h"
#include "GIron1D.h"

#define INF 100
#define SCALE 1 //change unit from um
#define PI 3.14159265359
#define Vatom 1.165e-11 //iron atom volume um^3
#define Burgers 2.4734e-4 //burgers vector (um) (sqrt(3)/2*a0)
#define Boltz_const 8.6173315e-5 //boltzmann constant eV/K
#define Rvi 3.3e-4 //um capture radius term

/*** reference: Influence of the picosecond defect distribution on damage accumulation in irradiated α-Fe ***/

template<>
InputParameters validParams<GIron1D>()
{
  InputParameters params = validParams<GMaterialConstants>();
  params.addClassDescription( "Calculate specific material properties");
  return params;
}

GIron1D::GIron1D(const InputParameters & parameters)
: GMaterialConstants(parameters)
{
//printf("GIron1D constructed\n");
  atomic_vol = Vatom;
  Ei_formation = 3.77; //interstitial formation energy eV
  Ev_formation = 2.07; //vacancy formation energy eV
  Eib2 = 0.8; // binding energy for interstitial cluster size 2
  Evb2 = 0.3; // binding energy for vacancy cluster size 2
  Ei_binding_factor = (Eib2-Ei_formation)/ (std::pow(2.0,2.0/3)-1);
  Ev_binding_factor = (Evb2-Ev_formation)/(std::pow(2.0,2.0/3)-1);
}

void GIron1D::initialize()
{
}

void GIron1D::execute()
{
}

void GIron1D::finalize()
{
}

double GIron1D::energy(int s,std::string species, std::string Etype) const{//unit:eV
    double E=0.0;
    double dE=0.0;//shift value for vacancy migration energies
    double idE= 0.0;//shift value for interstitial migration energies
    if ((species == "V") && (Etype == "migration")){
        switch(s){
            case 1:
                E = 0.67 + dE;
                break;
            case 2:
                E = 0.62 + dE;
                break;
            case 3:
                E = 0.35 + dE;
                break;
            case 4:
                E = 0.48 + dE;
                break;
            default:
                E = INF;
        }
        //if ( s>4 && s<=20)//added for test of vacancy mobility
        //    E = 0.48 + dE;
    }
    else if ((species == "I") && (Etype == "migration")){
        /*
        switch(s){
            case 1:
                //E = 0.34 + idE;
                E = 0.04 + idE;
                break;
            case 2:
                //E = 0.42 + idE;
                E = 0.04 + idE;
                break;
            case 3:
                //E = 0.43 + idE;
                E = 0.04 + idE;
                break;
            case 4:
                //E = 0.43 + idE;
                E = 0.04 + idE;
                break;
            default:
                E = INF;
        }
        */
        E = INF;
        if (s<=10)//added for test of SIA mobility
            E = 0.17;
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

double GIron1D::D_prefactor(int s, std::string species) const{
    double D0 = 8.2e5;//um^2/s
    return D0;
}


//vv reaction; flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile
double GIron1D::absorbVV(int S1, int S2, int flag, double T) const{
    double result = 0.0;
    double r1 = _v_bias*(std::pow(S1*Vatom*3/4/PI,1.0/3)); //cluster effective radius
    double r2 = _v_bias*(std::pow(S2*Vatom*3/4/PI,1.0/3)); //cluster effective radius
    //double r1 = _v_bias*(std::pow(S1*Vatom*3/4/PI,1.0/3)+Rvi); //cluster effective radius
    //double r2 = _v_bias*(std::pow(S2*Vatom*3/4/PI,1.0/3)+Rvi); //cluster effective radius
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
double GIron1D::absorbVI(int S1, int S2, int flag, double T) const{
    double result = 0.0;
    double r1 = _v_bias*(std::pow(S1*Vatom*3/4/PI,1.0/3));
    double r2 = _i_bias*(std::pow(S2*Vatom/Burgers/PI,1.0/2)); //cluster effective radius: loop
    if (flag == 1){
        double D_s1 = D_prefactor(S1,"V")*exp(-energy(S1,"V","migration")/Boltz_const/T);
        result = 4.0*PI*D_s1*(r1+r2);
    }
    else if (flag = 2 || flag ==3){//mobile 1D SIA (dominate); need add additional coefficient to be reaction coefficient k_i,j
        result = PI*pow(r1+r2+Rvi,2);
    }
    return result;
}

//ii reaction; flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile
double GIron1D::absorbII(int S1, int S2, int flag, double T) const{
    double result = 0.0;
    double r1 = _i_bias*(std::pow(S1*Vatom/Burgers/PI,1.0/2)); //cluster effective radius
    double r2 = _i_bias*(std::pow(S2*Vatom/Burgers/PI,1.0/2)); //cluster effective radius
    if (r1+r2<Rvi)
       result = PI*pow(r1+r2+Rvi,2);
    else
       result = PI*pow(r1+r2+Rvi,2)-PI*pow(r1+r2-Rvi,2);
    return result;
}

double GIron1D::diff(int S1, std::string C1,double T) const {
	return D_prefactor(S1,C1)*std::exp(-energy(S1,C1,"migration")/Boltz_const/T);
}//in um^2/s

double GIron1D::emit(int S1, int S2, double T, std::string C1, std::string C2, int tag1, int tag2) const{
    //for now only consider self species emmision, S1 emits S2, S1==1
    if (C1 == "I") return 0.0;//intersitial cluster doesnt' emit.
    double emit_c = 0.0;
    if (S1 > S2 && S2==1)
        emit_c = absorbVV(S1-1,S2,T,tag1+2*tag2)/(Vatom) *exp(-energy(S1,C1,"binding")/Boltz_const/T);//unit:/s only emit point defect of the same species 
    return emit_c;
}

//dislocation sink rate k^2*Cj*Dj, return k^2*Dj in this function, where k^2= z*rho_d, P230/839 Was book
double GIron1D::disl_ksq(int S1, std::string C1, double T, int tag) const {
   double bias = (! C1.compare("V"))? _v_bias : _i_bias;
   return tag * diff(S1,C1,T) * _rho_d * bias;
}


