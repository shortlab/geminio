/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

//*******************make pair to loop up all the group constants   '+': vacancy; '-': intersitial************************//
// calculate group constant based on hypothetical shape functions

#include "GGroup.h"
#include<math.h>
#include<algorithm>
#define DEBUG 0

template<>
InputParameters validParams<GGroup>()
{
  InputParameters params = validParams<GeneralUserObject>();
  MooseEnum GroupScheme("Uniform RSpace","Uniform");
  params.addRequiredParam<MooseEnum>("GroupScheme",GroupScheme, "Group method to use. Choices are: "+GroupScheme.getRawNames());
  params.addParam<Real>("dr_coef",0.2,"dr*(36*pi/Vatom)**(1./3) grouping size in r-space, eg. dr=1.0e-11m");
  params.addParam<int>("max_defect_v_size",0,"largest cluster size");
  params.addParam<int>("max_defect_i_size",0,"largest cluster size");
  params.addRequiredParam<int>("number_v","Total number of groups, count single size as a group with group size 1");
  params.addRequiredParam<int>("number_i","Total number of groups, count single size as a group with group size 1");
  params.addRequiredParam<int>("max_mobile_v", "A vector of mobile species sizes");
  params.addRequiredParam<int>("max_mobile_i", "A vector of mobile species sizes");
  params.addRequiredParam<int>("number_single_v","largest cluster size using group size of 1");
  params.addRequiredParam<int>("number_single_i","largest cluster size using group size of 1");
  params.addParam<Real>("temperature","[K], system temperature");
  params.addParam<FunctionName>("T_func","[K], system temperature as a function");
  params.addParam<bool>("update",false,"Update grouping scheme or not");
  params.addParam<UserObjectName>("material","","name of the userobject that provide material constants, i.e. emit, abosrb");
  params.addClassDescription("User object using shape functions to calculate group constants");
  return params;
}

GGroup::GGroup(const InputParameters & parameters) :
    GeneralUserObject(parameters),
    _GroupScheme(getParam<MooseEnum>("GroupScheme")),
    _dr_coef(getParam<Real>("dr_coef")),
    _Ng_v(getParam<int>("number_v")),
    _Ng_i(getParam<int>("number_i")),
    _num_v(getParam<int>("max_defect_v_size")),
    _num_i(getParam<int>("max_defect_i_size")),
    _v_size(getParam<int>("max_mobile_v")),
    _i_size(getParam<int>("max_mobile_i")),
    _single_v_group(getParam<int>("number_single_v")),
    _single_i_group(getParam<int>("number_single_i")),
    _T(isParamValid("temperature")?getParam<Real>("temperature"):0.0),
    _T_func(isParamValid("T_func")? &getFunction("T_func"):NULL),
    _update(getParam<bool>("update")),
    _has_material(getParam<UserObjectName>("material") != ""),
    _material(_has_material? &getUserObject<GMaterialConstants>("material"):NULL)
{

    _atomic_vol = _material->atomic_vol;
  // test input correctness
    if (!isParamValid("temperature") && !_T_func)
        mooseError("Temperature should be provided");
    if(_v_size > 0 && _single_v_group < _v_size)
        mooseError("max_single_group should be larger than the largest mobile size, here");
    if( _single_v_group > _Ng_v  || _single_i_group > _Ng_i)
        mooseError("max_single_group should be samller than total groups");
    if( _i_size > 0 && _single_i_group < _i_size){
        mooseError("max_single_group should be larger than the largest mobile size, there");
    }
    GroupScheme_v.reserve(_Ng_v+1);
    GroupScheme_i.reserve(_Ng_i+1);
     
    GroupScheme_v_sq = new Real[_Ng_v];
    GroupScheme_v_avg = new Real[_Ng_v];
    GroupScheme_v_del = new int[_Ng_v];
    GroupScheme_i_sq = new Real[_Ng_i];
    GroupScheme_i_avg = new Real[_Ng_i];
    GroupScheme_i_del = new int[_Ng_i];
  
    setGroupScheme();
}

GGroup::~GGroup(){
  delete[] GroupScheme_v_sq;
  delete[] GroupScheme_v_del;
  delete[] GroupScheme_i_sq;
  delete[] GroupScheme_i_del;
  delete[] GroupScheme_v_avg;
  delete[] GroupScheme_i_avg;
}

void
GGroup::initialize()
{
//print grouping info
  if (DEBUG){
    int num_groups = GroupScheme_v.size();
    for(int i=0;i<num_groups;i++){
        printf("group ID: %d; value: %d\n",i+1,GroupScheme_v[i]); 
    }
    num_groups = GroupScheme_i.size();
    for(int i=0;i<num_groups;i++){
        printf("group ID: %d; value: %d\n",-i-1,GroupScheme_i[i]); 
    }
  }
}

void
GGroup::setGroupScheme(){//total _Ng group, _Ng+1 node
  if(_GroupScheme=="Uniform"){
    if(_num_v < _Ng_v || _num_i < _Ng_i)
        mooseError("Size setting not correct");
//add vacancy group scheme
    if(_Ng_v>0){
        int single_v_group = _single_v_group+1;//1. 2. 3. each as a group, [1 2) [2 3) [3 4)
        for(int i=1;i<=single_v_group;i++){
            GroupScheme_v.push_back(i);
//            printf("add %d\n",GroupScheme_v.back());
        }
        if(_single_v_group<_Ng_v){
          double interval = 1.0*(_num_v-single_v_group)/(_Ng_v-_single_v_group);
          for(int i=1;i<(_Ng_v-_single_v_group+1);i++){
              int next_size = (int)(single_v_group+i*interval);
              GroupScheme_v.push_back(next_size);
          }
        }
        if(GroupScheme_v.back() != _num_v) GroupScheme_v[GroupScheme_v.size()-1] = _num_v;
        if((int)(GroupScheme_v.size()) != _Ng_v+1)
          mooseError("Group number " << GroupScheme_v.size() << " not correct");
        //shift to left by one (x1,x2],consistent with Golubov's paper
        for(int i=0;i<_Ng_v+1;i++)
          GroupScheme_v[i] -= 1;
    }
//append intersitial group scheme
    if(_Ng_i>0){
        int single_i_group = _single_i_group+1;//1. 2. 3. each as a group, [1 2) [2 3) [3 4)
        for(int i=1;i<=single_i_group;i++){
            GroupScheme_i.push_back(i);//make negative to distinguish from vacancy type
        }
        if(_single_i_group<_Ng_i){
          double interval = 1.0*(_num_i-single_i_group)/(_Ng_i-_single_i_group);
          for(int i=1;i<(_Ng_i-_single_i_group+1);i++){
              int next_size = (int)(single_i_group+i*interval);
              GroupScheme_i.push_back(next_size);
          }
        }
        if(GroupScheme_i.back() != _num_i) GroupScheme_i[GroupScheme_i.size()-1] = _num_i;
        if((int)(GroupScheme_i.size()) != _Ng_i+1)
          mooseError("Group number " << GroupScheme_i.size() << " not correct");
        //shift to left by one (x1,x2],consistent with Golubov's paper
        for(int i=0;i<_Ng_i+1;i++)
          GroupScheme_i[i] -= 1;
    }
  }
  else if(_GroupScheme=="RSpace"){
//add vacancy group scheme
    if(_Ng_v>0){
        //1. 2. 3. each as a group, (0,1],(1,2], (2,3]
        for(int i=0;i<=_single_v_group;i++)
            GroupScheme_v.push_back(i);
        int tmp = _single_v_group;
        while(tmp++<_Ng_v){
            int delta_size = (int)(_dr_coef*std::pow(GroupScheme_v.back(),2.0/3));
            int next_size = GroupScheme_v.back()+((delta_size>1)? delta_size:1);
            GroupScheme_v.push_back(next_size);
        }
        printf("maximum v size: %d\n",GroupScheme_v.back());
    }

    if(_Ng_i>0){
        //1. 2. 3. each as a group, (0,1],(1,2], (2,3]
        for(int i=0;i<=_single_i_group;i++)
            GroupScheme_i.push_back(i);
        int tmp = _single_i_group;
        while(tmp++<_Ng_i){
            int delta_size = (int)(_dr_coef*std::pow(GroupScheme_i.back(),2.0/3));
            int next_size = GroupScheme_i.back()+((delta_size>1)? delta_size:1);
            GroupScheme_i.push_back(next_size);
        }
        printf("maximum i size: %d\n",GroupScheme_i.back());
    } 
  }
  else
    mooseError("Group shceme: " << _GroupScheme << " not correct");


  //calculate the dispersion of each group
  int del;
  for(int i=1;i<=_Ng_v;i++){
    Real minu = 0.0, subt = 0.0;
    del = GroupScheme_v[i]-GroupScheme_v[i-1]; 
    for(int j=GroupScheme_v[i-1]+1;j<=GroupScheme_v[i];j++){
     minu += j*j;
     subt += j; 
    }
    GroupScheme_v_sq[i-1] = (minu-subt*subt/del)/del;
    GroupScheme_v_avg[i-1]= GroupScheme_v[i]-(del-1)/2.0;
    GroupScheme_v_del[i-1] = del;
    //printf("scheme: %d %f %f %d \n",GroupScheme_v[i-1],GroupScheme_v_sq[i-1],GroupScheme_v_avg[i-1],GroupScheme_v_del[i-1]);
  } 
  for(int i=1;i<=_Ng_i;i++){
    Real minu = 0.0, subt = 0.0;
    del = (GroupScheme_i[i])-(GroupScheme_i[i-1]);
    for(int j=GroupScheme_i[i-1]+1;j<=GroupScheme_i[i];j++){
     minu += j*j;
     subt += j; 
    }
    GroupScheme_i_sq[i-1]= (minu-subt*subt/del)/del;
    GroupScheme_i_avg[i-1]= GroupScheme_i[i]-(del-1)/2.0;
    GroupScheme_i_del[i-1] = del;
  } 

}

void
GGroup::updateGroupScheme(){
//adaptively update scheme based on the distribution profile
    GroupScheme_v.clear();
    GroupScheme_i.clear();

    setGroupScheme();//change to new one
}

void
GGroup::execute()
{
  if(_update){
    updateGroupScheme();
  }
}

void GGroup::finalize()
{}



Real
GGroup::_emit(int clustersize) const //[cr_start,cr_end)
{
  const char* species = (clustersize>0)?"V":"I";
  int tagi = 0;//denote mobility
  if(clustersize>0){
      if(clustersize<=_v_size)
          tagi = 1;
  }
  else{
      if(-clustersize<=_i_size)
          tagi = 1;
  }
  Real T = _T_func? _T_func->value(_t,dummy):_T;
  Real val = _material->emit((int)std::abs(clustersize),1,T,species,species,tagi,1);
  //printf("emit of clustersize (%d): %f\n",clustersize,val);
  return val;
}

Real
GGroup::_disl(int clustersize) const //[cr_start,cr_end)
{
  const char* species = (clustersize>0)?"V":"I";
  int tagi = 0;//denote mobility
  Real val = 0.0;
  Real T = _T_func? _T_func->value(_t,dummy):_T;
  if(clustersize>0){
    if(clustersize>_v_size) return 0.0;
    tagi = 1;
    val = _material->disl_ksq(clustersize,species,T,tagi);
  }
  else{
    if(-clustersize>_i_size) return 0.0;
    tagi = 1;
    val = _material->disl_ksq(-clustersize,species,T,tagi);
  }
  //printf("dislocation of clustersize (%d): %f\n",clustersize,val);
  return val;
}

Real
GGroup::_diff(int clustersize) const //[cr_start,cr_end)
{
  const char* species = (clustersize>0)?"V":"I";
  int tagi = 0;//denote mobility
  Real val = 0.0;
  Real T = _T_func? _T_func->value(_t,dummy):_T;
  if(clustersize>0){
    if(clustersize>_v_size) return 0.0;
    tagi = 1;
    val = _material->diff(clustersize,species,T);
  }
  else{
    if(-clustersize>_i_size) return 0.0;
    tagi = 1;
    val = _material->diff(clustersize,species,T);
  }
  //printf("diffusion of clustersize (%d): %f\n",clustersize,val);
  return val;
}

Real
GGroup::_absorb(int clustersize1, int clustersize2) const //[ot_start,ot_end),[cr_start,cr_end)
{
  Real val = 0.0;
  Real T = _T_func? _T_func->value(_t,dummy):_T;
  int i = std::abs(clustersize1);
  int j = std::abs(clustersize2);
  int tagi = 0,tagj = 0;//denote mobility: 0, imobile, 1, mobile
  int flag;
  if(clustersize1>0 && clustersize2>0){//vv 
    if(i<=_v_size)
        tagi = 1;
    if(j<=_v_size)
        tagj = 1;
    //printf("absorb v %d with v %d: %lf \n",i,j,_material->absorb(i,j,"V","V",_T,tagi,tagj));//,absorb(i,j,"V","V",_T,tagi,tagj));
    //return _material->absorb(i,j,"V","V",_T,tagi,tagj);//absorption between i and j
    flag = tagi + 2*tagj; 
    return _material->absorbVV(i,j,flag,T);//flag=0: i j immobile; flag=1: i mobile; flag=2: j mobile; flag=3: i j mobile
  }
  else if(clustersize1>0 && clustersize2<0){//vi
    if(i<=_v_size)
        tagi = 1;
    if(j<=_i_size)
        tagj = 1;
    //printf("absorb v %d with i %d: %lf \n",i,j,_material->absorb(i,j,"V","I",_T,tagi,tagj));//,absorb(i,j,"V","I",_T,tagi,tagj));
    //return _material->absorb(i,j,"V","I",_T,tagi,tagj);//absorption between i and j
    flag = tagi + 2*tagj; 
    return _material->absorbVI(i,j,flag,T);//flag=0: i j immobile; flag=1: i mobile; flag=2: j mobile; flag=3: i j mobile
  }
  else if(clustersize1<0 && clustersize2>0){//iv
    if(i<=_i_size)
        tagi = 1;
    if(j<=_v_size)
        tagj = 1;
    //printf("absorb i %d with v %d: %lf\n",i,j,_material->absorb(i,j,"I","V",_T,tagi,tagj));
    //return _material->absorb(i,j,"I","V",_T,tagi,tagj);//absorption between i and j
    flag = tagj + 2*tagi; 
    return _material->absorbVI(j,i,flag,T);//flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile
  }
  else{//ii
    if(i<=_i_size)
        tagi = 1;
    if(j<=_i_size)
        tagj = 1;
    //printf("absorb i %d with i %d: %lf\n",i,j,_material->absorb(i,j,"I","I",_T,tagi,tagj));
    //return _material->absorb(i,j,"I","I",_T,tagi,tagj);//absorption between i and j
    flag = tagi + 2*tagj; 
    return _material->absorbII(i,j,flag,T);//flag=0: i j immobile; flag=1: i mobile; flag=2: j mobile; flag=3: i j mobile
  }
}


int
GGroup::CurrentGroupV(int i) const{
    std::vector<int>::const_iterator it=std::lower_bound(GroupScheme_v.begin(),GroupScheme_v.end(),i); 
    return it-GroupScheme_v.begin();
}

int
GGroup::CurrentGroupI(int i) const{
    std::vector<int>::const_iterator it=std::lower_bound(GroupScheme_i.begin(),GroupScheme_i.end(),i); 
    return it-GroupScheme_i.begin();
}


