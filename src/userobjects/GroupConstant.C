/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

//*******************make pair to loop up all the group constants   '+': vacancy; '-': intersitial************************//
// calculate group constant based on hypothetical shape functions

#include "GroupConstant.h"
#include<math.h>
#include<algorithm>

template<>
InputParameters validParams<GroupConstant>()
{
  InputParameters params = validParams<GeneralUserObject>();
  MooseEnum GroupScheme("Uniform","Uniform");
  params.addRequiredParam<MooseEnum>("GroupScheme",GroupScheme, "Group method to use. Choices are: "+GroupScheme.getRawNames());
  params.addRequiredParam<int>("number_v","Total number of groups, count single size as a group with group size 1");
  params.addRequiredParam<int>("number_i","Total number of groups, count single size as a group with group size 1");
  params.addParam<std::vector<int> >("mobile_v_size", "A vector of mobile species sizes");
  params.addParam<std::vector<int> >("mobile_i_size", "A vector of mobile species sizes");
  params.addParam<int>("max_mobile_v", "maximum size of mobile vacancy cluster");
  params.addParam<int>("max_mobile_i", "maximum size of mobile intersitial cluster");
  params.addRequiredParam<Real>("temperature","[K], system temperature");
  params.addParam<bool>("update",false,"Update grouping scheme or not");
  params.addRequiredParam<UserObjectName>("material","name of the userobject that provide material constants, i.e. emit, abosrb");
  params.addClassDescription("User object using shape functions to calculate group constants");
  return params;
}

GroupConstant::GroupConstant(const InputParameters & parameters) :
    GeneralUserObject(parameters),
    _GroupScheme(getParam<MooseEnum>("GroupScheme")),
    _Ng_v(isParamValid("number_v")?getParam<int>("number_v"):(getParam<int>("max_defect_v_size")-1)),
    _Ng_i(isParamValid("number_i")?getParam<int>("number_i"):(getParam<int>("max_defect_i_size")-1)),
    _T(getParam<Real>("temperature")),
    _update(getParam<bool>("update")),
    _material(&getUserObject<GMaterialConstants>("material"))
   // _emit_array(NULL),
   // _absorb_matrix(NULL)
{
    if ((isParamValid("mobile_v_size")^isParamValid("max_mobile_v") != 1) || (isParamValid("mobile_v_size")^isParamValid("max_mobile_v") != 1)){
      mooseError("Check the setting of mobile species, only one definition is allowed");
    }
    std::vector<int> v_size,i_size;
    if (isParamValid("mobile_v_size"))
      _v_size = getParam<std::vector<int> >("mobile_v_size");
    else {
      int max_mobile_v = getParam<int>("max_mobile_v");
      for(int i=0;i<max_mobile_v;i++)
        _v_size.push_back(i+1);
    }
    if (isParamValid("mobile_i_size"))
      _i_size = getParam<std::vector<int> >("mobile_i_size");
    else {
      int max_mobile_i = getParam<int>("max_mobile_i");
      for(int i=0;i<max_mobile_i;i++)
        _i_size.push_back(i+1);
    }

    _atomic_vol = _material->atomic_vol;
  // test input correctness
    if(_v_size.size() > 0 && _Ng_v < *std::max_element(_v_size.begin(),_v_size.end()))
        mooseError("max_single_group should be larger than the largest mobile size, here");
    if( _i_size.size() > 0 && _Ng_i < *std::max_element(_i_size.begin(),_i_size.end()))
        mooseError("max_single_group should be larger than the largest mobile size, there");

    GroupScheme_v.reserve(_Ng_v+1);
    GroupScheme_i.reserve(_Ng_i+1);
}

GroupConstant::~GroupConstant(){
  //  delete [] _emit_array;
  //  _emit_array = NULL;
   // for(int i=0;i<_Ng_v+_Ng_i;i++)
   //     delete [] _absorb_matrix[i];
   // delete [] _absorb_matrix;
   // _absorb_matrix = NULL;
/*
    if(_Ng_i>0){
        delete [] _emit_i_array;
        for(int i=0;i<_Ng_i;i++)
            delete [] _absorb_ii_matrix[i];
        for(int i=0;i<_Ng_v;i++)
            delete [] _absorb_vi_matrix[i];
        delete [] _absorb_vi_matrix;
        _emit_i_array = NULL;
        _absorb_ii_matrix = NULL;
    }
 */   
}

void
GroupConstant::initialize()
{
  setGroupScheme();
  setGroupConstant();
//print grouping info
/*
  int num_groups = GroupScheme_v.size();
  for(int i=0;i<num_groups;i++){
      printf("group ID: %d; value: %d\n",i+1,GroupScheme_v[i]); 
  }
  num_groups = GroupScheme_i.size();
  for(int i=0;i<num_groups;i++){
      printf("group ID: %d; value: %d\n",i+1,-GroupScheme_i[i]); 
  }
//print group constants
  for(SingleKey::const_iterator it=_emit_array.begin(); it!=_emit_array.end();it++){
      printf("group: %d; emission: %10.8e\n",it->first,it->second); 
  }
  for(DoubleKey::const_iterator it=_absorb_matrix.begin();it!=_absorb_matrix.end();it++){
      printf("group1: %d; group2: %d; absorption: %10.8e\n",it->first.first,it->first.second,it->second);
  }
*/
}

void
GroupConstant::setGroupScheme(){//total _Ng group, _Ng+1 node
  if(_GroupScheme=="Uniform"){
//add vacancy group scheme
    if(_Ng_v>0){
        //1. 2. 3. each as a group, [1 2) [2 3) [3 4)
        for(int i=1;i<=_Ng_v+1;i++){
            GroupScheme_v.push_back(i);
//            printf("add %d\n",GroupScheme_v.back());
        }
        if((int)(GroupScheme_v.size()) != _Ng_v+1)
          mooseError("Group number ", GroupScheme_v.size(), " not correct");
    }
//append intersitial group scheme
    if(_Ng_i>0){
        //1. 2. 3. each as a group, [1 2) [2 3) [3 4)
        for(int i=1;i<=_Ng_i+1;i++){
            GroupScheme_i.push_back(-i);//make negative to distinguish from vacancy type
 //           printf("add %d\n",-i);
        }
        if((int)(GroupScheme_i.size()) != _Ng_i+1)
          mooseError("Group number ", GroupScheme_i.size(), " not correct");
    }
  }
  else
    mooseError("Group shceme: ", _GroupScheme, " not correct");
  

}

void
GroupConstant::updateGroupScheme(){
//adaptively update scheme based on the distribution profile
    GroupScheme_v.clear();
    GroupScheme_i.clear();
    setGroupScheme();//change to new one
}

void
GroupConstant::setGroupConstant(){
  Real val = 0.0;

//emission coefs
  for(int i=1;i<=_Ng_v;i++)
    _emit_array.insert(std::make_pair(i,emit_gc(GroupScheme_v[i-1],GroupScheme_v[i])));
  for(int i=1;i<=_Ng_i;i++)
    _emit_array.insert(std::make_pair(-i,emit_gc(GroupScheme_i[i-1],GroupScheme_i[i])));

  
  int total_mobile_v =(int)_v_size.size();
  int total_mobile_i =(int)_i_size.size();

//dislocation absorption coefs
  for(int i=1;i<=total_mobile_v;i++)
    _disl_array.insert(std::make_pair(i,disl_gc(GroupScheme_v[i-1],GroupScheme_v[i])));
  for(int i=1;i<=total_mobile_i;i++)
    _disl_array.insert(std::make_pair(-i,disl_gc(GroupScheme_i[i-1],GroupScheme_i[i])));

//diffusion coefs
  for(int i=1;i<=total_mobile_v;i++)
    _diff_array.insert(std::make_pair(i,disl_gc(GroupScheme_v[i-1],GroupScheme_v[i])));
  for(int i=1;i<=total_mobile_i;i++)
    _diff_array.insert(std::make_pair(-i,disl_gc(GroupScheme_i[i-1],GroupScheme_i[i])));

//absorption coefficients
  for(int i=1;i<=_Ng_v;i++){
    for(int j=1;j<=total_mobile_v;j++){
      val = absorb_gc(GroupScheme_v[i-1],GroupScheme_v[i],GroupScheme_v[j-1],GroupScheme_v[j]);
      _absorb_matrix.insert(std::make_pair(std::make_pair(i,j),val));
    }
    for(int j=1;j<=total_mobile_i;j++){
      val = absorb_gc(GroupScheme_v[i-1],GroupScheme_v[i],GroupScheme_i[j-1],GroupScheme_i[j]);
      _absorb_matrix.insert(std::make_pair(std::make_pair(i,-j),val));
    }
  }

  for(int i=1;i<=_Ng_i;i++){
    for(int j=1;j<=total_mobile_v;j++){
      val= absorb_gc(GroupScheme_i[i-1],GroupScheme_i[i],GroupScheme_v[j-1],GroupScheme_v[j]);
     // printf("iv: %d %d %f\n",i,j,val);
      _absorb_matrix.insert(std::make_pair(std::make_pair(-i,j),val));
    }
    for(int j=1;j<=total_mobile_i;j++){
      val = absorb_gc(GroupScheme_i[i-1],GroupScheme_i[i],GroupScheme_i[j-1],GroupScheme_i[j]);
     // printf("ii: %d %d %f\n",i,j,val);
      _absorb_matrix.insert(std::make_pair(std::make_pair(-i,-j),val));
    }
  }


  //print out coefficients for test
/*
  for(int i=0;i<_Ng_v;i++)
      printf("emission from size %d : %f\n",i+1,_emit_array[i]);
  for(int i=0;i<_Ng_i;i++)
      printf("emission from size %d : %f\n",-i-1,_emit_array[-i]);
  */
//  for(DoubleKey::iterator it=_absorb_matrix.begin();it!=_absorb_matrix.end();it++){
//      std::cout << (it->first).first << " "<< (it->first).second << " "  << it->second << std::endl;
//  }
}

Real
GroupConstant::emit_gc(int cr_start, int cr_end) //[cr_start,cr_end)
{
  int pos_start = abs(cr_start);
  int pos_end = abs(cr_end);
  
  //std::cout<< "max long int: " << std::numeric_limits<long int>::max() << std::endl;
  
  int multiply = ((cr_start>0)&&(cr_end>0)) || ((cr_start<0)&&(cr_end<0));
  if(!(pos_start<pos_end) || multiply <= 0){
      //printf("what %d %d %d %d %d\n", cr_start,cr_end,pos_start,pos_end,multiply);
      mooseError("Wrong interval terminal point in emit_gc function");
  }
  
  const char* species = (cr_start>0)?"V":"I";
  int tagi = 0;//denote mobility
  Real sum = 0.0;
  Real counter = 0.0;
  for(int j=pos_start;j<pos_end;j++){//current group
      tagi = 0;
      if(cr_start>0){
          if(std::find(_v_size.begin(),_v_size.end(),j) != _v_size.end())
              tagi = 1;
      }
      else{
          if(std::find(_i_size.begin(),_i_size.end(),j) != _i_size.end())
              tagi = 1;
      }
      sum += _material->emit(j,1,_T,species,species,tagi,1);//emission of vacancy from vacancy cluster, from function in IronProperty.C
      counter += 1.0;
  }
  //printf("start, end, value: %d %d %f\n",cr_start,cr_end,sum/counter);
  return sum/counter;
}

Real
GroupConstant::disl_gc(int cr_start, int cr_end) //[cr_start,cr_end)
{
  int pos_start = abs(cr_start);
  int pos_end = abs(cr_end);
  //double multiply = (double)cr_start * (double)cr_end;
  int multiply = ((cr_start>0)&&(cr_end>0)) || ((cr_start<0)&&(cr_end<0));
  if(!(pos_start<pos_end) || multiply <= 0)
      mooseError("Wrong interval terminal point in disl_gc function");
  
  const char* species = (cr_start>0)?"V":"I";
  int tagi = 0;//denote mobility
  Real sum = 0.0;
  Real counter = 0.0;
  for(int j=pos_start;j<pos_end;j++){//current group
      tagi = 0;
      if(cr_start>0){
          if(std::find(_v_size.begin(),_v_size.end(),j) != _v_size.end())
              tagi = 1;
      }
      else{
          if(std::find(_i_size.begin(),_i_size.end(),j) != _i_size.end())
              tagi = 1;
      }
      sum += _material->disl_ksq(j,species,_T,tagi);
      counter += 1.0;
  }
  return sum/counter;
}

Real
GroupConstant::diff_gc(int cr_start, int cr_end) //[cr_start,cr_end)
{
  int pos_start = abs(cr_start);
  int pos_end = abs(cr_end);
  //double multiply = (double)cr_start * (double)cr_end;
  int multiply = ((cr_start>0)&&(cr_end>0)) || ((cr_start<0)&&(cr_end<0));
  if(!(pos_start<pos_end) || multiply <= 0)
      mooseError("Wrong interval terminal point in diff_gc function");
  
  const char* species = (cr_start>0)?"V":"I";
  int tagi = 0;//denote mobility
  Real sum = 0.0;
  Real counter = 0.0;
  for(int j=pos_start;j<pos_end;j++){//current group
      tagi = 0;
      if(cr_start>0){
          if(std::find(_v_size.begin(),_v_size.end(),j) != _v_size.end())
              tagi = 1;
      }
      else{
          if(std::find(_i_size.begin(),_i_size.end(),j) != _i_size.end())
              tagi = 1;
      }
      sum += _material->diff(j,species,_T)*tagi;
      counter += 1.0;
  }
  return sum/counter;
}

Real
GroupConstant::absorb_gc(int ot_start, int ot_end, int cr_start, int cr_end) //[ot_start,ot_end),[cr_start,cr_end)
{
  //[cr_start,cr_end) should be size 1 and mobile, i.e. cr_end=cr_start+1
  int pos_ot_start = abs(ot_start);
  int pos_ot_end = abs(ot_end);
  int pos_cr_start = abs(cr_start);
  int pos_cr_end = abs(cr_end);

  //printf("start end: %d %d %d %d\n",ot_start,ot_end,cr_start,cr_end);
  if(pos_cr_end != pos_cr_start+1)
      mooseError("Wrong interval terminal point in absorb_gc mobile range");

  if(!(pos_ot_start<pos_ot_end) || !(pos_cr_start<pos_cr_end) || (double)ot_start*(double)ot_end<=0 || (double)cr_start*(double)cr_end<=0)
      mooseError("Wrong interval terminal point in absorb_gc function");
   
  Real sum = 0.0;
  Real counter = 0.0;
  int tagi = 0,tagj = 0;//denote mobility: 0, imobile, 1, mobile
  int total_mobile_v = (int)_v_size.size();
  int total_mobile_i = (int)_i_size.size();

 // if(ot_start >= _single_v_group && cr_start >= _single_v_group) return 0.0;
  if(ot_start>0 && cr_start>0){//vv reaction
      for(int i=pos_ot_start;i<pos_ot_end;i++){//other group
          tagi = 0;
          if(std::find(_v_size.begin(),_v_size.end(),i) != _v_size.end())
              tagi = 1;
          for(int j=pos_cr_start;j<pos_cr_end;j++){//current group//redundant
              tagj = 0;
              if(std::find(_v_size.begin(),_v_size.end(),j) != _v_size.end())
                  tagj = 1;
              sum += _material->absorb(i,j,"V","V",_T,tagi,tagj);//absorption between i and j
              counter += 1.0;
              //printf("absorb v %d with v %d: %lf %lf\n",i,j,_material->absorb(i,j,"V","V",_T,tagi,tagj));//,absorb(i,j,"V","V",_T,tagi,tagj));
          }
      }
  }
  else if(ot_start>0 && cr_start<0){//vi reaction
      for(int i=pos_ot_start;i<pos_ot_end;i++){//other group
          tagi = 0;
          if(std::find(_v_size.begin(),_v_size.end(),i) != _v_size.end())
              tagi = 1;
          for(int j=pos_cr_start;j<pos_cr_end;j++){//current group//redundant
              tagj = 0;
              if(std::find(_i_size.begin(),_i_size.end(),j) != _i_size.end())
                  tagj = 1;
              sum += _material->absorb(i,j,"V","I",_T,tagi,tagj);//absorption between i and j
              counter += 1.0;
              //printf("absorb v %d with i %d: %lf %lf\n",i,j,_material->absorb(i,j,"V","I",_T,tagi,tagj),absorb(i,j,"V","I",_T,tagi,tagj));
          }
      }
  }
  else if(ot_start<0 && cr_start>0){//iv reaction
      for(int i=pos_ot_start;i<pos_ot_end;i++){//other group
          tagi = 0;
          if(std::find(_i_size.begin(),_i_size.end(),i) != _i_size.end())
              tagi = 1;
          for(int j=pos_cr_start;j<pos_cr_end;j++){//current group//redundant
              tagj = 0;
              if(std::find(_v_size.begin(),_v_size.end(),j) != _v_size.end())
                  tagj = 1;
              sum += _material->absorb(i,j,"I","V",_T,tagi,tagj);//absorption between i and j
              counter += 1.0;
              //printf("absorb i %d with v %d: %lf\n",i,j,absorb(i,j,"I","V",_T,tagi,tagj));
          }
      }
  }
  else{//ii reaction
      for(int i=pos_ot_start;i<pos_ot_end;i++){//other group
          tagi = 0;
          if(std::find(_i_size.begin(),_i_size.end(),i) != _i_size.end())
              tagi = 1;
          for(int j=pos_cr_start;j<pos_cr_end;j++){//current group//redundant
              tagj = 0;
              if(std::find(_i_size.begin(),_i_size.end(),j) != _i_size.end())
                  tagj = 1;
              sum += _material->absorb(i,j,"I","I",_T,tagi,tagj);//absorption between i and j
              counter += 1.0;
              //printf("absorb i %d with i %d: %lf\n",i,j,absorb(i,j,"I","I",_T,tagi,tagj));
          }
      }
  }
  //printf("sum: %f counter: %d\n",sum,counter);
  return sum/counter;
}

void
GroupConstant::execute()
{
  if(_update){
    updateGroupScheme();
    setGroupConstant();
  }
}

void GroupConstant::finalize()
{}



Real
GroupConstant::_emit(int groupid) const //[cr_start,cr_end)
{
    SingleKey::const_iterator it = _emit_array.find(groupid);
    if(it != _emit_array.end()) 
        return it->second;
    return 0.0;
}

Real
GroupConstant::_disl(int groupid) const //[cr_start,cr_end)
{
    SingleKey::const_iterator it = _disl_array.find(groupid);
    if(it != _disl_array.end()) 
        return it->second;
    return 0.0;
}

Real
GroupConstant::_diff(int groupid) const //[cr_start,cr_end)
{
    SingleKey::const_iterator it = _diff_array.find(groupid);
    if(it != _diff_array.end()) 
        return it->second;
    return 0.0;
}

Real
GroupConstant::_absorb(int groupid1, int groupid2) const //[ot_start,ot_end),[cr_start,cr_end)
{
    DoubleKey::const_iterator it = _absorb_matrix.find(std::make_pair(groupid1,groupid2));
    if(it != _absorb_matrix.end())
        return it->second;//_absorb_matrix[std::make_pair(groupid1,groupid2)]; 
    it = _absorb_matrix.find(std::make_pair(groupid2,groupid1));
    if(it != _absorb_matrix.end())
        return it->second;//_absorb_matrix[std::make_pair(groupid2,groupid1)]; 

    printf("Absorption constant not found and return 0: %d vs %d\n",groupid1,groupid2);
    return 0.0;
}


