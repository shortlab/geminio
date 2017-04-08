/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ImmobileDefects.h"
#include "Conversion.h"

template<>
InputParameters validParams<ImmobileDefects>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<int>("number_v","Maximum vacancy cluster size");
  params.addRequiredParam<int>("number_i","Maximum interstitial cluster size");
  params.addCoupledVar("coupled_v_vars","coupled vacancy type variables");
  params.addCoupledVar("coupled_i_vars","coupled intersitial type variables");
  params.addRequiredParam<std::vector<int> >("mobile_v_size", "A vector of mobile species, sorted by custom action");
  params.addRequiredParam<std::vector<int> >("mobile_i_size", "A vector of mobile species, sorted by custom action");
  params.addRequiredParam<UserObjectName>("user_object","The name of user object providing interaction constants");
  return params;
}

ImmobileDefects::ImmobileDefects(const InputParameters & parameters)
     :Kernel(parameters),
     _number_v(getParam<int>("number_v")),
     _number_i(getParam<int>("number_i")),
     _v_size(getParam<std::vector<int> >("mobile_v_size")),
     _i_size(getParam<std::vector<int> >("mobile_i_size")),
     _gc(getUserObject<GroupConstant>("user_object"))
{
  NonlinearVariableName cur_var_name = getParam<NonlinearVariableName>("variable");
  _cur_size = getGroupNumber(cur_var_name);
  if(_cur_size>0 && _cur_size< _v_size.size()){
    mooseError("Please check the mobile vacancy cluster size range, should be immobile>mobile");
  }
  else if(_cur_size<0 && -_cur_size<_i_size.size()){
    mooseError("Please check the mobile intersitial cluster size range, should be immobile>mobile");
  }
  unsigned int num_v_coupled = coupledComponents("coupled_v_vars");
  unsigned int num_i_coupled = coupledComponents("coupled_i_vars");
  if(num_v_coupled>0){
    _no_v_vars.resize(num_v_coupled);
    _val_v_vars.resize(num_v_coupled);
  }
  if(num_i_coupled>0){
    _no_i_vars.resize(num_i_coupled);
    _val_i_vars.resize(num_i_coupled);
  }
 

  for (unsigned int i=0; i < num_v_coupled; ++i){
    _no_v_vars[i] = coupled("coupled_v_vars",i);
    _val_v_vars[i] = &coupledValue("coupled_v_vars",i);
  }
  for (unsigned int i=0; i < num_i_coupled; ++i){
    _no_i_vars[i] = coupled("coupled_i_vars",i);
    _val_i_vars[i] = &coupledValue("coupled_i_vars",i);
  }
    
  //printf("immobile constructed: %s\n",cur_var_name.c_str());
}

Real
ImmobileDefects::computeQpResidual()
{
  int v_size = _v_size.size();
  int i_size = _i_size.size();
  int cur_size;
  Real res_sum = 0.0;
  //printf("return immobile initial: %f %d\n",res_sum,_cur_size);     

  if(_cur_size>0){//v type
    cur_size = _cur_size;

    //vv reaction loss(-)
    int tmp_size = std::min(v_size,_number_v-cur_size);
    for(int i=0;i<tmp_size;i++){
      res_sum += (*_val_v_vars[i])[_qp]*_u[_qp]*_gc._absorb(_v_size[i],cur_size);
      //printf("reaction %d (-): %d %d\n",cur_size,_v_size[i],cur_size);
    }

    //vi reaction loss(-)
    for(int i=0;i<i_size;i++){
      res_sum += (*_val_i_vars[i])[_qp]*_u[_qp]*_gc._absorb(-_i_size[i],cur_size);
    }

    //vv reaction gain(+)
    int tmp = std::min((int)(cur_size/2),v_size);//_v_size.back()==v_size, avoid double count
    for(int i=0;i<tmp;i++){
        res_sum -= (*_val_v_vars[v_size+i])[_qp]*(*_val_v_vars[i])[_qp]*_gc._absorb(_v_size[i],cur_size-_v_size[i]);
        //printf("reaction %d (+): %d %d\n",cur_size,_v_size[i],cur_size-_v_size[i]);
    }

    //vi reaction gain(+)
    if(i_size!=0 && cur_size != _number_v){
      tmp = _no_v_vars.size()-v_size-v_size;
      for(int i=0;i<tmp;i++){
      //  printf("return imm : %d %d %d\n",tmp,v_size,_val_v_vars.size());     
        res_sum -= (*_val_v_vars[v_size+v_size+i])[_qp]*(*_val_i_vars[i])[_qp]*_gc._absorb(-_i_size[i],cur_size+_i_size[i]);
      }
    }

    //v+1 emission gain(+)  
    if(cur_size<_number_v){
      res_sum -= (*_val_v_vars[v_size+v_size])[_qp]*_gc._emit(cur_size+1);    
      //printf("emission %d (+): %d\n",cur_size,cur_size+1);
    }

    //v emission loss(-)
    res_sum += _u[_qp]*_gc._emit(cur_size);
    //printf("emission %d (-): %d\n",cur_size,cur_size);
  
  }

  else{
    cur_size = -_cur_size;

    //ii reaction loss(-)
    int tmp_size = std::min(i_size,_number_i-cur_size);
    for(int i=0;i<tmp_size;i++){
      res_sum += (*_val_i_vars[i])[_qp]*_u[_qp]*_gc._absorb(-_i_size[i],-cur_size);
      //printf("reaction %d (-): %d %d\n",_cur_size,_cur_size,-_i_size[i]);     
    }

    //iv reaction loss(-)
    for(int i=0;i<v_size;i++){
      res_sum += (*_val_v_vars[i])[_qp]*_u[_qp]*_gc._absorb(_v_size[i],-cur_size);
      //printf("reaction %d (-): %d %d\n",_cur_size,_cur_size,_v_size[i]);     
    }

    //ii reaction gain(+)
    int tmp = std::min((int)(cur_size/2),i_size);//_i_size.back()==i_size, avoid double count
    for(int i=0;i<tmp;i++){
        res_sum -= (*_val_i_vars[i_size+i])[_qp]*(*_val_i_vars[i])[_qp]*_gc._absorb(-_i_size[i],-cur_size+_i_size[i]);
        //printf("reaction %d (+): %d %d\n",_cur_size,_i_size[i]-cur_size,-_i_size[i]);     
    }

    //iv reaction gain(+)
    if(v_size!=0 && cur_size != _number_i){
      tmp = _no_i_vars.size()-i_size-i_size;
      for(int i=0;i<tmp;i++){
        res_sum -= (*_val_i_vars[i_size+i_size+i])[_qp]*(*_val_v_vars[i])[_qp]*_gc._absorb(_v_size[i],-(cur_size+_v_size[i]));
        //printf("reaction %d (+): %d %d\n",_cur_size,-_v_size[i]-cur_size,_v_size[i]);     
      }
    }

    //i+1 emission gain(+)  
    if(cur_size<_number_i){
      res_sum -= (*_val_i_vars[i_size+i_size])[_qp]*_gc._emit(-cur_size-1);    
      //printf("emission %d (+): %d\n",_cur_size,-cur_size-1);     
    }

    //i emission loss(-)
    res_sum += _u[_qp]*_gc._emit(-cur_size);
    //printf("emission %d (-): %d\n",_cur_size,-cur_size);     
  }

  return res_sum *_test[_i][_qp];
}

Real
ImmobileDefects::computeQpJacobian()
{
  int v_size = _v_size.size();
  int i_size = _i_size.size();
  int cur_size;
  Real jac_sum = 0.0;

  if(_cur_size>0){//v type
    cur_size = _cur_size;

    //vv reaction loss(-)
    int tmp_size = std::min(v_size,_number_v-cur_size);
    for(int i=0;i<tmp_size;i++){
      jac_sum += (*_val_v_vars[i])[_qp]*_gc._absorb(_v_size[i],cur_size);
    }

    //vi reaction loss(-)
    for(int i=0;i<i_size;i++){
      jac_sum += (*_val_i_vars[i])[_qp]*_gc._absorb(-_i_size[i],cur_size);
    }

    //v emission loss(-)
    jac_sum += _gc._emit(cur_size);
  
  }

  else{
    cur_size = -_cur_size;
    //ii reaction loss(-)
    int tmp_size = std::min(i_size,_number_i-cur_size);
    for(int i=0;i<tmp_size;i++){
      jac_sum += (*_val_i_vars[i])[_qp]*_gc._absorb(-_i_size[i],-cur_size);
    }

    //iv reaction loss(-)
    for(int i=0;i<v_size;i++){
      jac_sum += (*_val_v_vars[i])[_qp]*_gc._absorb(_v_size[i],-cur_size);
    }

    //i emission loss(-)
    jac_sum += _gc._emit(-cur_size);
  }

  return jac_sum *_test[_i][_qp]*_phi[_j][_qp];
}

Real 
ImmobileDefects::computeQpOffDiagJacobian(unsigned int jvar){
      return 0.0;
}


int
ImmobileDefects::getGroupNumber(std::string str)
{
  int len=str.length(),i=len;
  while(std::isdigit(str[i-1])) i--;
  int no = std::atoi((str.substr(i)).c_str());
  while(i>=0){
      i--;
      if(str[i]=='v'){no = no;break;}
      if(str[i]=='i'){no = -no;break;}
  }
  return no;
}
