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

#include "MobileDefects.h"
#include "Conversion.h"

template<>
InputParameters validParams<MobileDefects>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<int>("number_v","Maximum vacancy cluster size");
  params.addRequiredParam<int>("number_i","Maximum interstitial cluster size");
  params.addCoupledVar("coupled_v_vars","coupled vacancy type variables");
  params.addCoupledVar("coupled_i_vars","coupled intersitial type variables");
  params.addRequiredParam<std::vector<int> >("mobile_v_size", "A vector of mobile species");
  params.addRequiredParam<std::vector<int> >("mobile_i_size", "A vector of mobile species");
  params.addRequiredParam<UserObjectName>("user_object","The name of user object providing interaction constants");
  return params;
}

MobileDefects::MobileDefects(const InputParameters & parameters)
     :Kernel(parameters),
     _number_v(getParam<int>("number_v")),
     _number_i(getParam<int>("number_i")),
     _v_size(getParam<std::vector<int> >("mobile_v_size")),
     _i_size(getParam<std::vector<int> >("mobile_i_size")),
     _gc(getUserObject<GroupConstant>("user_object"))
{
   int nvcoupled = coupledComponents("coupled_v_vars");
   int nicoupled = coupledComponents("coupled_i_vars");
  if(_number_v>0){
    _no_v_vars.resize(nvcoupled);
    _val_v_vars.resize(nvcoupled);
  }
  if(_number_i>0){
    _no_i_vars.resize(nicoupled);
    _val_i_vars.resize(nicoupled);
  }

  std::string var_name;
  for (int i=0; i < nvcoupled; ++i)
  {
 //  printf("current v: %d\n",i);
    _no_v_vars[i] = coupled("coupled_v_vars",i);
    _val_v_vars[i] = &coupledValue("coupled_v_vars",i);
  }
  for (int i=0; i < nicoupled; ++i)
  {
//   printf("current i: %d\n",i);
    _no_i_vars[i] = coupled("coupled_i_vars",i);
    _val_i_vars[i] = &coupledValue("coupled_i_vars",i);
  }
  NonlinearVariableName cur_var_name = getParam<NonlinearVariableName>("variable");
  _cur_size = getGroupNumber(cur_var_name.c_str());
  
  //printf("mobile constructed: %s\n",cur_var_name.c_str());
}

Real
MobileDefects::computeQpResidual()
{
  Real res_sum = 0.0;
  int cur_size;//should be positive value
  int ii = _i_size.size();
  int vv = _v_size.size();
  int max_vi;

  //printf("return mobile initial: %f %d\n",res_sum,_cur_size);     

  if(_cur_size>0){//v type
    cur_size = _cur_size; 
    max_vi =  std::min(_cur_size+ii,_number_v);

    //vi reaction loss(-)
    for(int i=1;i<=_number_i;i++){
      res_sum += (*_val_i_vars[i-1])[_qp]*_u[_qp]*_gc._absorb(cur_size,-i);
    }

    //vv reaction loss(-)
    for(int i=1;i <= _number_v-cur_size;i++){//garantee the largest size doesn't exceed _number_v
      //printf("reaction %d (-): %d %d\n",cur_size,cur_size,i);     
      res_sum += (*_val_v_vars[i-1])[_qp]*_u[_qp]*_gc._absorb(cur_size,i);
    }
    if(cur_size*2<=_number_v){
      //printf("reaction %d (-): %d %d\n",cur_size,cur_size,cur_size);     
      res_sum += (*_val_v_vars[cur_size-1])[_qp]*_u[_qp]*_gc._absorb(cur_size,cur_size);
    }
      
  
    //vv reaction gain(+)
    for(int i=1;i <= (int)(cur_size/2);i++){
      //if(std::find(_v_size.begin(),_v_size.end(),cur_size-i) != _v_size.end() || std::find(_v_size.begin(),_v_size.end(),i) != _v_size.end()){//TODO: can optimize
        //printf("reaction %d (+): %d %d\n",cur_size,cur_size-i,cur_size);     
        res_sum -= (*_val_v_vars[cur_size-i-1])[_qp]*(*_val_v_vars[i-1])[_qp]*_gc._absorb(cur_size-i,i);
     // }
    }
  
    //vi reaction gain(+)
    for(int i=cur_size+1;i<=max_vi;i++){
      //if(std::find(_i_size.begin(),_i_size.end(),i-cur_size) != _i_size.end() || std::find(_v_size.begin(),_v_size.end(),i) != _v_size.end()){
      if(i-cur_size <= ii || i <= vv ){//make sure one is mobile
        res_sum -= (*_val_i_vars[i-cur_size-1])[_qp]*(*_val_v_vars[i-1])[_qp]*_gc._absorb(i,cur_size-i);
      }
    }

    //v emission loss(-)
    if(cur_size!=1){
      res_sum += _u[_qp]*_gc._emit(cur_size);
      //printf("emission %d (-): %d\n",cur_size,cur_size);     
    }
    
    //v+1 emission gain(+)
    if(cur_size<_number_v){
      res_sum -= (*_val_v_vars[cur_size])[_qp]*_gc._emit(cur_size+1);
      //printf("emission gain %d (+): %d\n",cur_size,cur_size+1);     
    }
    if(cur_size==1){
      for(int i=2;i<=_number_v;i++){
        res_sum -= (*_val_v_vars[i-1])[_qp]*_gc._emit(i);
        //printf("emission gain %d (+): %d\n",cur_size,i);     
      }
    }

    //dislocation loss(-)
    res_sum += _u[_qp]*_gc._disl(cur_size);
  }

  else{//i type
   
    cur_size = -_cur_size;//make it positive
    max_vi = std::min(cur_size+vv,_number_i);

    //iv reaction loss(-)
    for(int i=1;i<=_number_v;i++){
      res_sum += (*_val_v_vars[i-1])[_qp]*_u[_qp]*_gc._absorb(-cur_size,i);
      //printf("reaction %d (-): %d %d\n",_cur_size,_cur_size,i);     
    }

    //ii reaction loss(-)
    for(int i=1;i <= _number_i-cur_size;i++){//garantee the largest size doesn't exceed _number_v
      //printf("reaction %d (-): %d %d\n",_cur_size,_cur_size,-i);     
      res_sum += (*_val_i_vars[i-1])[_qp]*_u[_qp]*_gc._absorb(-cur_size,-i);
    }
    if(cur_size*2<=_number_i){
      res_sum += (*_val_i_vars[cur_size-1])[_qp]*_u[_qp]*_gc._absorb(-cur_size,-cur_size);
      //printf("reaction %d (-): %d %d\n",_cur_size,_cur_size,_cur_size);     
    } 

    //ii reaction gain(+)
    for(int i=1;i <= (int)(cur_size/2);i++){
      //if(std::find(_i_size.begin(),_i_size.end(),cur_size-i) != _i_size.end() || std::find(_i_size.begin(),_i_size.end(),i) != _i_size.end()){//TODO: can optimize
        res_sum -= (*_val_i_vars[i-1])[_qp]*(*_val_i_vars[cur_size-i-1])[_qp]*_gc._absorb(i-cur_size,-i);
        //printf("reaction %d (+): %d %d\n",_cur_size,i-cur_size,-i);     
      //}
    }
  
    //iv reaction gain(+)
    for(int i=cur_size+1;i<=max_vi;i++){
      //if(std::find(_v_size.begin(),_v_size.end(),i-cur_size) != _v_size.end() || std::find(_i_size.begin(),_i_size.end(),i) != _i_size.end()){
      if(i-cur_size <= vv || i <= ii){//make sure one is mobile
        res_sum -= (*_val_v_vars[i-cur_size-1])[_qp]*(*_val_i_vars[i-1])[_qp]*_gc._absorb(-i,i-cur_size);
        //printf("reaction %d (+): %d %d\n",_cur_size,i-cur_size,-i);     
      }
    }

    //i emission loss(-)
    if(cur_size!=1){
      res_sum += _u[_qp]*_gc._emit(-cur_size);
      //printf("emission %d (-): %d\n",_cur_size,_cur_size);     
    }
    
    //i+1 emission gain(+)
    if(cur_size<_number_i){
      res_sum -= (*_val_i_vars[cur_size])[_qp]*_gc._emit(-cur_size-1);
      //printf("emission %d (+): %d\n",_cur_size,-cur_size-1);     
    }
    if(cur_size==1){
      for(int i=2;i<=_number_i;i++){
        res_sum -= (*_val_i_vars[i-1])[_qp]*_gc._emit(-i);
       // printf("emission %d (+): %d\n",_cur_size,-i);     
      }
    }

    //dislocation loss(-)
    res_sum += _u[_qp]*_gc._disl(-cur_size);
  }
  return res_sum*_test[_i][_qp];
}

Real
MobileDefects::computeQpJacobian()
{
  Real jac_sum = 0.0;
  int cur_size;//should be positive value
  if(_cur_size>0){//v type
    
    cur_size = _cur_size; 

    //vi reaction loss(-)
    for(int i=1;i<=_number_i;i++){
      jac_sum += (*_val_i_vars[i-1])[_qp]*_gc._absorb(cur_size,-i);
    }

    //vv reaction loss(-)
    for(int i=1;i <= _number_v-cur_size;i++){//garantee the largest size doesn't exceed _number_v
      jac_sum += (*_val_v_vars[i-1])[_qp]*_gc._absorb(cur_size,i);//10 is test
    }
    if(cur_size*2<=_number_v)//2*u^2->4*u*phi
      jac_sum += 3.0*_u[_qp]*_gc._absorb(cur_size,cur_size);
//(*_val_v_vars[cur_size-1])[_qp]
  
    //v emission loss(-)
    jac_sum += _gc._emit(cur_size);
    
    //dislocation loss(-)
    jac_sum += _gc._disl(cur_size);
  }

  else{//i type
   
    cur_size = -_cur_size;//make it positive

    //iv reaction loss(-)
    for(int i=1;i<=_number_v;i++){
      jac_sum += (*_val_v_vars[i-1])[_qp]*_gc._absorb(-cur_size,i);
    }

    //ii reaction loss(-)
    for(int i=1;i <= _number_i-cur_size;i++){//garantee the largest size doesn't exceed _number_v
      jac_sum += (*_val_i_vars[i-1])[_qp]*_gc._absorb(-cur_size,-i);
    }
    if(cur_size*2<=_number_i)
      jac_sum += 3.0*_u[_qp]*_gc._absorb(-cur_size,-cur_size);// *(*_val_i_vars[cur_size-1])[_qp]
  
    //i emission loss(-)
    jac_sum += _gc._emit(-cur_size);
    
    //dislocation loss(-)
    jac_sum += _gc._disl(-cur_size);
  }
  return jac_sum*_test[_i][_qp] * _phi[_j][_qp];
}

Real 
MobileDefects::computeQpOffDiagJacobian(unsigned int jvar){
      return 0.0;//not provided now
}


int
MobileDefects::getGroupNumber(std::string str)
{
  int len=str.length(),i=len;
  while(std::isdigit(str[i-1])) i--;
  int no = std::atoi((str.substr(i)).c_str());
  while(i>=0){
      i--;
      if(str[i]=='v'){no = no;break;}//v type: "+"
      if(str[i]=='i'){no = -no;break;}//- type: "-"

  }
  return no;
}
