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

#include "GMobile.h"
#include "Conversion.h"
#define DEBUG 0

template<>
InputParameters validParams<GMobile>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<int>("number_v","Maximum vacancy cluster size");
  params.addRequiredParam<int>("number_i","Maximum interstitial cluster size");
  params.addCoupledVar("coupled_v_vars","coupled vacancy type variables");
  params.addCoupledVar("coupled_i_vars","coupled intersitial type variables");
  params.addRequiredParam<int>("max_mobile_v", "A vector of mobile species");
  params.addRequiredParam<int>("max_mobile_i", "A vector of mobile species");
  params.addRequiredParam<UserObjectName>("user_object","The name of user object providing interaction constants");
  return params;
}

GMobile::GMobile(const InputParameters & parameters)
     :Kernel(parameters),
     _number_v(getParam<int>("number_v")),
     _number_i(getParam<int>("number_i")),
     _max_mobile_v(getParam<int>("max_mobile_v")),
     _max_mobile_i(getParam<int>("max_mobile_i")),
     _gc(getUserObject<GGroup>("user_object"))
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
    _no_v_vars[i] = coupled("coupled_v_vars",i);
    _val_v_vars[i] = &coupledValue("coupled_v_vars",i);
  }
  for (int i=0; i < nicoupled; ++i)
  {
    _no_i_vars[i] = coupled("coupled_i_vars",i);
    _val_i_vars[i] = &coupledValue("coupled_i_vars",i);
  }
  NonlinearVariableName cur_var_name = getParam<NonlinearVariableName>("variable");
  _cur_size = getGroupNumber(cur_var_name.c_str());

  max_i = (_gc.GroupScheme_i.size()>0?_gc.GroupScheme_i.back():0);
  max_v = (_gc.GroupScheme_v.size()>0?_gc.GroupScheme_v.back():0);
 
  if(DEBUG){
    std::vector<VariableName> coupled_v_vars = getParam<std::vector<VariableName> >("coupled_v_vars");
    std::cout << "GMobile: current variable => " << cur_var_name << std::endl;
    std::cout << "coupled with: " << std::endl;
    for (int i=0; i < nvcoupled; ++i){
      std::cout << coupled_v_vars[i] << "  ";  
    }
    std::cout << std::endl;
  } 
}

Real
GMobile::computeQpResidual()
{
  Real res_sum = 0.0;
  int cur_size;//should be positive value
  int ii = _max_mobile_i;
  int vv = _max_mobile_v;
  int max_vi;
  Real conc,conci,concj;

  //printf("return mobile initial: %f %d\n",res_sum,_cur_size);     

  if(_cur_size>0){//v type
    cur_size = _cur_size; 
    max_vi =  std::min(_cur_size+ii,max_v);


    //vi reaction loss(-)
    for(int i=1;i<=max_i;i++){
      conc = getConcBySize(-i);
      res_sum += conc * _u[_qp] * _gc._absorb(cur_size,-i);
      //printf("vi reaction %d (-): %d %d\n",cur_size,cur_size,-i);     
    }

    //vv reaction loss(-)
    for(int i=1;i <= max_v-cur_size;i++){//garantee the largest size doesn't exceed _number_v
      //printf("vv reaction %d (-): %d %d\n",cur_size,cur_size,i);     
      conc = getConcBySize(i);
      res_sum += conc*_u[_qp]*_gc._absorb(cur_size,i);
    }
    if(cur_size*2 <= max_v){
      //printf("vv reaction %d (-): %d %d\n",cur_size,cur_size,cur_size);     
      res_sum += _u[_qp]*_u[_qp]*_gc._absorb(cur_size,cur_size);
    }
      
  
    //vv reaction gain(+)
    for(int i=1;i <= (int)(cur_size/2);i++){
        //printf("vv reaction %d (+): %d %d\n",cur_size,cur_size-i,cur_size);     
        conci = getConcBySize(cur_size-i);
        concj = getConcBySize(i);
        res_sum -= conci * concj *_gc._absorb(cur_size-i,i);
    }

    //vi reaction gain(+)
    for(int i=cur_size+1;i<=max_vi;i++){
      if(i-cur_size <= ii || i <= vv ){//make sure one is mobile
        conci = getConcBySize(cur_size-i);
        concj = getConcBySize(i);
        res_sum -= conci * concj * _gc._absorb(i,cur_size-i);
        //printf("vi reaction %d (+): %d %d\n",cur_size,cur_size-i,i);     
      }
    }

  

    //v emission loss(-)
    if(cur_size!=1){
      res_sum += _u[_qp]*_gc._emit(cur_size);
      //printf("emission loss %d (-): %d\n",cur_size,cur_size);     
    }
    
    //v+1 emission gain(+)
    if(cur_size<max_v){
      //printf("emission gain %d (+): %d\n",cur_size,cur_size+1);     
      conc = getConcBySize(cur_size+1);
      res_sum -= conc *_gc._emit(cur_size+1);
    }
    if(cur_size==1){
      for(int i=2;i<=max_v;i++){
        //printf("emission gain %d (+): %d\n",cur_size,i);     
        conc = getConcBySize(i);
        res_sum -= conc * _gc._emit(i);
        //printf("res: %f %d %f\n",res_sum,cur_size,_gc._emit(i));
      }
    }

    //dislocation loss(-)
    res_sum += _u[_qp]*_gc._disl(cur_size);

  }

  else{//i type
   
    cur_size = -_cur_size;//make it positive
    max_vi = std::min(cur_size+vv,max_i);

    //iv reaction loss(-)
    for(int i=1;i<=max_v;i++){
      conc = getConcBySize(i);
      res_sum += conc *_u[_qp]*_gc._absorb(-cur_size,i);
      //printf("reaction %d (-): %d %d\n",_cur_size,_cur_size,i);     
    }

    //ii reaction loss(-)
    for(int i=1;i <= max_i-cur_size;i++){//garantee the largest size doesn't exceed _number_v
      //printf("reaction %d (-): %d %d\n",_cur_size,_cur_size,-i);     
      conc = getConcBySize(-i);
      res_sum += conc *_u[_qp]*_gc._absorb(-cur_size,-i);
    }
    if(cur_size*2<=max_i){
      res_sum += _u[_qp]*_u[_qp]*_gc._absorb(-cur_size,-cur_size);
      //printf("reaction %d (-): %d %d\n",_cur_size,_cur_size,_cur_size);     
    } 

    //ii reaction gain(+)
    for(int i=1;i <= (int)(cur_size/2);i++){
        conci = getConcBySize(-i);
        concj = getConcBySize(i-cur_size);
        res_sum -= conci * concj *_gc._absorb(i-cur_size,-i);
        //printf("reaction %d (+): %d %d\n",_cur_size,i-cur_size,-i);     
      //}
    }
  
    //iv reaction gain(+)
    for(int i=cur_size+1;i<=max_vi;i++){
      if(i-cur_size <= vv || i <= ii){//make sure one is mobile
        conci = getConcBySize(-i);
        concj = getConcBySize(i-cur_size);
        res_sum -= conci * concj *_gc._absorb(-i,i-cur_size);
        //printf("reaction %d (+): %d %d\n",_cur_size,i-cur_size,-i);     
      }
    }

    //i emission loss(-)
    if(cur_size!=1){
      res_sum += _u[_qp]*_gc._emit(-cur_size);
      //printf("emission %d (-): %d\n",_cur_size,_cur_size);     
    }
    
    //i+1 emission gain(+)
    if(cur_size<max_i){
      conc = getConcBySize(-cur_size-1);
      res_sum -= conc *_gc._emit(-cur_size-1);
      //printf("emission %d (+): %d\n",_cur_size,-cur_size-1);     
    }
    if(cur_size==1){
      for(int i=2;i<=max_i;i++){
        conc = getConcBySize(-i);
        res_sum -= conc *_gc._emit(-i);
       // printf("emission %d (+): %d\n",_cur_size,-i);     
      }
    }

    //dislocation loss(-)
    res_sum += _u[_qp]*_gc._disl(-cur_size);
  }
  return res_sum*_test[_i][_qp];
}

Real
GMobile::computeQpJacobian()
{
  Real jac_sum = 0.0;
  int cur_size;//should be positive value
  Real conc;
  if(_cur_size>0){//v type
    
    cur_size = _cur_size; 

    //vi reaction loss(-)
    for(int i=1;i<=max_i;i++){
      conc = getConcBySize(-i);
      jac_sum += conc *_gc._absorb(cur_size,-i);
    }

    //vv reaction loss(-)
    for(int i=1;i <= max_v-cur_size;i++){//garantee the largest size doesn't exceed _number_v
      conc = getConcBySize(i);
      jac_sum += conc *_gc._absorb(cur_size,i);//10 is test
    }
    if(cur_size*2<=max_v)//2*u^2->4*u*phi
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
    for(int i=1;i<=max_v;i++){
      conc = getConcBySize(i);
      jac_sum += conc *_gc._absorb(-cur_size,i);
    }

    //ii reaction loss(-)
    for(int i=1;i <= max_i-cur_size;i++){//garantee the largest size doesn't exceed _number_v
      conc = getConcBySize(-i);
      jac_sum += conc*_gc._absorb(-cur_size,-i);
    }
    if(cur_size*2<=max_i)
      jac_sum += 3.0*_u[_qp]*_gc._absorb(-cur_size,-cur_size);// *(*_val_i_vars[cur_size-1])[_qp]
  
    //i emission loss(-)
    jac_sum += _gc._emit(-cur_size);
    
    //dislocation loss(-)
    jac_sum += _gc._disl(-cur_size);
  }
  return jac_sum*_test[_i][_qp] * _phi[_j][_qp];
}

Real 
GMobile::computeQpOffDiagJacobian(unsigned int jvar){
      return 0.0;//not provided now
}


int
GMobile::getGroupNumber(std::string str)
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

double
GMobile::getConcBySize(int i)
{
  if(i>0){
    int g = _gc.CurrentGroupV(i);
    return (*_val_v_vars[2*(g-1)])[_qp]+(*_val_v_vars[2*(g-1)+1])[_qp]*(i-_gc.GroupScheme_v_avg[g-1]);
  }
 else{
    int g = _gc.CurrentGroupI(-i);
    return (*_val_i_vars[2*(g-1)])[_qp]+(*_val_i_vars[2*(g-1)+1])[_qp]*(-i-_gc.GroupScheme_i_avg[g-1]);
 } 
}
