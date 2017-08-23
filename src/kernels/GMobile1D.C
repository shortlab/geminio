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

#include "GMobile1D.h"
#include "Conversion.h"
#define DEBUG 0

template<>
InputParameters validParams<GMobile1D>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<int>("number_v","Maximum vacancy cluster size");
  params.addRequiredParam<int>("number_i","Maximum interstitial cluster size");
  params.addCoupledVar("coupled_i_auxvars","coupled reciprocle mean free path ");
  params.addCoupledVar("coupled_v_vars","coupled vacancy type variables");
  params.addCoupledVar("coupled_i_vars","coupled intersitial type variables");
  params.addRequiredParam<int>("max_mobile_v", "A vector of mobile species");
  params.addRequiredParam<int>("max_mobile_i", "A vector of mobile species");
  params.addRequiredParam<UserObjectName>("user_object","The name of user object providing interaction constants");
  return params;
}

GMobile1D::GMobile1D(const InputParameters & parameters)
     :Kernel(parameters),
     _number_v(getParam<int>("number_v")),
     _number_i(getParam<int>("number_i")),
     _max_mobile_v(getParam<int>("max_mobile_v")),
     _max_mobile_i(getParam<int>("max_mobile_i")),
     _gc(getUserObject<GGroup>("user_object"))
{
  int nvcoupled = coupledComponents("coupled_v_vars");
  int nicoupled = coupledComponents("coupled_i_vars");
  int nicoupledaux = coupledComponents("coupled_i_auxvars");
  if(_number_v>0){
    _no_v_vars.resize(nvcoupled);
    _val_v_vars.resize(nvcoupled);
  }
  if(_number_i>0){
    _no_i_vars.resize(nicoupled);
    _val_i_vars.resize(nicoupled);
  }
  _val_i_auxvars.resize(nicoupledaux);

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

  for (int i=0;i<coupledComponents("coupled_i_auxvars");i++){
    _val_i_auxvars[i] = &coupledValue("coupled_i_auxvars",i);
  }
 
  if(DEBUG){
    std::vector<VariableName> coupled_v_vars = getParam<std::vector<VariableName> >("coupled_v_vars");
    std::cout << "GMobile1D: current variable => " << cur_var_name << std::endl;
    std::cout << "coupled with: " << std::endl;
    for (int i=0; i < nvcoupled; ++i){
      std::cout << coupled_v_vars[i] << "  ";  
    }
    std::cout << std::endl;
  } 
}

/* use approximation: kinetic parameter is constant for each size within one group */
Real
GMobile1D::computeQpResidual()
{
  Real res_sum = 0.0, partial_res = 0.0;
  int cur_size;//should be positive value
  int ii = _max_mobile_i;
  int vv = _max_mobile_v;
  int max_vi;
  Real conc,conci,concj;
  Real group_conc_sum;
  Real absorb_coef, absorb_coef_mod;

  if(_cur_size>0){//v type
    cur_size = _cur_size; 
    max_vi =  std::min(_cur_size+ii,max_v);


    //vi reaction loss(-)
    for(int i=1;i<=_max_mobile_i;i++){
      group_conc_sum = (*_val_i_vars[2*(i-1)])[_qp]*_gc.GroupScheme_i_del[i-1];
      absorb_coef = getAbsorbCoef(i); 
      partial_res += group_conc_sum *_gc._absorb(cur_size,-i)*absorb_coef;//group cur_size and group i; sigma_i,j*Di*2/lambda_i (assume interstitial mobility dominiate)
    }
    for(int i=_max_mobile_i+1;i<=_number_i;i++){
      group_conc_sum = (*_val_i_vars[2*(i-1)])[_qp]*_gc.GroupScheme_i_del[i-1];
      partial_res += group_conc_sum *_gc._absorb(cur_size,-i);//group cur_size and group i; 
    }
    res_sum += partial_res * _u[_qp];

    //vv reaction loss(-)
    partial_res = 0.0;
    int limit_group = _gc.CurrentGroupV(max_v-cur_size);
    for(int i=1;i <= limit_group;i++){//garantee the largest size doesn't exceed GroupScheme_v.back()
      group_conc_sum = (*_val_v_vars[2*(i-1)])[_qp]*_gc.GroupScheme_v_del[i-1];
      partial_res += group_conc_sum *_gc._absorb(cur_size,i);//group cur_size and group i; 
    }
    for(int i=max_v-cur_size+1;i<=_gc.GroupScheme_v[limit_group];i++){//deduct over added terms
      partial_res -= _gc._absorb(cur_size,limit_group)*((*_val_v_vars[2*(limit_group-1)])[_qp]+(*_val_v_vars[2*(limit_group-1)+1])[_qp]*(i-_gc.GroupScheme_v_avg[limit_group-1]));
    }
    res_sum += partial_res * _u[_qp];
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
        int ig = _gc.CurrentGroupI(i-cur_size);
        int vg = _gc.CurrentGroupV(i);
        absorb_coef = (ig>_max_mobile_i)? 1.0:getAbsorbCoef(ig);
        res_sum -= conci * concj * _gc._absorb(vg,-ig)*absorb_coef;
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
      res_sum -= conc *_gc._emit(_gc.CurrentGroupV(cur_size+1));
    }
    if(cur_size==1){
      for(int i=2;i<=_number_v;i++){
        //printf("emission gain %d (+): %d\n",cur_size,i);     
        group_conc_sum = (*_val_v_vars[2*(i-1)])[_qp]*_gc.GroupScheme_v_del[i-1];
        res_sum -= group_conc_sum * _gc._emit(i);
        //printf("res: %f %d %f\n",res_sum,cur_size,_gc._emit(i));
      }
    }

    //dislocation loss(-)
    res_sum += _u[_qp]*_gc._disl(cur_size);
  }

  else{//i type 
   
    cur_size = -_cur_size;//make it positive
    max_vi = std::min(cur_size+vv,max_i);
    absorb_coef = getAbsorbCoef(cur_size);

    //iv reaction loss(-)
    for(int i=1;i<=_number_v;i++){
      group_conc_sum = (*_val_v_vars[2*(i-1)])[_qp] * _gc.GroupScheme_v_del[i-1];
      partial_res += group_conc_sum * _gc._absorb(-cur_size,i)*absorb_coef;
    }
    res_sum += partial_res *_u[_qp];

    //ii reaction loss(-)
    partial_res = 0.0;
    int limit_group = _gc.CurrentGroupI(max_i-cur_size);
    for(int i=1;i <= limit_group;i++){//garantee the largest size doesn't exceed GroupScheme_i.back()
      group_conc_sum = (*_val_i_vars[2*(i-1)])[_qp] * _gc.GroupScheme_i_del[i-1];
      absorb_coef_mod = absorb_coef; 
      if (i<=_max_mobile_i)
        absorb_coef_mod += getAbsorbCoef(i);//account for both mobile
      partial_res += group_conc_sum * _gc._absorb(-cur_size,-i)*absorb_coef_mod;
    }

    absorb_coef_mod = absorb_coef; 
    if (limit_group<=_max_mobile_i)
      absorb_coef_mod += getAbsorbCoef(limit_group);//account for both mobile
    for(int i=max_i-cur_size+1;i<=_gc.GroupScheme_i[limit_group];i++){//deduct over added terms
      partial_res -= _gc._absorb(-cur_size,-limit_group)*absorb_coef_mod*((*_val_i_vars[2*(limit_group-1)])[_qp]+(*_val_i_vars[2*(limit_group-1)+1])[_qp]*(i-_gc.GroupScheme_v_avg[limit_group-1]));
    }
    res_sum += partial_res * _u[_qp];
    if(cur_size*2<=max_i){
      res_sum += _u[_qp]*_u[_qp]*_gc._absorb(-cur_size,-cur_size)*absorb_coef*2;
      //printf("reaction %d (-): %d %d\n",_cur_size,_cur_size,_cur_size);     
    } 

    //ii reaction gain(+)
    for(int i=1;i <= (int)(cur_size/2);i++){
        conci = getConcBySize(-i);
        concj = getConcBySize(i-cur_size);
        res_sum -= conci * concj *_gc._absorb(i-cur_size,-i)*(getAbsorbCoef(cur_size-i)+getAbsorbCoef(i));
        //printf("reaction %d (+): %d %d\n",_cur_size,i-cur_size,-i);     
      //}
    }
  
    //iv reaction gain(+)
    for(int i=cur_size+1;i<=max_vi;i++){
      if(i-cur_size <= vv || i <= ii){//make sure one is mobile
        conci = getConcBySize(-i);
        concj = getConcBySize(i-cur_size);
        int ig = _gc.CurrentGroupI(i);
        int vg = _gc.CurrentGroupV(i-cur_size);
        absorb_coef_mod = (ig>_max_mobile_i)? 1.0:getAbsorbCoef(ig);
        res_sum -= conci * concj *_gc._absorb(-ig,vg)*absorb_coef_mod;
        //printf("reaction %d (+): %d %d\n",_cur_size,i-cur_size,-i);     
      }
    }

    //dislocation loss(-)
    res_sum += _u[_qp]*_gc._disl(-cur_size);
  }
  return res_sum*_test[_i][_qp];
}

Real
GMobile1D::computeQpJacobian()
{
  Real jac_sum = 0.0;
  int cur_size;//should be positive value
  Real conc,group_conc_sum;
  Real absorb_coef,absorb_coef_mod;
  if(_cur_size>0){//v type
    
    cur_size = _cur_size; 

    //vi reaction loss(-)
    for(int i=1;i<=_number_i;i++){
      group_conc_sum = (*_val_i_vars[2*(i-1)])[_qp]*_gc.GroupScheme_i_del[i-1];
      absorb_coef = (i<=_max_mobile_i)?getAbsorbCoef(i):1.0;
      jac_sum += group_conc_sum *_gc._absorb(cur_size,-i)*absorb_coef;//group cur_size and group i; 
    }

    //vv reaction loss(-)
    int limit_group = _gc.CurrentGroupV(max_v-cur_size);
    for(int i=1;i <= limit_group;i++){//garantee the largest size doesn't exceed GroupScheme_v.back()
      group_conc_sum = (*_val_v_vars[2*(i-1)])[_qp]*_gc.GroupScheme_v_del[i-1];
      jac_sum += group_conc_sum *_gc._absorb(cur_size,i);//group cur_size and group i; 
    }
    for(int i=max_v-cur_size+1;i<=_gc.GroupScheme_v[limit_group];i++){//deduct over added terms
      jac_sum -= _gc._absorb(cur_size,limit_group)*((*_val_v_vars[2*(limit_group-1)])[_qp]+(*_val_v_vars[2*(limit_group-1)+1])[_qp]*(i-_gc.GroupScheme_v_avg[limit_group-1]));
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
    absorb_coef = getAbsorbCoef(cur_size);

    //iv reaction loss(-)
    for(int i=1;i<=_number_v;i++){
      group_conc_sum = (*_val_v_vars[2*(i-1)])[_qp] * _gc.GroupScheme_v_del[i-1];
      jac_sum += group_conc_sum * _gc._absorb(-cur_size,i)*absorb_coef;
    }

    //ii reaction loss(-)
    int limit_group = _gc.CurrentGroupI(max_i-cur_size);
    for(int i=1;i <= limit_group;i++){//garantee the largest size doesn't exceed GroupScheme_i.back()
      group_conc_sum = (*_val_i_vars[2*(i-1)])[_qp] * _gc.GroupScheme_i_del[i-1];
      absorb_coef_mod = absorb_coef; 
      if (i<=_max_mobile_i)
        absorb_coef_mod += getAbsorbCoef(i);//account for both mobile
      jac_sum += group_conc_sum * _gc._absorb(-cur_size,-i)*absorb_coef_mod;
    }

    absorb_coef_mod = absorb_coef; 
    if (limit_group<=_max_mobile_i)
      absorb_coef_mod += getAbsorbCoef(limit_group);//account for both mobile
    for(int i=max_i-cur_size+1;i<=_gc.GroupScheme_i[limit_group];i++){//deduct over added terms
      jac_sum -= _gc._absorb(-cur_size,-limit_group)*absorb_coef_mod*((*_val_i_vars[2*(limit_group-1)])[_qp]+(*_val_i_vars[2*(limit_group-1)+1])[_qp]*(i-_gc.GroupScheme_v_avg[limit_group-1]));
    }
    if(cur_size*2<=max_i)
      jac_sum += 3.0*_u[_qp]*_gc._absorb(-cur_size,-cur_size)*absorb_coef*2;// *(*_val_i_vars[cur_size-1])[_qp]
  
    //dislocation loss(-)
    jac_sum += _gc._disl(-cur_size);
  }
  return jac_sum*_test[_i][_qp] * _phi[_j][_qp];
}



Real 
GMobile1D::computeQpOffDiagJacobian(unsigned int jvar){
      return 0.0;//not provided now
}


int
GMobile1D::getGroupNumber(std::string str)
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
GMobile1D::getConcBySize(int i)
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

double
GMobile1D::getAbsorbCoef(int i)//for 1D diffusers(valid for i<=_max_mobile_i
{
  if (i>_max_mobile_i)
    mooseError("Call reciprocal mean free path out of range");
  Real absorb_coef = 2*_gc._diff(-i)*(*_val_i_auxvars[i-1])[_qp];
  return absorb_coef;
}
