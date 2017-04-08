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

#include "GImmobileL0.h"
#include "Conversion.h"
#define DEBUG 0

template<>
InputParameters validParams<GImmobileL0>()
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

GImmobileL0::GImmobileL0(const InputParameters & parameters)
     :Kernel(parameters),
     _number_v(getParam<int>("number_v")),
     _number_i(getParam<int>("number_i")),
     _max_mobile_v(getParam<int>("max_mobile_v")),
     _max_mobile_i(getParam<int>("max_mobile_i")),
     _gc(getUserObject<GGroup>("user_object"))
{
  NonlinearVariableName cur_var_name = getParam<NonlinearVariableName>("variable");
  _cur_size = getGroupNumber(cur_var_name);
  if(_cur_size>0 && _cur_size< _max_mobile_v){
    mooseError("Please check the mobile vacancy cluster size range, should be immobile>mobile");
  }
  else if(_cur_size<0 && -_cur_size<_max_mobile_i){
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
    
  if(DEBUG){
    std::vector<VariableName> coupled_v_vars = getParam<std::vector<VariableName> >("coupled_v_vars");
    std::cout << "GImmobileL0: current variable => " << cur_var_name << std::endl;
    std::cout << "coupled with: " << std::endl;
    for (int i=0; i < num_v_coupled; ++i){
      std::cout << i << ":" << coupled_v_vars[i] << "  ";  
    }
    std::cout << std::endl;
  } 
}

Real
GImmobileL0::computeQpResidual()
{
  int cur_size,other_size,group_num;
  Real res_sum = 0.0;
  Real conc1,conc2,conc;
  //printf("return immobile initial: %f %d\n",res_sum,_cur_size);     

  if(_cur_size>0){//v type
    //printf("L0: Start; current V group: %d\n",_cur_size);
    int index = 2*_max_mobile_v;//current variable index (start from 0)
    cur_size = _cur_size;

    //left boundary x_{i-1}+1, absorb the same species
    for(int i=0;i<=_max_mobile_v-1;i++){
      int tmp_size = std::min(_max_mobile_v-1-i,_gc.GroupScheme_v_del[cur_size-1]-1);
      tmp_size = std::min(tmp_size, _gc.GroupScheme_v[cur_size-1]-1-2*i);//prevent duplicating pair from mobile ones.
      group_num =  _gc.CurrentGroupV(_gc.GroupScheme_v[cur_size-1]-i);
      other_size = _max_mobile_v+_max_mobile_v-(cur_size-group_num);
      conc1 = (*_val_v_vars[2*other_size])[_qp]+(*_val_v_vars[2*other_size+1])[_qp]*(_gc.GroupScheme_v[cur_size-1]-i-_gc.GroupScheme_v_avg[group_num-1]);
      for(int j=0;j<=tmp_size;j++){
        conc2 = (*_val_v_vars[2*(i+j)])[_qp];
        res_sum -= conc1 * conc2 * _gc._absorb(_gc.GroupScheme_v[cur_size-1]-i,i+j+1);
        //printf("absorb %d (vv gain): %d and %d; var: %d %d\n",_cur_size,_gc.GroupScheme_v[cur_size-1]-i,i+j+1,2*(i+j),2*other_size);
      }//vv (gain)
    }

    //left boundary x_{i-1}+1, emission
    conc = _u[_qp]+(_gc.GroupScheme_v[cur_size-1]+1-_gc.GroupScheme_v_avg[cur_size-1])*(*_val_v_vars[2*index+1])[_qp];
    res_sum += conc * _gc._emit(_gc.GroupScheme_v[cur_size-1]+1);//v emit (loss)
    //printf("emit %d (v loss): %d; var: %d\n",_cur_size,_gc.GroupScheme_v[cur_size-1]+1,2*index);

    //left boundary x_{i-1}+1, absorb the opposite species
    for(int i=0;i<=_max_mobile_i-1;i++){
      int tmp_size = std::min(_max_mobile_i-1-i,_gc.GroupScheme_v_del[cur_size-1]-1);
      for(int j=0;j<=tmp_size;j++){
        conc1 = _u[_qp] + (_gc.GroupScheme_v[cur_size-1]+j+1-_gc.GroupScheme_v_avg[cur_size-1])*(*_val_v_vars[2*index+1])[_qp];
        conc2 = (*_val_i_vars[2*(i+j)])[_qp];
        res_sum += conc1 * conc2 * _gc._absorb(_gc.GroupScheme_v[cur_size-1]+j+1,-(i+j+1));
        //printf("absorb %d (vi loss): %d and %d\n",_cur_size,_gc.GroupScheme_v[cur_size-1]+j+1,-(i+j+1));
      }//vi (loss)
    }
    

    if(cur_size != (int)(_gc.GroupScheme_v.size()-1)){

      //right boundary x_{i}+1, absorb the same species
      int tmp_size = std::min(_max_mobile_v-1,_gc.GroupScheme_v_del[cur_size-1]-1);
      for(int i=0;i<=tmp_size;i++){
        conc1 = _u[_qp] + (_gc.GroupScheme_v[cur_size]-i-_gc.GroupScheme_v_avg[cur_size-1])*(*_val_v_vars[2*index+1])[_qp];
        for(int j=0;j<=_max_mobile_v-1-i;j++){
          conc2 = (*_val_v_vars[2*(i+j)])[_qp];
          res_sum += conc1 * conc2 * _gc._absorb(_gc.GroupScheme_v[cur_size]-i,i+j+1);
          //printf("absorb %d (vv loss): %d and %d; var: %d %d\n",_cur_size,_gc.GroupScheme_v[cur_size]-i,i+j+1,2*(i+j),2*index);
        }//vv (loss)
      }
  
      //right boundary x_{i}+1, absorb the opposite species
      tmp_size = std::min(_max_mobile_i-1,_gc.GroupScheme_v_del[cur_size-1]-1);
      for(int i=0;i<=tmp_size;i++){
        int tmp2 = std::min(_max_mobile_i-1-i,_gc.GroupScheme_v.back()-_gc.GroupScheme_v[cur_size]-1);
        for(int j=0;j<=tmp2;j++){
          group_num =  _gc.CurrentGroupV(_gc.GroupScheme_v[cur_size]+j+1);
          //other_size = index+i+j+1;
          other_size = index+(group_num-cur_size);
          conc1 = (*_val_v_vars[2*other_size])[_qp]+(*_val_v_vars[2*other_size+1])[_qp]*(_gc.GroupScheme_v[cur_size]+j+1-_gc.GroupScheme_v_avg[group_num-1]);
          conc2 = (*_val_i_vars[2*(i+j)])[_qp]; 
          res_sum -= conc1 * conc2 * _gc._absorb(_gc.GroupScheme_v[cur_size]+j+1,-(i+j+1));
          //printf("absorb %d (vi gain): %d and %d\n",_cur_size,_gc.GroupScheme_v[cur_size]+j+1,-(i+j+1));
        }//vi (gain)
      }

      //right boundary x_{i}+1, emission
      conc = (*_val_v_vars[2*index+2])[_qp]+(_gc.GroupScheme_v[cur_size]+1-_gc.GroupScheme_v_avg[cur_size])*(*_val_v_vars[2*index+3])[_qp];
      res_sum -= conc * _gc._emit(_gc.GroupScheme_v[cur_size]+1);//v emit (gain)
      //printf("emit %d (v gain): %d; var: %d\n",_cur_size,_gc.GroupScheme_v[cur_size]+1,2*(index+1));
    }  

    //printf("L0 END\n");
    return 1.0/(_gc.GroupScheme_v_del[cur_size-1])*res_sum *_test[_i][_qp];
  
  }

  else{
    int index = 2*_max_mobile_i;//-1;//current variable index (start from 0)
    cur_size = -_cur_size;

    //left boundary x_{i-1}+1, absorb the same species
    for(int i=0;i<=_max_mobile_i-1;i++){
      int tmp_size = std::min(_max_mobile_i-1-i,_gc.GroupScheme_i_del[cur_size-1]-1);
      tmp_size = std::min(tmp_size, _gc.GroupScheme_i[cur_size-1]-1-2*i);//prevent duplicating pair from mobile ones.
      group_num =  _gc.CurrentGroupI(_gc.GroupScheme_i[cur_size-1]-i);
      other_size = _max_mobile_i+_max_mobile_i-(cur_size-group_num);
      conc1 = (*_val_i_vars[2*other_size])[_qp]+(*_val_i_vars[2*other_size+1])[_qp]*(_gc.GroupScheme_i[cur_size-1]-i-_gc.GroupScheme_i_avg[group_num-1]);
      for(int j=0;j<=tmp_size;j++){
        conc2 = (*_val_i_vars[2*(i+j)])[_qp];
        res_sum -= conc1 * conc2 * _gc._absorb(-(_gc.GroupScheme_i[cur_size-1]-i),-(i+j+1));
      }//ii (gain)
    }

    //left boundary x_{i-1}+1, emission
    conc = _u[_qp]+(_gc.GroupScheme_i[cur_size-1]+1-_gc.GroupScheme_i_avg[cur_size-1])*(*_val_i_vars[2*index+1])[_qp];
    res_sum += conc * _gc._emit(-(_gc.GroupScheme_i[cur_size-1]+1));//i emit (loss)

    //left boundary x_{i-1}+1, absorb the opposite species
    for(int i=0;i<=_max_mobile_v-1;i++){
      int tmp_size = std::min(_max_mobile_v-1-i,_gc.GroupScheme_i_del[cur_size-1]-1);
      for(int j=0;j<=tmp_size;j++){
        conc1 = _u[_qp] + (_gc.GroupScheme_i[cur_size-1]+j+1-_gc.GroupScheme_i_avg[cur_size-1])*(*_val_i_vars[2*index+1])[_qp];
        conc2 = (*_val_v_vars[2*(i+j)])[_qp];
        res_sum += conc1 * conc2 * _gc._absorb(-(_gc.GroupScheme_i[cur_size-1]+j+1),(i+j+1));
      }//iv (loss) 
    }
    

    if(cur_size != (int)(_gc.GroupScheme_i.size()-1)){//last grid point

      //right boundary x_{i}+1, absorb the same species
      int tmp_size = std::min(_max_mobile_i-1,_gc.GroupScheme_i_del[cur_size-1]-1);
      for(int i=0;i<=tmp_size;i++){
        conc1 = _u[_qp] + (_gc.GroupScheme_i[cur_size]-i-_gc.GroupScheme_i_avg[cur_size-1])*(*_val_i_vars[2*index+1])[_qp];
        for(int j=0;j<=_max_mobile_i-1-i;j++){
          conc2 = (*_val_i_vars[2*(i+j)])[_qp];
          res_sum += conc1 * conc2 * _gc._absorb(-(_gc.GroupScheme_i[cur_size]-i),-(i+j+1));
        }//ii (loss)
      }
  
      //right boundary x_{i}+1, absorb the opposite species
      tmp_size = std::min(_max_mobile_v-1,_gc.GroupScheme_i_del[cur_size-1]-1);
      for(int i=0;i<=tmp_size;i++){
        int tmp2 = std::min(_max_mobile_v-1-i,_gc.GroupScheme_i.back()-_gc.GroupScheme_i[cur_size]-1);
        for(int j=0;j<=tmp2;j++){
          group_num =  _gc.CurrentGroupI(_gc.GroupScheme_i[cur_size]+j+1);
          //other_size = index+i+j+1;
          other_size = index+(group_num-cur_size);
          conc1 = (*_val_i_vars[2*other_size])[_qp]+(*_val_i_vars[2*other_size+1])[_qp]*(_gc.GroupScheme_i[cur_size]+j+1-_gc.GroupScheme_i_avg[group_num-1]);
          conc2 = (*_val_v_vars[2*(i+j)])[_qp]; 
          res_sum -= conc1 * conc2 * _gc._absorb(-(_gc.GroupScheme_i[cur_size]+j+1),(i+j+1));
        }//iv (gain)
      }

      //right boundary x_{i}+1, emission
      conc = (*_val_i_vars[2*index+2])[_qp]+(_gc.GroupScheme_i[cur_size]+1-_gc.GroupScheme_i_avg[cur_size])*(*_val_i_vars[2*index+3])[_qp];
      res_sum -= conc * _gc._emit(-(_gc.GroupScheme_i[cur_size]+1));//i emit (gain) 
    } 

 
    return 1.0/(_gc.GroupScheme_i_del[cur_size-1])*res_sum *_test[_i][_qp];
  }
}

Real
GImmobileL0::computeQpJacobian()
{
  int cur_size;
  Real jac_sum = 0.0;
  Real conc,conc1,conc2;

  if(_cur_size>0){//v type
    cur_size = _cur_size;

    //left boundary x_{i-1}+1, emission
    conc = 1.0;
    jac_sum += conc * _gc._emit(_gc.GroupScheme_v[cur_size-1]+1);//v emit (loss)


    //left boundary x_{i-1}+1, absorb the opposite species
    for(int i=0;i<=_max_mobile_i-1;i++){
      int tmp_size = std::min(_max_mobile_i-1-i,_gc.GroupScheme_v_del[cur_size-1]-1);
      conc1 = 1.0;
      for(int j=0;j<=tmp_size;j++){
        conc2 = (*_val_i_vars[2*(i+j)])[_qp];
        jac_sum += conc1 * conc2 * _gc._absorb(_gc.GroupScheme_v[cur_size-1]+j+1,-(i+j+1));
      }//vi (loss)
    }

    if(cur_size != (int)(_gc.GroupScheme_v.size()-1)){

      //right boundary x_{i}+1, absorb the same species
      int tmp_size = std::min(_max_mobile_v-1,_gc.GroupScheme_v_del[cur_size-1]-1);
      for(int i=0;i<=tmp_size;i++){
        conc1 = 1.0;
        for(int j=0;j<=_max_mobile_v-1-i;j++){
          conc2 = (*_val_v_vars[2*(i+j)])[_qp];
          jac_sum += conc1 * conc2 * _gc._absorb(_gc.GroupScheme_v[cur_size]-i,i+j+1);
        }//vv (loss)
      }
    } 
    return 1.0/(_gc.GroupScheme_v_del[cur_size-1])*jac_sum *_test[_i][_qp]*_phi[_j][_qp];
  }

  else{
    cur_size = -_cur_size;


    //left boundary x_{i-1}+1, emission
    conc = 1.0;
    jac_sum += conc * _gc._emit(-(_gc.GroupScheme_i[cur_size-1]+1));//i emit (loss)

    //left boundary x_{i-1}+1, absorb the opposite species
    for(int i=0;i<=_max_mobile_v-1;i++){
      int tmp_size = std::min(_max_mobile_v-1-i,_gc.GroupScheme_i_del[cur_size-1]-1);
      conc1 = 1.0;
      for(int j=0;j<=tmp_size;j++){
        conc2 = (*_val_v_vars[2*(i+j)])[_qp];
        jac_sum += conc1 * conc2 * _gc._absorb(-(_gc.GroupScheme_i[cur_size-1]+j+1),(i+j+1));
      }//iv (loss) 
    }
    
    if(cur_size != (int)(_gc.GroupScheme_i.size()-1)){

      //right boundary x_{i}+1, absorb the same species
      int tmp_size = std::min(_max_mobile_i-1,_gc.GroupScheme_i_del[cur_size-1]-1);
      for(int i=0;i<=tmp_size;i++){
        conc1 = 1.0;
        for(int j=0;j<=_max_mobile_i-1-i;j++){
          conc2 = (*_val_i_vars[2*(i+j)])[_qp];
          jac_sum += conc1 * conc2 * _gc._absorb(-(_gc.GroupScheme_i[cur_size]-i),-(i+j+1));
        }//ii (loss)
      }
    } 
    return 1.0/(_gc.GroupScheme_i_del[cur_size-1])*jac_sum *_test[_i][_qp]*_phi[_j][_qp];
  }

}

Real 
GImmobileL0::computeQpOffDiagJacobian(unsigned int jvar){
      return 0.0;
}


int
GImmobileL0::getGroupNumber(std::string str)
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
