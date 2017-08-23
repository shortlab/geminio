
//reference: Cluster Dynamics Models of Irradiation Damage Accumulation in Ferritic Iron Part II: Effects of Reaction Dimensionality
//calculate according to Eq (8)
#include "ReciprocalMeanFreePath1D.h"

template<>
InputParameters validParams<ReciprocalMeanFreePath1D>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("coupled_v_vars","coupled vacancy type variables");
  params.addRequiredCoupledVar("coupled_i_vars","coupled sia type variables");
  params.addParam<int>("mobile_SIA_size",1,"calculate results for this size");
  params.addRequiredParam<UserObjectName>("user_object","The name of user object providing interaction constants");
  return params;
}


ReciprocalMeanFreePath1D::ReciprocalMeanFreePath1D(const
                                   InputParameters & parameters)
  :AuxKernel(parameters),
  _gc(getUserObject<GGroup>("user_object")),
  _cur_size(getParam<int>("mobile_SIA_size"))
{
  int nvcoupled = coupledComponents("coupled_v_vars");
  _val_v_vars.resize(nvcoupled);

  for (int i=0; i < nvcoupled; ++i)
  {
    _val_v_vars[i] = &coupledValue("coupled_v_vars",i);
  }

  int nicoupled = coupledComponents("coupled_i_vars");
  _val_i_vars.resize(nicoupled);

  for (int i=0; i < nicoupled; ++i)
  {
    _val_i_vars[i] = &coupledValue("coupled_i_vars",i);
  }

  if (_cur_size > _gc._i_size)//larger than maximum mobile SIA size
    mooseError("Current kernel is for mobile SIA clusters");

}

Real
ReciprocalMeanFreePath1D::computeValue()
{
  Real recip = 0.0;

  for(int i=0;i<_gc._Ng_v;i++){//loop each group
     Real L0 = (*_val_v_vars[i])[_qp];
     Real del = _gc.GroupScheme_v_del[i];
     recip += L0*del*_gc._absorb(-_cur_size,i+1);//means sigma_i,j
  }
  for(int i=0;i<_gc._Ng_i;i++){//loop each group
     Real L0 = (*_val_i_vars[i])[_qp];
     Real del = _gc.GroupScheme_i_del[i];
     recip += L0*del*_gc._absorb(-_cur_size,-(i+1));
  }

  return recip;
}

