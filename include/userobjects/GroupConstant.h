/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/


#ifndef GROUPCONSTANT_H
#define GROUPCONSTANT_H

#include "GMaterialConstants.h"


#include "GeneralUserObject.h"

class GroupConstant;


template<>
InputParameters validParams<GroupConstant>();

/**
 * Saturation of a phase as a function of
 * effective saturation of that phase,
 * and its derivatives wrt effective saturation
 */
class GroupConstant : public GeneralUserObject
{
public:
  GroupConstant(const InputParameters & parameters);
  ~GroupConstant();

  void initialize();
  void execute();
  void finalize();

  void setGroupScheme();
  void setGroupConstant();
  void updateGroupScheme();

  Real emit_gc(int,int);//emission group constant
  Real disl_gc(int,int);//dislocation group constant
  Real diff_gc(int,int);//diffusion group constant
  Real absorb_gc(int,int,int,int);//absorption group constant
  Real _emit(int) const;//return kth group constant based on single shape function
  Real _disl(int) const;//return dislocation sink strenght based on shape function
  Real _diff(int) const;//return diffusion coefficient based on shape function
  Real _absorb(int,int) const;//return kth,jth group constant based on double shape functions

protected:
  MooseEnum _GroupScheme;
  Real _sigma;
  Real _boosting_factor;
  typedef std::map<std::pair<int,int> , Real> DoubleKey;// +: vacancy -:interstitial
  typedef std::map<int,Real> SingleKey;// +: vacancy -:interstitial
  int _Ng_v;
  int _Ng_i;
  int _num_v;
  int _num_i;
  std::vector<int> GroupScheme_v;
  std::vector<int> GroupScheme_i;
  Real _atomic_vol;

  std::vector<int> _v_size;
  std::vector<int> _i_size;
  int _single_v_group;
  int _single_i_group;
  Real _T;
  bool _update;
//  Real* _emit_array;//total _Ng_v + _Ng_i
//  Real** _absorb__matrix;//(_Ng_v+_Ng_i)xtotal_no_of_mobile_species
  SingleKey _emit_array;
  SingleKey _disl_array;
  SingleKey _diff_array;
  DoubleKey _absorb_matrix;  

  bool _has_material;
  const GMaterialConstants * const _material;
};

#endif // 
