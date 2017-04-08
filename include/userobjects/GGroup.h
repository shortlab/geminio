/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/


#ifndef GGROUP_H
#define GGROUP_H

#include "GMaterialConstants.h"
#include "Function.h"
#include "GeneralUserObject.h"

class GGroup;
class Function;


template<>
InputParameters validParams<GGroup>();

/**
 * Saturation of a phase as a function of
 * effective saturation of that phase,
 * and its derivatives wrt effective saturation
 */
class GGroup : public GeneralUserObject
{
public:
  GGroup(const InputParameters & parameters);
  ~GGroup();

  void initialize();
  void execute();
  void finalize();

  void setGroupScheme();
  void updateGroupScheme();
  int CurrentGroupV(int) const;
  int CurrentGroupI(int) const;

  Real _emit(int) const;//return kth group constant based on single shape function
  Real _disl(int) const;//return dislocation sink strenght based on shape function
  Real _diff(int) const;//return diffusion coefficient based on shape function
  Real _absorb(int,int) const;//return kth,jth group constant based on double shape functions

  std::vector<int> GroupScheme_v;
  std::vector<int> GroupScheme_i;
  Real* GroupScheme_v_sq;//dispersion
  Real* GroupScheme_i_sq;//dispersion
  Real* GroupScheme_v_avg;//group avg
  Real* GroupScheme_i_avg;//group avg
  int* GroupScheme_v_del;//group del
  int* GroupScheme_i_del;//group del
  Real _atomic_vol;


protected:
  MooseEnum _GroupScheme;
  Real _dr_coef;
  int _Ng_v;
  int _Ng_i;
  int _num_v;
  int _num_i;
  int _v_size;
  int _i_size;
  int _single_v_group;
  int _single_i_group;
  Real _T;
  Function * const _T_func;
  bool _update;
//  Real* _emit_array;//total _Ng_v + _Ng_i
//  Real** _absorb__matrix;//(_Ng_v+_Ng_i)xtotal_no_of_mobile_species
  bool _has_material;
  const GMaterialConstants * const _material;
  Point dummy;
};

#endif // 
