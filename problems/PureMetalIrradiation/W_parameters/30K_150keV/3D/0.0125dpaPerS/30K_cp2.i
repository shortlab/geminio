#UNITS: um,s,/um^3
#consider only vacancy cluster for tungsten
# implement grouping method

[GlobalParams]
#set the largest size for vacancy clusters and interstitial clusters. Also defined in blocks to be clearer.

  number_v = 50    #number of vacancy variables i.e. total_groups
  max_defect_v_size = 1001  #put in [Global] largest total_groups=max_defect_size-1
  number_single_v = 20  #max size with group size 1
  max_mobile_v = 1

  number_i = 200      #number of interstitial variables, set to 0
  max_defect_i_size = 1001 #put in [Global] largest total_groups=max_defect_size-1
  number_single_i = 45  #max size with group size 1
  max_mobile_i = 10

  temperature = 30  #temperature [K]
  #T_func = T_func
[]

[Mesh]
  type = GeneratedMesh
  xmin = 0
  xmax = 1 #uniform source for simplicity, no spatical dependence
  dim = 1
  nx = 2
[]

# define defect variables, set variables and boundadry condition as 0 where appropriate
[GVariable]
  [./groups]
#    boundary_value = 0.0
    scaling = 1.0  #important factor, crucial to converge
    bc_type = neumann
    #IC_v_size = '1 2 3 4'
    #IC_v = '2000.0 4000.0 1000.0 250.0' #'3.9            2.323' #thermal equil
    IC_v_size = ''
    IC_v = '' #'3.9            2.323' #thermal equil
    #initial concentration for species with value NON-ZERO
    IC_i_size = ''
    IC_i = '' #thermal equil
  [../]
[]

[GTimeDerivative]
  [./groups]
  [../]
[]

[GMobile]
  [./groups]
    group_constant = group_constant
  [../]
[]

[GImmobile]
  [./groups]
    group_constant = group_constant
  [../]
[]

[Sources]
  [./groups]
    source_v_size = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17'
    source_v_value = '287571929 87265746 33890565 17504954 10147668 4276530 3516732 2767298 1047524 785346 677338 284194 146631 57330 21322 94918 21322'
    source_i_size = '1 2 3 4 5 6 7 8 9 10'
    source_i_value = '627673201 55872900 10429603 2427409 700564 522434 306630 20129 67260 3231'
    scaling_factor = 1.0 
  [../]
[]

[AuxVariables]
  [./void_swelling]
  [../]
  [./SIA_density]
  [../]
[]
[GVoidSwelling]
  [./groups]
    aux_var = void_swelling
    group_constant = group_constant
  [../]
[]
[GSumSIAClusterDensity]
#sum up of SIA cluster density in range [lower_bound,upper_bound]
  [./groups]
    aux_var = SIA_density 
    group_constant = group_constant
    lower_bound = 60
  [../]
[]
[Functions]
  [./T_func]
    type = ParsedFunction
    value = '363.0*(t<131579)+773.0*(t>=131579)'
  [../]
[]

[UserObjects]
  [./material]
    type = GTungsten   #definition should be in front of the usage
    i_disl_bias = 1.15
    v_disl_bias = 1.0
    dislocation = 1 #dislocation density 1.0 /um^2
  [../]

  [./group_constant]
    type = GGroup
    material = 'material'
    #GroupScheme = Uniform
    GroupScheme = RSpace
    dr_coef = 0.5
    update = false
    execute_on = initial
  [../]
[]

[Postprocessors]
  [./FluxChecker-V]
    type = NodalVariableValue
    nodeid = 1
    variable = groups0v1
  [../]
  [./FluxChecker-I]
    type = NodalVariableValue
    nodeid = 1
    variable = groups0i1
  [../]
  [./Swelling]
    type = NodalVariableValue
    nodeid = 1
    variable = void_swelling
  [../]
  [./SIADensity]
    type = NodalVariableValue
    nodeid = 1
    variable = SIA_density 
  [../]
[]


#[Preconditioning]
#  active = smp
#  [./smp]
#    type = SMP
#    full = true
#  [../]
#[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient
  solve_type = 'PJFNK'
#  petsc_options =  '-snes_mf_operator'
#  petsc_options_iname =  '-pc_type -pc_hypre_type -ksp_gmres_restart'
#  petsc_options_value =  'hypre    boomeramg  81'
  petsc_options_iname =  '-pc_type -sub_pc_type -ksp_gmres_restart'
  petsc_options_value =  'bjacobi ilu  81'
  #trans_ss_check = true
  #ss_check_tol = 1.0e-14
  l_max_its =  30
  nl_max_its =  40
  nl_abs_tol=  1e-10  #Question: why change to 1e-12 not work!!!
  nl_rel_tol =  1e-7
  l_tol =  1e-8
  num_steps = 500
  start_time = 0
  end_time = 1.116
  #dt = 1.0e-2
  dtmin = 1.0e-10 
  dtmax = 0.01
  active = 'TimeStepper'
  [./TimeStepper]
      cutback_factor = 0.4
      dt = 1e-9
      growth_factor = 2
      type = IterationAdaptiveDT
  [../]
[]

[Debug]
#    show_top_residuals=1
#    show_var_residual_norms=1
[]


[Outputs]
  #output_linear = true
  #file_base = out
  exodus = true
  csv = true
  console = false
[]
