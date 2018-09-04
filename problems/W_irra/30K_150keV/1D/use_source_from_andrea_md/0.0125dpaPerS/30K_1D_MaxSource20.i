#UNITS: um,s,/um^3
# implement grouping method

[GlobalParams]
#set the largest size for vacancy clusters and interstitial clusters. Also defined in blocks to be clearer.

  number_v = 150    #number of vacancy variables i.e. total_groups
  number_single_v = 101  #max size with group size 1
  max_mobile_v = 1

  number_i = 250      #number of interstitial variables, set to 0
  number_single_i = 101  #max size with group size 1
  max_mobile_i = 5

  temperature = 30  #temperature [K]
  SIAMotionDim = 1D
  aux_prefix = rLambda
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
    source_i_size = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20'
    source_i_value = '304254095 57645371 21784390 10921755 6393054 4127370 2851025 2069285 1559747 1857541 1564702 1337864 1158351 1013696 895311 797117 714710 644834 585035 533438'
    source_v_size = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20'
    source_v_value = '485628689 60703586 17986247 7587948 3885029 2248280 1415827 948493 666157 480533 408654 352464 307623 271210 241195 216132 194966 176913 161377 147901'
    scaling_factor = 0.5
  [../]
[]

[RecipMeanFreePath]
  [./groups]
    group_constant = group_constant
  [../]
[]

[UserObjects]
  [./material]
    type = GTungsten1D   #definition should be in front of the usage
    i_disl_bias = 1.15
    v_disl_bias = 1.0
    dislocation = 1 #dislocation density 1.0 /um^2
  [../]

  [./group_constant]
    type = GGroup
    material = 'material'
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
