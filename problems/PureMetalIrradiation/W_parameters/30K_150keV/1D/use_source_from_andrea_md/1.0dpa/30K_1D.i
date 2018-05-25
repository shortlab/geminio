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
    source_i_size = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50'
    source_i_value = '239836283.83 45440478.58 17172118.43 8609360.77 5039493.34 3253508.05 2247395.93 1631168.84 1229512.26 1464255.82 1233417.55 1054607.16 913101.22 799073.07 705752.93 628349.09 563389.7 508307.42 461169.39 420497.06 385143.19 354206.18 326969.21 302856.38 281400.55 262219.47 244997.8 229473.48 215427.21 202674.31 191058.3 180445.89 170722.93 161791.2 153565.79 145972.97 138948.5 132436.13 126386.49 120756.07 115506.39 110603.34 106016.57 101719.01 97686.47 93897.25 90331.83 86972.66 83803.86 80811.09' 
    source_v_size = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50'
    source_v_value = '450394700.71 56299337.59 16681285.21 7037417.2 3603157.61 2085160.65 1313104.08 879677.15 617825.38 445669.31 379005.19 326892.24 285304.77 251533.05 223695.67 200451.24 180821.29 164077.67 149668.91 137170.82 126252.48 116652.53 108162.17 100612.89 93867.38 87812.84 82355.83 77418.38 72935.02 68850.42 65117.52 61696.11 58551.64 55654.28 52978.18 50500.83 48202.53 46066.01 44076.07 42219.27 40483.7 38858.76 37335.02 35904.03 34558.23 33290.81 32095.66 30967.25 29900.56 28891.08' 
    scaling_factor = 1.0 
  [../]
[]

[RecipMeanFreePath]
  [./groups]
    group_constant = group_constant
  [../]
[]

[AuxVariables]
  [./SIA_density]
  [../]
[]
[GSumSIAClusterDensity]
#sum up of SIA cluster density in range [lower_bound,upper_bound]
  [./groups]
    aux_var = SIA_density 
    group_constant = group_constant
    lower_bound = 2
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
  end_time = 80.0
  #dt = 1.0e-2
  dtmin = 1.0e-10 
  dtmax = 0.5
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
