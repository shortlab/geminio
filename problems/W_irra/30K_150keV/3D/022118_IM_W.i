#UNITS: um,s,/um^3
#solve discrete cluster dynamics equation and compare with grouping method in GMIC++
#Two main actions for adding mobile and immobile defects

[GlobalParams]
#set the largest size for vacancy clusters and interstitial clusters. Also defined in blocks to be clearer.
  number_v = 1000
  max_mobile_v = 1
  max_defect_v_size = 1001 #put in [Global] largest total_groups=max_defect_size-1
  mobile_v_size = '1 2 3 4 5 6 7 8 9 10'

  number_i = 10000      #number of interstitial variables, set to 0
  max_mobile_i = 10
  max_defect_i_size = 10001 #put in [Global] largest total_groups=max_defect_size-1
  mobile_i_size = '1 2 3 4 5 6 7 8 9 10'

  temperature = 30  #temperature [K]

[]

[Mesh]
  type = GeneratedMesh
# 1D ion source
  xmin = 0
  xmax = 1 #2000 change to 1, uniform source for simplicity, no spatical dependence
  dim = 1
  nx = 2
[]


[LotsOfVariables]
# define defect variables, set variables and boundadry condition as 0
  [./defects]
    boundary_value = 0.0
    scaling = 1.0  #important factor, crucial to converge
    bc_type = neumann
  [../]
[]


[ClusterIC]
  [./defects]
    IC_v_size = '1'
    IC_i_size = ''
    #initial concentration for species with value NON-ZERO
    IC_v = '0.0'        #'3.9            2.323'
    IC_i = ''   #'2.79            2.34        0.57' #corresponds to IC_i_size
  [../]
[]

[MobileDefects]
  [./defects]
    group_constant = group_constant
  [../]
[]
[ImmobileDefects]
  [./defects]
    group_constant = group_constant
  [../]
[]

[LotsOfTimeDerivative]
  [./defects]
  [../]
[]

[LotsOfSource]
  [./defects]
    custom_input = true
    source_i_size = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50'
    source_i_value = '239836283 45440478 17172118 8609360 5039493 3253508 2247395 1631168 1229512 1464255 1233417 1054607 913101 799073 705752 628349 563389 508307 461169 420497 385143 354206 326969 302856 281400 262219 244997 229473 215427 202674 191058 180445 170722 161791 153565 145972 138948 132436 126386 120756 115506 110603 106016 101719 97686 93897 90331 86972 83803 80811'
    source_v_size = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50'
    source_v_value = '450394700 56299337 16681285 7037417 3603157 2085160 1313104 879677 617825 445669 379005 326892 285304 251533 223695 200451 180821 164077 149668 137170 126252 116652 108162 100612 93867 87812 82355 77418 72935 68850 65117 61696 58551 55654 52978 50500 48202 46066 44076 42219 40483 38858 37335 35904 34558 33290 32095 30967 29900 28891'
    # either provide data file to construct functions or direct input of constant source for each size.
    scaling_factor = 0.5
  [../]
[]

[LotsOfCoeffDiffusion]
  [./defects]
    custom_input = false
  [../]
[]

[UserObjects]
  [./material]
    type = GTungsten   #definition should be in front of the usage
    i_disl_bias = 1.15
    v_disl_bias = 1.0
    dislocation = 1 #dislocation density  /um^2
  [../]

  [./group_constant]
    type = GroupConstant
    material = 'material'
    update = false
    execute_on = initial
  [../]
[]


[Postprocessors]
  [./FluxChecker-V]
    type = NodalVariableValue
    nodeid = 1
    variable = defectsv1
  [../]
  [./FluxChecker-I]
    type = NodalVariableValue
    nodeid = 1
    variable = defectsi1
  [../]
#  [./Total_Number]
#    type = NodalConservationCheck
#    nodeid = 1
#    var_prefix = 'defectsv'
#    size_range = '1 1000'
#  [../]
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

#[Debug]
#    show_top_residuals=1
#    show_var_residual_norms=1
#[]

[Outputs]
  #output_linear = true
#  file_base = out
  exodus = true
  csv = true
  console = true
[]
