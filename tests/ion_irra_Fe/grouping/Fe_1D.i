#UNITS: um,s,/um^3
#IRON
# implement grouping method

[GlobalParams]
#set the largest size for vacancy clusters and interstitial clusters. Also defined in blocks to be clearer.

  number_v = 120    #number of vacancy variables i.e. total_groups
  number_single_v = 20  #max size with group size 1
  max_mobile_v = 2

  number_i = 120      #number of interstitial variables i.e. total_groups
  number_single_i = 20  #max size with group size 1
  max_mobile_i = 4

  temperature = 723  #temperature [K]
  source_v_size = '1 2 3 4 5 6 7 8 9 10' 
  source_i_size = '1 2 3 4 5 6 7 8 9 10'
  SIAMotionDim = 1D
  aux_prefix = rLambda
[]

[Mesh]
  type = GeneratedMesh
  xmin = 0
  xmax = 1.6 #um with spatical dependence
  dim = 1
  nx = 50
[]

# define defect variables, set variables and boundadry condition as 0 where appropriate
[GVariable]
  [./groups]
    boundary_value = 0.0
    scaling = 1.0  #important factor, crucial to converge
    bc_type = dirichlet
    #IC_v_size = '1 2 3 4'
    #IC_v = '2000.0 4000.0 1000.0 250.0' #'3.9            2.323' #thermal equil
    IC_v_size = ''
    IC_v = '' #'3.9            2.323' #thermal equil
    #initial concentration for species with value NON-ZERO
    IC_i_size = ''
    IC_i = '' #thermal equil
  [../]
[]

[AuxVariables]
  [./void_swelling]
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

[LotsOfFunction]
  [./func_defect_] 
#for now only 1D, and data file should be in column, postion+v_size+i_size; defectsi1 should include injected interstitials. Make sure unit is in um^3 or um
#data should be the same with [sources] block and the sub_block name [func_defect_]
    type = PiecewiseLinearTimeLimit
    tlimit = 2.0e8 #limit the time with sources [s]
    data_file = spatial_defect_cluster_production_powerlaw_v10i10.txt
    axis = 0 #x axis
    format = columns
    scale_factor = 1.0e-2 #efficiency
  [../]
[]
[LotsOfSource]
  [./groups0] 
  #Attention: add 0 (e.g. groups0) after defects for L0 term
    func_pre_name = func_defect_
    # either provide data file to construct functions or direct input of constant source for each size.
  [../]
[]

[LotsOfUserObjectDiffusion]
  [./groups0]
  #Attention: add 0 (e.g. groups0) after defects for L0 term
    group_constant = group_constant
  [../]
[]

[RecipMeanFreePath]
  [./groups]
    group_constant = group_constant
  [../]
[]

[UserObjects]
  [./material]
    type = GIron1D   #definition should be in front of the usage
    i_disl_bias = 1.15
    v_disl_bias = 1.0
    dislocation = 100 #dislocation density 1.0 /um^2
  [../]

  [./group_constant]
    type = GGroup
    material = 'material'
    GroupScheme = RSpace
    SIAMotionDim = 1D
    dr_coef = 0.5
    update = false
    execute_on = initial
  [../]
[]

[GVoidSwelling]
  [./groups]
    aux_var = void_swelling
    group_constant = group_constant
    #lower_bound = 152 #1.5nm in diameter
    lower_bound = 360 #2.0nm in diameter
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
  nl_abs_tol=  1e-10  
  nl_rel_tol =  1e-7
  l_tol =  1e-8
  num_steps = 1000
  start_time = 0
  end_time = 1.0e2  #certain dpa with dose rate defined by source term
  dtmin = 1.0e-10 
  dtmax = 50
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
  #file_base = out
  exodus = true
  csv = true
  console = true
[]
