#UNITS: um,s,/um^3
#solve discrete cluster dynamics equation and compare with grouping method in GMIC++
#Two main actions for adding mobile and immobile defects

[GlobalParams]
#set the largest size for vacancy clusters and interstitial clusters. Also defined in blocks to be clearer.
  number_v = 2000
  #max_mobile_v = 2
  mobile_v_size = '1 2'

  number_i = 2000      #number of interstitial variables, set to 0
  max_mobile_i = 3 
  #mobile_i_size = '1 2 3'
  temperature = 573  #temperature [K]
  source_v_size = '1 8'
  source_i_size = '1 4'

[]

[Mesh]
  type = GeneratedMesh
# 1D ion source
  xmin = 0
  xmax = 1 #2000 change to 1, uniform source for simplicity, no spatical dependence
  dim = 1
  nx = 1
[]


[LotsOfVariables]
# define defect variables, set variables and boundadry condition as 0
  [./defects]
#    boundary_value = 0.0
    scaling = 1.0  #important factor, crucial to converge
    bc_type = neumann
  [../]
[]

[LotsOfFunction]
#checked: no problem in adding many functions here
  [./func_defect_] 
#for now only 1D, and data file should be in column, postion+v_size+i_size; defectsi1 should include injected interstitials. Make sure unit is in um^3 or um
#data should be the same with [sources] block and the sub_block name [func_defect_]
    type = PiecewiseLinearTimeLimit
    tlimit = 2.0e8 #limit the time with sources [s]
    data_file = data.txt
    axis = 0 #x axis
    format = columns
    scale_factor = 1.0
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
    func_pre_name = func_defect_
    # either provide data file to construct functions or direct input of constant source for each size.
  [../]
[]

[UserObjects]
  [./material]
    type = TestProperty   #definition should be in front of the usage
    i_disl_bias = 1.1
    v_disl_bias = 1.0
    dislocation = 1.0 #dislocation density 1.0 /um^2
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
#  [./Total_Number_v]
#    type = NodalConservationCheck
#    nodeid = 1
#    var_prefix = 'defectsv'
#    size_range = '1 8'
#  [../]
#  [./Total_Number]
#    type = NodalConservationCheck
#    nodeid = 1
#    var_prefix = 'defectsv'
#    size_range = '1 1000'
#  [../]
[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient
  scheme = bdf2
  solve_type = 'PJFNK'
#  petsc_options =  '-snes_mf_operator'
  petsc_options_iname =  '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value =  'hypre    boomeramg  41'
#  petsc_options_iname =  '-pc_type -sub_pc_type -ksp_gmres_restart'
#  petsc_options_value =  'bjacobi ilu  81'
  l_max_its = 40
  nl_max_its =  40
  nl_rel_tol =  1e-6
  l_tol =  1e-6
  num_steps = 400
  start_time = 0
  end_time = 500
 # dt = 1e0
  dtmin = 1.0e-10 
  #dtmax = 100
  #dt = 1.0e-2
  active = 'TimeStepper'

  [./TimeStepper]
      cutback_factor = 0.4
      dt = 1e-9
      growth_factor = 2
      type = IterationAdaptiveDT
  [../]
  [./TimeStepper0]
      type                     = DT2
      dt                       = 1e-9                        # The initial time step size.
      e_max                    = 1e7                  # Maximum acceptable error.
      e_tol                    = 1e6                  # Target error tolerance.
      max_increase             = 1.5                       # Maximum ratio that the time step can increase.
  [../]
[]

[Debug]
  #show_top_residuals=1
  #show_var_residual_norms=1
[]

[Outputs]
  #output_linear = true
  #file_base = out
  exodus = true
  csv = true
  console = true
  print_perf_log = true
[]
