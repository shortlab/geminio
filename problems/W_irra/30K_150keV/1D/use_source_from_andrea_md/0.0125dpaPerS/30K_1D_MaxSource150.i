#UNITS: um,s,/um^3
# implement grouping method

[GlobalParams]
#set the largest size for vacancy clusters and interstitial clusters. Also defined in blocks to be clearer.

  number_v = 200    #number of vacancy variables i.e. total_groups
  number_single_v = 151  #max size with group size 1
  max_mobile_v = 1

  number_i = 300      #number of interstitial variables, set to 0
  number_single_i = 151  #max size with group size 1
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
    source_v_size ='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150' 
    source_v_value = '402611279 50326409 14911528 6290801 3220890 1863941 1173793 786350 552278 398387 338795 292211 255036 224847 199963 179184 161637 146670 133790 122618 112858 104276 96686 89938 83908 78496 73618 69204 65197 61545 58209 55150 52339 49749 47357 45143 43088 41178 39399 37740 36188 34736 33374 32094 30891 29758 28690 27681 26728 25825 24971 24160 23390 22658 21962 21300 20668 20066 19492 18942 18418 17915 17435 16974 16533 16109 15702 15312 14936 14576 14228 13894 13572 13262 12962 12674 12395 12126 11866 11615 11373 11138 10911 10691 10478 10272 10072 9878 9690 9508 9331 9159 8992 8830 8673 8520 8371 8226 8085 7948 7815 7685 7559 7436 7316 7199 7085 6974 6865 6759 6656 6555 6457 6361 6267 6176 6086 5999 5913 5830 5748 5668 5590 5514 5439 5366 5294 5224 5155 5088 5022 4958 4895 4833 4772 4712 4654 4597 4541 4486 4432 4379 4327 4276 4226 4177 4129 4081 4035 3989'
    source_i_size ='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150' 
    source_i_value = '182629572 34601833 13076155 6555821 3837453 2477468 1711338 1242095 936244 1114995 939217 803058 695304 608474 537413 478472 429007 387063 351169 320198 293277 269719 248979 230617 214279 199673 186559 174738 164042 154331 145486 137405 130001 123200 116936 111154 105805 100846 96240 91952 87955 84221 80729 77456 74385 71500 68785 66227 63814 61535 59380 57341 55408 53575 51834 50180 48606 47108 45681 44320 43020 41779 40593 39459 38373 37333 36336 35379 34462 33581 32734 31920 31138 30384 29659 28960 28287 27637 27011 26406 25822 25258 24713 24186 23676 23183 22705 22243 21795 21361 20941 20533 20137 19753 19380 19018 18667 18325 17994 17671 17357 17052 16755 16466 16185 15911 15645 15385 15132 14885 14645 14410 14181 13958 13740 13528 13321 13118 12920 12727 12538 12354 12174 11998 11825 11657 11492 11331 11174 11019 10868 10721 10576 10434 10296 10160 10027 9896 9768 9643 9520 9400 9282 9166 9053 8942 8832 8725 8620 8517' 
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
