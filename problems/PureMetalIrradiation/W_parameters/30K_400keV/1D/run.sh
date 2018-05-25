
######### directory 0.0125dpaPerS
dir=0.014dpa0.0125dpaPerS
cd ${dir}
mobiles=(5 6 7)
for size in ${mobiles[@]};
do
  mpirun ~/projects/geminio-dev/Geminio-opt -i 30K_1D.i GlobalParams/max_mobile_i=${size} Outputs/file_base=30K_1D_mobile${size}_out
done
cd ..

######### directory 0.025dpaPerS
dir=0.014dpa0.016dpaPerS
cd ${dir}
mobiles=(5 6 7)
for size in ${mobiles[@]};
do
  mpirun ~/projects/geminio-dev/Geminio-opt -i 30K_1D.i GlobalParams/max_mobile_i=${size} Outputs/file_base=30K_1D_mobile${size}_out
done
cd ..

######### directory 0.05dpaPerS
dir=0.018dpa0.016dpaPerS
cd ${dir}
mobiles=(6)
for size in ${mobiles[@]};
do
  mpirun ~/projects/geminio-dev/Geminio-opt -i 30K_1D.i GlobalParams/max_mobile_i=${size} Outputs/file_base=30K_1D_mobile${size}_out
done
cd ..

