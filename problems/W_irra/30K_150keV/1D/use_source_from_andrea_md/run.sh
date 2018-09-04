
######### directory 0.0125dpaPerS
dir=0.0125dpaPerS
cd ${dir}
mobiles=(3 5 6 7 8 10 15)
for size in ${mobiles[@]};
do
  mpirun ~/projects/geminio-dev/Geminio-opt -i 30K_1D.i GlobalParams/max_mobile_i=${size} Outputs/file_base=30K_1D_mobile${size}_out
done
cd ..

######### directory 0.025dpaPerS
dir=0.025dpaPerS
cd ${dir}
mobiles=(5 6 7)
for size in ${mobiles[@]};
do
  mpirun ~/projects/geminio-dev/Geminio-opt -i 30K_1D.i GlobalParams/max_mobile_i=${size} Outputs/file_base=30K_1D_mobile${size}_out
done
cd ..

######### directory 0.05dpaPerS
dir=0.05dpaPerS
cd ${dir}
mobiles=(5 6 7)
for size in ${mobiles[@]};
do
  mpirun ~/projects/geminio-dev/Geminio-opt -i 30K_1D.i GlobalParams/max_mobile_i=${size} Outputs/file_base=30K_1D_mobile${size}_out
done
cd ..


######### directory 1.0dpa
dir=1.0dpa
cd ${dir}
mobiles=(5 6 7)
for size in ${mobiles[@]};
do
  mpirun ~/projects/geminio-dev/Geminio-opt -i 30K_1D.i GlobalParams/max_mobile_i=${size} Outputs/file_base=30K_1D_mobile${size}_out
done
cd ..

