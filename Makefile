main.x : master.o slave.o
	mpiCC -hybrid *.o -o main.x -L/usr/sw-mpp/apps/Logiciels/netcdf-4.1.2/lib/ -lnetcdf_c++ -L/usr/sw-mpp/sw-mpp/apps/bin/WRF/netcdf/lib/ -lnetcdf -L/home/export/online3/lyf1527/yanxu2/shallow_water2D-master/gptl/lib -lgptl

master.o: Grid.o slave.o main.cpp
	sw5CC -host -c -o master.o main.cpp Grid.o -I/usr/sw-mpp/mpi2/include -I/usr/sw-mpp/apps/bin/WRF/netcdf/include -I/home/export/online3/lyf1527/yanxu2/shallow_water2D-master/gptl/include

Grid.o: Grid.cpp
	sw5CC -host -c Grid.cpp -OPT:IEEE_arith=1 -I/usr/sw-mpp/mpi2/include -I/usr/sw-mpp/apps/bin/WRF/netcdf/include -I/home/export/online3/lyf1527/yanxu2/shallow_water2D-master/gptl/include

slave.o: slave.c
	sw5cc -slave -c -msimd -O0 slave.c

run:
	bsub -b -I -q q_sw_expr -n 9 -cgsp 64 -share_size 7000 -host_stack 128 ./main.x 

clean:
	rm -f *.o *.x
	rm -f ./eta/*
	rm -f ./png/*
