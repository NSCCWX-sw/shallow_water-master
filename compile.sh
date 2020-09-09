	sw5CC -host -c  main.cpp Grid.cpp -I/usr/sw-mpp/apps/bin/WRF/netcdf/include -L/usr/sw-mpp/apps/Logiciels/netcdf-4.1.2/lib/ -lnetcdf_c++ -L/usr/sw-mpp/sw-mpp/apps/bin/WRF/netcdf/lib/ -lnetcdf -lslave
	# mpiCC master.o slave.o -lstdc++ -o main.x
	#mpiCC -hybrid *.o -o main.x
