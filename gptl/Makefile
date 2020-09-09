all:
	mpicc test.c -o a.out -I./include -L./lib -lgptl  
run:
	bsub -I -b -q q_sw_expr  -n 4 ./a.out
clean:
	rm *.o
