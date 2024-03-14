###########################################################
#
#	Makefile for ECE 588 Final Project - Winter '24
#
#	Comments: 
#	 
#	$ make clean    --->   command also creates 
#	placeholder.md, an empty file that ensure the
#	empty directory bin/ does not disappear
#	
#
###########################################################

all:
	cc -Wall -lpthread src/sequential_solver.c -o bin/serial_solver
	cc -Wall -lpthread src/parallel_solver.c -o bin/parallel_solver
	cc -Wall -lpthread src/paralex.c -o bin/paralex

thread:
	cc -Wall -lpthread -fsanitize=thread src/parallel_solver.c -o bin/parallel_solver
	cc -Wall -lpthread -fsanitize=thread src/paralex.c -o bin/paralex

address:
	cc -Wall -lpthread -fsanitize=address src/parallel_solver.c -o bin/parallel_solver
	cc -Wall -lpthread -fsanitize=address src/paralex.c -o bin/paralex

parallel:
	cc -Wall -lpthread src/parallel_solver.c -o bin/parallel_solver

sequential:
	cc -Wall src/sequential_solver.c -o bin/sequential_solver

alex:
	cc -Wall -lpthread src/paralex.c -o bin/paralex

clean:
	rm -r bin/
	mkdir bin
	touch bin/placeholder.md
