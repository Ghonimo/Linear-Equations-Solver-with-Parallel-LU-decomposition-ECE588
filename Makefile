###########################################################
#
#	Makefile for ECE 588 Final Project - Winter '24
#
#	Comments: Compile commands for parallel2.c are
#	commented out because it throws warning and 
#	prematurely exits. 'make clean' creates 
#	placeholder.md, an empty file that ensure the
#	empty directory bin/ does not disappear
#	
#
###########################################################

all:
	cc src/main.c -Wall -lpthread -o bin/serial_solver
	cc src/main_parallel.c -Wall -lpthread -o bin/parallel_solver
	cc src/parallel1.c  -Wall -lpthread -o bin/parallel1
	cc src/paralex.c -Wall -lpthread -o bin/paralex
	cc src/lu_parallel.c  -Wall -lpthread -o bin/lu_parallel

thread:
	cc src/main_parallel.c -Wall -lpthread -fsanitize=thread -o bin/parallel_solver
	cc src/parallel1.c  -Wall -lpthread -fsanitize=thread -o bin/parallel1
	cc src/paralex.c -Wall -lpthread -fsanitize=thread -o bin/paralex
	cc src/lu_parallel.c  -Wall -lpthread -fsanitize=thread -o bin/lu_parallel

address:
	cc src/main_parallel.c -Wall -lpthread -fsanitize=address -o bin/parallel_solver
	cc src/parallel1.c  -Wall -lpthread -fsanitize=address -o bin/parallel1
	cc src/paralex.c -Wall -lpthread -fsanitize=address -o bin/paralex
	cc src/lu_parallel.c  -Wall -lpthread -fsanitize=address -o bin/lu_parallel

parallel:
	cc src/main_parallel.c -Wall -lpthread -o bin/parallel_solver

parallel_lu:
	cc src/lu_parallel.c -Wall -lpthread -o bin/lu_parallel

p1:
	cc src/parallel1.c  -Wall -lpthread  -o bin/parallel1

alex:
	cc src/paralex.c -Wall -lpthread  -o bin/paralex

serial:
	cc src/main.c -Wall -lpthread  -o bin/serial_solver

clean:
	rm -r bin/
	mkdir bin
	touch bin/placeholder.md
