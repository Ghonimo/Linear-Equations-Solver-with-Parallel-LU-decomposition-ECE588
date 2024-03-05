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
	cc -Wall -lpthread src/main.c -o bin/serial_solver
	cc -Wall -lpthread src/main_parallel.c -o bin/parallel_solver
	cc -Wall -lpthread src/parallel1.c -o bin/parallel1
	#cc -Wall -lpthread src/parallel2.c -o bin/parallel2
	cc -Wall -lpthread src/parallel3.c -o bin/parallel3
	cc -Wall -lpthread src/paralex.c -o bin/paralex

thread:
	cc -Wall -lpthread -fsanitize=thread src/main_parallel.c -o bin/parallel_solver
	cc -Wall -lpthread -fsanitize=thread src/parallel1.c -o bin/parallel1
	#cc -Wall -lpthread -fsanitize=thread src/parallel2.c -o bin/parallel2
	cc -Wall -lpthread -fsanitize=thread src/parallel3.c -o bin/parallel3
	cc -Wall -lpthread -fsanitize=thread src/paralex.c -o bin/paralex

address:
	cc -Wall -lpthread -fsanitize=address src/main_parallel.c -o bin/parallel_solver
	cc -Wall -lpthread -fsanitize=address src/parallel1.c -o bin/parallel1
	#cc -Wall -lpthread -fsanitize=address src/parallel2.c -o bin/parallel2
	cc -Wall -lpthread -fsanitize=address src/parallel3.c -o bin/parallel3
	cc -Wall -lpthread -fsanitize=address src/paralex.c -o bin/paralex

parallel:
	cc -Wall -lpthread src/main_parallel.c -o bin/parallel_solver

p1:
	cc -Wall -lpthread src/parallel1.c -o bin/parallel1

p2:
	cc -Wall -lpthread src/parallel2.c -o bin/parallel2

p3:
	cc -Wall -lpthread src/parallel3.c -o bin/parallel3

alex:
	cc -Wall -lpthread src/paralex.c -o bin/paralex

serial:
	cc -Wall src/main.c -o bin/serial_solver

clean:
	rm -r bin/
	mkdir bin
	touch bin/placeholder.md
