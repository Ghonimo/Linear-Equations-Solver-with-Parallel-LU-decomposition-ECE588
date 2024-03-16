###########################################################
#
#	Makefile for ECE 588 Final Project - Winter '24
#
#	Comments: 
#	 
#	$ make clean    --->   command also creates 
#	placeholder.md, an blank file to ensure the
#	empty directory bin/ does not disappear
#	
#
###########################################################

all:
	cc -Wall src/sequential_solver.c -o bin/sequential_solver
	cc -Wall -lpthread src/parallel_solver.c -o bin/parallel_solver
	cc -Wall -lpthread src/parallel_solver_with_pivoting.c -o bin/parallel_solver_with_pivoting
	cc -Wall -lpthread src/paralex.c -o bin/paralex

linux:
	cc src/sequential_solver.c -o bin/sequential_solver -Wall -std=gnu11
	cc src/parallel_solver.c -o bin/parallel_solver -Wall -lpthread -std=gnu11
	cc src/parallel_solver_with_pivoting.c -o bin/parallel_solver_with_pivoting -Wall -lpthread -std=gnu11
	cc src/paralex.c -o bin/paralex -Wall -lpthread -std=gnu11

thread:
	cc -Wall -lpthread -fsanitize=thread src/parallel_solver.c -o bin/parallel_solver
	cc -Wall -lpthread -fsanitize=thread src/parallel_solver_with_pivoting.c -o bin/parallel_solver_with_pivoting
	cc -Wall -lpthread -fsanitize=thread src/paralex.c -o bin/paralex

address:
	cc -Wall -lpthread -fsanitize=address src/parallel_solver.c -o bin/parallel_solver
	cc -Wall -lpthread -fsanitize=address src/parallel_solver_with_pivoting.c -o bin/parallel_solver_with_pivoting
	cc -Wall -lpthread -fsanitize=address src/paralex.c -o bin/paralex

parallel:
	cc -Wall -lpthread src/parallel_solver.c -o bin/parallel_solver

pivoting:
	cc -Wall -lpthread src/parallel_solver_with_pivoting.c -o bin/parallel_solver_with_pivoting

sequential:
	cc -Wall src/sequential_solver.c -o bin/sequential_solver

alex:
	cc src/paralex.c -Wall -lpthread  -o bin/paralex

clean:
	rm -r bin/
	mkdir bin
	touch bin/placeholder.md
