

all:
	gcc -Wall src/main.c -o bin/serial_solver
	gcc -Wall src/main_parallel.c -o bin/parallel_solver
	gcc -Wall src/parallel1.c -o bin/parallel1
	gcc -Wall src/parallel2.c -o bin/parallel2
	gcc -Wall src/parallel3.c -o bin/parallel3

parallel:
	gcc -Wall src/main_parallel.c -o bin/parallel_solver

p1:
	gcc -Wall src/parallel1.c -o bin/parallel1

p2:
	gcc -Wall src/parallel2.c -o bin/parallel2

p3:
	gcc -Wall src/parallel3.c -o bin/parallel3

serial:
	gcc -Wall src/main.c -o bin/serial_solver

clean:
	rm bin/*
