

all:
	gcc -Wall src/main.c -o bin/serial_solver

serial:
	gcc -Wall src/main.c -o bin/serial_solver

clean:
	rm bin/serial_solver
