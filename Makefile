all: bvm cvm

bvm: bvm.c hypocrisy.h
	gcc -Wall -o bvm bvm.c -lgsl

cvm: cvm.c hypocrisy.h
	gcc -Wall -o cvm cvm.c -lgsl


