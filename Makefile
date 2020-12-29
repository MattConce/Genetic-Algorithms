OP1 = -Wall -Wpedantic

LM = -lm

CC = gcc

VER = -std=c11

SRC = lib/*.c

INC = -I./lib

CC_OMP=-fopenmp

BIN = bin

.PHONY : run clear

default :
	@echo 'To complile prog.cc type make prog'
	@echo 'To run it type make run-prog'
	@echo ''

% : %.c
	$(CC) $(OP1) $(VER) $(CC_OMP) $< -o $(BIN)/$@ $(LM)

run% :
	./$(BIN)/$(subst run-,'',$@)

clean :
	rm -r bin/*

