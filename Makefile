all: poisson

# -g outputs debugging information
# -Wall enables all warnings
# -pthread configures threading
#CC = gcc
CFLAGS = -O3 -ffast-math -funroll-loops -pthread

poisson: poisson.c

.PHONY: disassembly
disassembly: poisson.s

poisson.s: poisson
	objdump -S --disassemble $< > $@

.PHONY: test
test: poisson
	./test.sh

.PHONY: clean
clean:
	rm -f poisson *.o *.s
