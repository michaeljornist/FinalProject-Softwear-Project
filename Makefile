CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LDFLAGS = -lm  # Flags for the linker

.PHONY: clean

# Main executable
symnmf: symnmf.o
	$(CC) -o $@ $^ $(LDFLAGS)

# Object files
symnmf.o: symnmf.c symnmf.h
	$(CC) -c $< $(CFLAGS)

# Clean up
clean:
	rm -f *.o symnmf
