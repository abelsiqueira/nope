#DEBUG = -ggdb
all:
	gcc -c -o preprocessor.o *.c $(DEBUG) -Wall -Wextra
