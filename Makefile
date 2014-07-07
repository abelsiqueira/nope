DEBUG = -ggdb

C = gcc
CFLAGS = -Wall -Wextra $(DEBUG)
OBJS = preprocessor.o functions.o

all: $(OBJS)
	ar rv libprep.a $(OBJS)

%.o: %.c
	$(C) -c -o $@ $< $(CFLAGS)

clean:
	rm -f $(OBJS) *.a
