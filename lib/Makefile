include ../make.inc

SRCS = $(wildcard ../src/*.c)
OBJS = $(notdir $(SRCS:.c=.o))

all: $(OBJS)
	$(AR) ../lib$(PACKAGE).a $(OBJS)

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	$(RM) $(OBJS) ../lib$(PACKAGE).a
