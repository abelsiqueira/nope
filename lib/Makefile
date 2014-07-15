include ../make.inc

SRCS = $(wildcard ../src/*.c)
OBJS = $(notdir $(SRCS:.c=.o))

all: $(OBJS)
	$(AR) lib$(PACKAGE).a $(OBJS)
	$(MV) lib$(PACKAGE).a ../

%.o: ../src/%.c
	$(CC) -c -o $@ $< $(CFLAGS) $(I)

clean:
	$(RM) $(OBJS) ../lib$(PACKAGE).a
