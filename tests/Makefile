include ../make.inc

SRCS = $(wildcard *.c)
TSTS = $(SRCS:.c=.tst)

unit-test: unit-test-trail-whitespaces unit-tsts

unit-test-trail-whitespaces:
	@./test-trailing-whitespaces.sh

unit-tsts: $(TSTS)
	for tst in $(TSTS); do ./$$tst || exit 1; done

%.tst: %.c
	$(CC) -o $@ $< ../lib$(PACKAGE).a $(CFLAGS) $(I)

clean:
	$(RM) $(TSTS)
