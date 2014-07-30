include ../make.inc

SRCS = $(wildcard *.c)
TSTS = $(SRCS:.c=.tst)

unit-test: library unit-test-trail-whitespaces unit-tsts unit-test-debug

library:
	$(MAKE) -C ../lib all

unit-test-trail-whitespaces:
	@./test-trailing-whitespaces.sh

unit-test-debug:
	@if [ ! -z "$(DEBUG)" ]; then \
		echo -e "\e[0;31mDon't commit with DEBUG\e[0m"; exit 1; fi

unit-tsts: $(TSTS)
	@for tst in $(TSTS); do \
		echo -e "\e[0;32mTesting $$tst\e[0m"; \
		if ! ./$$tst; then \
			echo -e "\e[0;31mFailed\e[0m"; exit 1;\
		fi; done

%.tst: %.c
	$(CC) -o $@ $< ../lib/lib$(PACKAGE).a $(CFLAGS) $(I)

clean:
	$(RM) $(TSTS)
