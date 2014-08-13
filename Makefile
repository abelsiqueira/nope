include make.inc

all:
	$(MAKE) -C lib

clean:
	$(MAKE) -C lib clean
	$(MAKE) -C tests clean

test:
	$(MAKE) -C tests

install:
	$(CP) lib/lib$(PACKAGE).a $(PREFIX)/lib/
	$(CP) include/* $(PREFIX)/include/
