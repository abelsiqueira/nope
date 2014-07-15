include make.inc

all:
	$(MAKE) -C lib

clean:
	$(MAKE) clean -C lib

test:
	$(MAKE) -C tests
