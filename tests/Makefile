include ../make.inc

unit-test: unit-test-trail-whitespaces

unit-test-trail-whitespaces:
	@[ -z "$$(grep -r '\s\+$$' ../*)" ] || exit 1
