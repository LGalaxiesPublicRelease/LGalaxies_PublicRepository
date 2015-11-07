#!/bin/sh

# beautifies code with GNU-indent

/afs/ipp/home/v/vrs/bin/indent_psi/bin/indent -gnu -npsl -npcs -nbs -nsaf -nsai -nsaw -nprs -bad -bap -pmt -nut  -l110 *.c

# indent -gnu -npsl -npcs -nbs -nsaf -nsai -nsaw -nprs -bad -bap -pmt -nut  -l110 *.c

rm -f *.c~
