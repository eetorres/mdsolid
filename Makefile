# Edmanuel Torres 21/10/2004

include makeinclude

DIRS =	src bin

all:   makeinclude
	@for dir in $(DIRS); do\
	    echo "==== Compiling in $$dir ====";\
	    (cd $$dir; make || break);\
	done
	
clean: makeinclude
	@for dir in $(DIRS); do\
	    echo "==== Cleaning in $$dir ====";\
	    (cd $$dir; make clean || break);\
	done
	