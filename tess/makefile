#----------------------------------------------------------------------------
#
# top-level makefile
#
# Tom Peterka
# Argonne National Laboratory
# 9700 S. Cass Ave.
# Argonne, IL 60439
# tpeterka@mcs.anl.gov
#
# (C) 2012 by Argonne National Laboratory.
# See COPYRIGHT in top-level directory.
#
#----------------------------------------------------------------------------

default: all

all clean: FORCE
	cd lib && $(MAKE) $@
	cd driver && $(MAKE) $@
	cd postprocessor && $(MAKE) $@

FORCE:
	cd lib && $(MAKE) clean
	cd driver && $(MAKE) clean
	cd postprocessor && $(MAKE) clean
