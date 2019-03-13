# Makefile for building the utility scripts in the CURP program
#

CURP_HOME = ../
# include $(CURP_HOME)/lib/Makefile.in

TARGETS = current twobody forcefield volume table \
		  parser/amber

##### Target definitions #####

install:
	@for subdir in $(TARGETS); do\
		$(MAKE) -C $$subdir || exit 1; \
	done

intel:
	@for subdir in $(TARGETS); do\
		$(MAKE) intel -C $$subdir || exit 1; \
	done

clean:
	@for subdir in $(TARGETS); do\
		$(MAKE) clean -C $$subdir; \
	done

