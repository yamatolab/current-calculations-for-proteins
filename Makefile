# Makefile for building the CURP program
#

MAKE = make -w
VENV_VERSION=16.4.3
export FCOMPILER = # '' or intelem, pass to sub make.

#######
CURP_HOME=$(PWD)
VIRTUAL_DIR=$(CURP_HOME)/curp-environ
VENV=virtualenv-$(VENV_VERSION)
ACTIVATE="$(VIRTUAL_DIR)/bin/activate"

PY_PACKAGES = numpy nose benchmarker setproctitle netCDF4 pygraphviz
PIP_DOWNLOAD_CACHE=$(VIRTUAL_DIR)/pylibs

all: serial success

environment: venv packages

serial: analysis curp

intel: venv netcdf-intel $(PY_PACKAGES) analysis-intel curp-intel success-intel
 
mpi: mpi4py

allclean:
	@rm -rf docs/_api
	@rm -rf curp-environ
	$(MAKE) clean -C script
	$(MAKE) clean -C src

clean:
	@rm -rf docs/_api
	$(MAKE) clean -C curp-environ
	$(MAKE) clean -C script
	$(MAKE) clean -C src
	$(MAKE) clean -C docs

venv: $(VIRTUAL_DIR)/$(VENV) $(VIRTUAL_DIR)/bin/activate

$(VIRTUAL_DIR)/$(VENV): 
	mkdir -p $(VIRTUAL_DIR)
	curl  -k -o $(VIRTUAL_DIR)/$(VENV).tar.gz \
		https://github.com/pypa/virtualenv/archive/1.11.6.tar.gz
	tar xzf $(VIRTUAL_DIR)/$(VENV).tar.gz -C $(VIRTUAL_DIR)

$(VIRTUAL_DIR)/bin/activate:
	mkdir -p $(VIRTUAL_DIR)/logs
	cp $(CURP_HOME)/lib/Makefile.environ $(VIRTUAL_DIR)/Makefile
	# create virtual environment
	python2.7 $(VIRTUAL_DIR)/$(VENV)/virtualenv.py $(VIRTUAL_DIR) \
		--no-site-packages
	@echo ""

netcdf:
	$(MAKE) -C curp-environ

netcdf-intel:
	$(MAKE) intel -C curp-environ

packages:
	@echo "Installing $@ (time-consuming)"
	. $(VIRTUAL_DIR)/bin/activate; \
		$(VIRTUAL_DIR)/bin/pip install -r requirements.txt >> $(VIRTUAL_DIR)/logs/$@.log 2>&1 || \
		pip install $(PIP_DOWNLOAD_CACHE)/$@* \
		>> $(VIRTUAL_DIR)/logs/$@.log 2>&1 && \
		exit


numpy:
	@echo "Installing $@ (may be time-consuming)..."
	. $(VIRTUAL_DIR)/bin/activate; \
		$(VIRTUAL_DIR)/bin/pip install $@==1.11.2 >> $(VIRTUAL_DIR)/logs/$@.log 2>&1 || \
		pip install $(PIP_DOWNLOAD_CACHE)/$@* \
		#pip install install $@==1.11.2 >> $(VIRTUAL_DIR)/logs/$@.log 2>&1 || 
		#pip install $(PIP_DOWNLOAD_CACHE)/$@* 
		>> $(VIRTUAL_DIR)/logs/$@.log 2>&1 && \
		exit

f2py:
	@echo making symbolic link: $(VIRTUAL_DIR)/bin/f2py-curp
	@if [ -f $(VIRTUAL_DIR)/bin/f2py2 ]; then\
		ln -sf $(VIRTUAL_DIR)/bin/f2py2 $(VIRTUAL_DIR)/bin/f2py-curp ;\
	elif [ -f $(VIRTUAL_DIR)/bin/f2py2.7 ]; then\
		ln -sf $(VIRTUAL_DIR)/bin/f2py2.7 $(VIRTUAL_DIR)/bin/f2py-curp ;\
	elif [ -f $(VIRTUAL_DIR)/bin/f2py-2.7 ]; then\
		ln -sf $(VIRTUAL_DIR)/bin/f2py-2.7 $(VIRTUAL_DIR)/bin/f2py-curp ;\
	elif [ -f $(VIRTUAL_DIR)/bin/f2py ]; then\
		ln -sf $(VIRTUAL_DIR)/bin/f2py $(VIRTUAL_DIR)/bin/f2py-curp ;\
	fi

	@echo "Processing dependencies for $@"
	@echo "Finished processing dependencies for $@"

mpi4py:
	@echo "Installing $@ (may be time-consuming)..."
	@. $(VIRTUAL_DIR)/bin/activate; \
		pip install $@==2.0.0 >> $(VIRTUAL_DIR)/logs/$@.log 2>&1 || \
		pip install $(PIP_DOWNLOAD_CACHE)/$@* \
		>> $(VIRTUAL_DIR)/logs/$@.log 2>&1 && \
		exit
	@echo "Processing dependencies for $@"
	@echo "Finished processing dependencies for $@"
	@echo

nose:
	@. $(VIRTUAL_DIR)/bin/activate; \
		pip install $@ || \
		pip install $(PIP_DOWNLOAD_CACHE)/$@* && \
		exit

benchmarker:
	. $(VIRTUAL_DIR)/bin/activate; \
		pip install $@ || \
		pip install $(PIP_DOWNLOAD_CACHE)/Benchmarker* && \
		exit

pygraphviz:
	. $(VIRTUAL_DIR)/bin/activate; \
		pip install $@ || \
		pip install $(PIP_DOWNLOAD_CACHE)/$@* && \
		exit

setproctitle:
	. $(VIRTUAL_DIR)/bin/activate; \
		pip install $@ || \
		pip install $(PIP_DOWNLOAD_CACHE)/$@* && \
		exit

export NETCDF4_DIR = $(VIRTUAL_DIR)
export HDF5_DIR    = $(VIRTUAL_DIR)
netCDF4:
	@echo "Installing $@ (may be time-consuming)..."
	@. $(VIRTUAL_DIR)/bin/activate; \
		pip install $@==1.2.4 >> $(VIRTUAL_DIR)/logs/$@.log 2>&1 || \
		pip install $(PIP_DOWNLOAD_CACHE)/$@* \
		>> $(VIRTUAL_DIR)/logs/$@.log 2>&1 && \
		exit
	@echo "Processing dependencies for $@"
	@echo "Finished processing dependencies for $@"
	@echo

epydoc:
	. $(VIRTUAL_DIR)/bin/activate; \
		pip install $@ || \
		pip install $(PIP_DOWNLOAD_CACHE)/$@* && \
		exit

sphinx:
	. $(VIRTUAL_DIR)/bin/activate; \
		pip install install $@ || \
		pip install $(PIP_DOWNLOAD_CACHE)/$@* && \
		exit

analysis:
	$(MAKE) -C script

curp:
	$(MAKE) -C src

analysis-intel:
	$(MAKE) intel -C script

curp-intel:
	$(MAKE) intel -C src

success:
	@echo ""
	@echo "    *************************************************"
	@echo "    **            Congratulations!                 **"
	@echo "    **     CURP has been successfully installed.   **"
	@echo "    *************************************************"

success-intel:
	@echo ""
	@echo "    ******************************************************"
	@echo "    **                Congratulations!                  **"
	@echo "    ** The CURP program has been compiled with intel    **"
	@echo "    ** fortran compiler and installed successfully.     **"
	@echo "    ******************************************************"

apihtml: epydoc
	. $(ACTIVATE) ;\
	epydoc --show-imports -o docs/_build/api --graph all src/*.py src/**/*.py ;\
	exit

apipdf: epydoc
	. $(ACTIVATE) ;\
	epydoc --pdf --show-imports -o docs/_build/api --graph all src/*.py \
	src/**/*.py ;\ exit

docs: sphinx
	. $(ACTIVATE) ; \
		$(MAKE) -C docs html ;\
		exit
