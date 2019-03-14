
## Variables ###################################################################
pdb_fp = system-renum.pdb
ec_fp  = 
pml_fp = ec.pml
thresh = 0.003
resolution = 820
movie_style = 'y-rotation'
movie_options =

################################################################################
all: view
################################################################################

-include config.mk

# PATH to command
cmdpath = $(CURP_HOME)/script
toolpath = $(CURP_HOME)/tool
PYMOL = pymol -S

depend_files = $(pdb_fp) $(ec_fp) $(pml_fp)

all: $(targets)

clean:
	@rm -rf sel.ec.dat _movie *.png movie.gif

% : %.ec.pdf
	@echo "Creating $^ ..."

# standard
view: $(pdb_fp) $(ec_fp) $(pml_fp)
	$(PYMOL) $(pdb_fp) $(pml_fp) -r $(cmdpath)/view_ec.py -- $(ec_fp) -t $(thresh) -d

view-inc: $(pdb_fp) $(ec_fp) $(pml_fp)
	$(PYMOL) $(pdb_fp) $(pml_fp) -r $(cmdpath)/view_ec.py -- $(ec_fp) -t $(thresh)

view-dec: $(pdb_fp) $(ec_fp) $(pml_fp)
	$(PYMOL) $(pdb_fp) $(pml_fp) -r $(cmdpath)/view_ec.py -- $(ec_fp) -t $(thresh) -I -d

images: een-3D-inc een-3D-dec een-3D-both

een-3D-inc: $(pdb_fp) $(ec_fp) $(pml_fp)
	$(PYMOL) $(pdb_fp) $(pml_fp) \
    -d "hide labels; turn y, -60" \
	-r $(cmdpath)/view_ec.py -- $(ec_fp) -t $(threst) -r 820 -o $@-60.png -L
	$(PYMOL) $(pdb_fp) $(pml_fp) \
    -d "hide labels;" \
	-r $(cmdpath)/view_ec.py -- $(ec_fp) -t $(threst) -r 820 -o $@.png -L
	$(PYMOL) $(pdb_fp) $(pml_fp) \
	-d "hide labels; turn y, 60" \
	-r $(cmdpath)/view_ec.py -- $(ec_fp) -t $(threst) -r 820 -o $@+60.png -L

een-3D-dec: $(pdb_fp) $(ec_fp) $(pml_fp)
	$(PYMOL) $(pdb_fp) $(pml_fp) \
    -d "hide labels; turn y, -60" \
	-r $(cmdpath)/view_ec.py -- $(ec_fp) -t $(threst) -r 820 -o $@-60.png -d -I -L
	$(PYMOL) $(pdb_fp) $(pml_fp) \
    -d "hide labels;" \
	-r $(cmdpath)/view_ec.py -- $(ec_fp) -t $(threst) -r 820 -o $@.png -d -I -L
	$(PYMOL) $(pdb_fp) $(pml_fp) \
	-d "hide labels; turn y, 60" \
	-r $(cmdpath)/view_ec.py -- $(ec_fp) -t $(threst) -r 820 -o $@+60.png -d -I -L

een-3D-both: $(pdb_fp) $(ec_fp) $(pml_fp)
	$(PYMOL) $(pdb_fp) $(pml_fp) \
    -d "turn y, -60" \
	-r $(cmdpath)/view_ec.py -- $(ec_fp) -t $(threst) -r 820 -o $@-60.png -d
	$(PYMOL) $(pdb_fp) $(pml_fp) \
	-d "" \
	-r $(cmdpath)/view_ec.py -- $(ec_fp) -t $(threst) -r 820 -o $@.png -d
	$(PYMOL) $(pdb_fp) $(pml_fp) \
	-d "turn y, 60" \
	-r $(cmdpath)/view_ec.py -- $(ec_fp) -t $(threst) -r 820 -o $@+60.png -d

movie: $(depend_files)
	rm -f movie.gif
	mkdir -p _movie
	$(PYMOL) $(pdb_fp) $(pml_fp) \
	-r $(cmdpath)/view_ec.py -- $(ec_fp) -t $(threst) -d \
	-ms $(movie_style) $(movie_options)
	convert -delay 5.0 -loop 0 _movie/0???.png  movie.gif

