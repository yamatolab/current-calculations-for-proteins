
## Variables ###################################################################
ec_fp = ../all.ec.dat

fix_resnums =

thresh =  0.008
thresh_weak = 0.003

line_values  = 0.015 0.008 0.003
line_colors  = red   blue  green
line_thicks  = 4.0   4.0   2.0
line_weights =   5     3     1
other_opts = --ratio 0.5 --direction TB -I

################################################################################
targets = strong weak
################################################################################

-include config.mk

# PATH to command
graph_een = $(CURP_HOME)/bin/graph-een
toolpath = $(CURP_HOME)/tool

# configuration files
conf_files = config.mk node.cfg cluster.cfg

all: $(targets)

clean:
	@rm -rf outdata

%: outdata/%.ec.pdf;

# standard
outdata/sel.ec.dat: $(ec_fp) $(conf_files)
	mkdir -p outdata
	cat $(ec_fp) \
		| sed -e "s/CYX/CYS/g" \
		| $(toolpath)/renum_residue.py $(fix_resnums) \
		| $(toolpath)/sel_noneighbor.py WAT $(ligand) \
		| $(toolpath)/sel_nocap.py \
		> $@
		# | ./sel_thresh_dist.py dist.dat 6 3 \
		# | $(cmdpath)/sel_noasa.py asa.dat \

outdata/strong.ec.pdf: other_opts += $(shell [ -n "$(thresh)" ] && echo "-r $(thresh)")
outdata/strong.ec.pdf: outdata/sel.ec.dat $(conf_files)
	@echo "Creating $^ ..."
	$(graph_een) $< \
		-f $@ $(other_opts) \
		--line-values $(line_values) --line-colors  $(line_colors) \
		--line-thicks $(line_thicks) --line-weights $(line_weights) \
		--with-one-letter \
		-n node.cfg \
		-c cluster.cfg \
		2>&1 | grep  -v "CoreText performance note"

outdata/weak.ec.pdf: other_opts += $(shell [ -n "$(thresh_weak)" ] && echo "-r $(thresh_weak)")
outdata/weak.ec.pdf: outdata/sel.ec.dat $(conf_files)
	@echo "Creating $^ ..."
	$(graph_een) $< \
		-f $@ $(other_opts) \
		--line-values $(line_values) --line-colors  $(line_colors) \
		--line-thicks $(line_thicks) --line-weights $(line_weights) \
		--with-one-letter \
		-n node.cfg \
		-c cluster.cfg \
		2>&1 | grep  -v "CoreText performance note"
