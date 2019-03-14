
## Variables ###################################################################
ec_fp1 = ../chart-een1/sel.ec.dat
ec_fp2 = ../chart-een2/sel.ec.dat
dval = 0.0001

thresh = 0.003
# tpl_thresh = 0.008
line_values  = 0.015 0.008 0.003
line_colors  = red   blue  green
line_thicks  = 4.0   4.0   2.0
line_weights =   5     3     1
other_opts = --ratio 0.3 --direction TB
		# --forced-output-nodes 350 380 \

################################################################################
targets = inc-fill dec-fill inc dec both
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

## Targets #####################################################################
outdata/inc.dat: $(ec_fp1) $(ec_fp2)
	mkdir -p outdata
	python $(toolpath)/cal-diffec.py $^ -o $@ > $*.log

outdata/dec.dat: $(ec_fp2) $(ec_fp1)
	mkdir -p outdata
	python $(toolpath)/cal-diffec.py $^ -o $@ > $*.log

outdata/inc-fill.dat: $(ec_fp1) $(ec_fp2) config.mk
	mkdir -p outdata
	python $(toolpath)/cal-diffec.py \
		--fill-deficit-pairs --deficit-value $(dval) \
		$^ -o $@ > $*.log

outdata/dec-fill.dat: $(ec_fp2) $(ec_fp1) config.mk
	mkdir -p outdata
	python $(toolpath)/cal-diffec.py \
		--fill-deficit-pairs --deficit-value $(dval) \
		$^ -o $@ > $*.log

other_opts += $(shell [ -n "$(thresh)" ] && echo " -r $(thresh)")
other_opts += $(shell [ -n "$(tpl_thresh)" ] && echo " -R $(tpl_thresh)")

outdata/inc-fill.ec.pdf: outdata/inc-fill.dat $(conf_files)
	$(graph_een) $(ec_fp1) $(ec_fp2) $< -f $@ $(other_opts) \
		--line-values $(line_values) --line-colors  $(line_colors) \
		--line-thicks $(line_thicks) --line-weights $(line_weights) \
		--with-one-letter \
		-n node.cfg -c cluster.cfg \
		2>&1 | grep  -v "CoreText performance note"

outdata/dec-fill.ec.pdf: outdata/dec-fill.dat $(conf_files)
	$(graph_een) $(ec_fp1) $(ec_fp2) $< -f $@ $(other_opts) \
		--line-values $(line_values) --line-colors  $(line_colors) \
		--line-thicks $(line_thicks) --line-weights $(line_weights) \
		--with-one-letter \
		-n node.cfg -c cluster.cfg \
		2>&1 | grep  -v "CoreText performance note"

outdata/inc.ec.pdf: outdata/inc.dat  $(conf_files)
	$(graph_een) $(ec_fp1) $(ec_fp2) $< -f $@ $(other_opts) \
		--line-values $(line_values) --line-colors  $(line_colors) \
		--line-thicks $(line_thicks) --line-weights $(line_weights) \
		--with-one-letter \
		-n node.cfg -c cluster.cfg \
		2>&1 | grep  -v "CoreText performance note"

outdata/dec.ec.pdf: outdata/dec.dat $(conf_files)
	$(graph_een) $(ec_fp1) $(ec_fp2) $< -f $@ $(other_opts) \
		--line-values $(line_values) --line-colors  $(line_colors) \
		--line-thicks $(line_thicks) --line-weights $(line_weights) \
		--with-one-letter \
		-n node.cfg -c cluster.cfg \
		2>&1 | grep  -v "CoreText performance note"

outdata/both.ec.pdf: outdata/inc.dat $(conf_files)
	$(graph_een) $(ec_fp1) $(ec_fp2) $< -f $@ $(other_opts) \
		--line-values $(line_values) --line-colors  $(line_colors) \
		--line-thicks $(line_thicks) --line-weights $(line_weights) \
		--with-one-letter --show-negative-values \
		-n node.cfg -c cluster.cfg \
		2>&1 | grep  -v "CoreText performance note"
