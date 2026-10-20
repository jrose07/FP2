all: build/tex/vXXX.pdf 

var_plots = $(addprefix content/plots/, $(addsuffix .pdf, plot ))
MODE = None

# hier Python-Skripte:
$(var_plots) &: $(wildcard *.py) | plots 
	mamba run -n toolbox sh -c 'for py in *.py; do python $$py; done'

# hier weitere Abhängigkeiten für build/vXXX.pdf deklarieren:

build/tex/vXXX.pdf: FORCE $(var_plots) | tex
# to find header and bib files in the main directory
ifeq ($(MODE), write)
	TEXINPUTS=..: \
	BIBINPUTS=..: \
	max_print_line=1048576 \
	latexmk \
	  -pvc \
	  --lualatex \
	  --output-directory=build/tex \
	  --interaction=nonstopmode \
	  --halt-on-error \
	  -f \
	vXXX.tex
else
	TEXINPUTS=..: \
	BIBINPUTS=..: \
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build/tex \
	  --interaction=nonstopmode \
	  --halt-on-error \
	  -f \
	vXXX.tex
endif

build:
	mkdir -p build

plot: plots
	mamba run -n toolbox sh -c 'for py in *.py; do python $$py; done'

plots: build
	mkdir -p content/plots
	mkdir -p content/tables
	touch build/ergebnis.txt

tex:
	mkdir -p build/tex

clean:
	rm -rf build
	rm -rf content/plots
	rm -rf content/tables

FORCE:

.PHONY: all clean
