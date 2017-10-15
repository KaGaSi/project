all:
	doxygen doxyconfig;
	cd latex/; pdflatex refman.tex
	cd latex/; pdflatex refman.tex
	cp latex/refman.pdf .
	if [ -d "doc" ]; then rm -r doc/*; else mkdir doc; fi
	mv latex html DoxyErrors.txt doc
