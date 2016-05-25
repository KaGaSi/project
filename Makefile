all:
	cd Common; make
	cd Linear; make
	cd OneTime; make

clean:
	rm -rf bin/* obj *.o *~ core
	cd Common; make clean
	cd Linear; make clean
	cd OneTime; make clean

doc:
	doxygen doxyconfig;
	cd latex/; pdflatex refman.tex
	cd latex/; pdflatex refman.tex
	cp latex/refman.pdf .
