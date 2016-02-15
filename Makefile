all:
	cd Common; make
	cd Linear; make
	cd Nanoparticles; make

clean:
	rm -rf obj *.o *~ core
	cd Common; make clean
	cd Linear; make clean
	cd Nanoparticles; make clean
