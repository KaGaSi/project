all:
	cd Common; make
	cd Linear; make

clean:
	rm -rf obj *.o *~ core
	cd Common; make clean
	cd Linear; make clean
