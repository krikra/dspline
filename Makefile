default:
	(cd storage;make)
	(cd solver;make)
	(cd dspline;make)
	(cd ippe;make)
	(cd fort;make)

clean:
	(cd storage;make clean)
	(cd solver;make clean)
	(cd dspline;make clean)
	(cd ippe;make clean)
	(cd fort;make clean)
