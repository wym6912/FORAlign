THREADS = 1

all:
	make -C foralign all -j$(THREADS)
	make -C swg all -j$(THREADS)
	make -C library_example all -j$(THREADS)
	-make -C wfa2-test all -j$(THREADS)
	make -C wfa2-test

clean:
	rm -rf bin/ foralign/*.o lib
	make -C wfa2-test/ clean

export
$(folders):
	$(MAKE) -C $@ all
