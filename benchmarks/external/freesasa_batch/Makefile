CXX      ?= c++
CXXFLAGS ?= -O3 -std=c++17
FREESASA ?= ../freesasa

freesasa_batch: freesasa_batch.cc
	$(CXX) $(CXXFLAGS) -I $(FREESASA)/src \
		-o $@ $< \
		$(FREESASA)/src/.libs/libfreesasa.a -lpthread

clean:
	rm -f freesasa_batch

.PHONY: clean
