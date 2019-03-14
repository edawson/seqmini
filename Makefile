CXX?=g++
CXXFLAGS:= -O3 -pipe -fPIC -march=native -mtune=native -std=c++11 -fopenmp -g -ggdb -Wa,-q
LD_INC_FLAGS:= -I. -Imkmh -Igfakluge/src -ItinyFA -Imkmh/murmur3 -Imkmh/xxHash -ItinyFA/pliib -Impakt
LD_LIB_FLAGS:= -Lmkmh/murmur3 -lmurmur3
seqmini: src/main.cpp src/seqmini.hpp tinyFA/tinyfa.hpp gfakluge/src/gfakluge.hpp mkmh/mkmh.hpp mkmh/murmur3/libmurmur3.a
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

mkmh/murmur3/libmurmur3.a:
	cd mkmh && cd murmur3 && make

clean:
	$(RM) seqmini

.PHONY: clean
