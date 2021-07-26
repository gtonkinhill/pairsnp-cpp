CXX = g++

ifdef ompoff
	CXXFLAGS = -Wall -Wextra -Ofast -std=c++11 
else
	CXXFLAGS = -Wall -Wextra -Ofast -std=c++11 -fopenmp
endif
 
LIBS = -lz -lm -I libs/CRoaringUnityBuild/ -I libs/

EXE = pairsnp
PREFIX = /usr/local
TESTDIR = test

.PHONY: all check clean format 
.DEFAULT: all

all: $(EXE)


$(EXE): src/main.cpp
	$(CXX) $(CXXFLAGS) -o $(EXE) $< $(OBJECTS) $(LIBS) 

install: $(EXE)
	install -v -t $(PREFIX)/bin $(EXE)

clean:
	$(RM) *~ *.o $(EXE)

check: $(EXE)
	./$(EXE) -v
	./$(EXE) /dev/null || true
	./$(EXE) $(TESTDIR)/empty.aln || true
	./$(EXE) $(TESTDIR)/singleton.aln | diff -bB - $(TESTDIR)/singleton.out
	./$(EXE) $(TESTDIR)/good.aln | diff -bB - $(TESTDIR)/good.out
	./$(EXE) $(TESTDIR)/good.aln.gz | diff -bB -  $(TESTDIR)/good.out
	./$(EXE) -s $(TESTDIR)/good.aln | diff -bB - $(TESTDIR)/sparse.out
	./$(EXE) -s -i $(TESTDIR)/good.aln | diff -bB - $(TESTDIR)/sparse_index.out
	./$(EXE) -s -d 2 $(TESTDIR)/good.aln | diff -bB - $(TESTDIR)/filter.out
	./$(EXE) $(TESTDIR)/lowercase.aln | diff -bB - $(TESTDIR)/lowercase.out
	./$(EXE) $(TESTDIR)/ambig.aln | diff -bB - $(TESTDIR)/ambig.out
	./$(EXE) $(TESTDIR)/bad.aln 2>&1 >/dev/null | tail -n1 | diff -bB - $(TESTDIR)/bad.out


format:
	clang-format -i src/main.cpp