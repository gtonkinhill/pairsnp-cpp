CXX = g++
CXXFLAGS = -Wall -Wextra -Ofast -std=c++11 
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
	./$(EXE) -qb $(TESTDIR)/singleton.aln | diff -bB - $(TESTDIR)/singleton.res
	./$(EXE) -qb $(TESTDIR)/good.aln | diff -bB - $(TESTDIR)/good.res
	./$(EXE) -qb $(TESTDIR)/gzip.aln.gz | diff -bB -  $(TESTDIR)/gzip.res
	./$(EXE) -qb -k $(TESTDIR)/lowercase.aln | diff -bB - $(TESTDIR)/lowercase-k.res
	./$(EXE) -qb    $(TESTDIR)/lowercase.aln | diff -bB - $(TESTDIR)/lowercase.res
	./$(EXE) -qb -c -q $(TESTDIR)/good.aln | diff -bB - $(TESTDIR)/good-c.res
	./$(EXE) -qb -a $(TESTDIR)/ambig.aln | diff -bB - $(TESTDIR)/ambig-a.res
	./$(EXE) -qbcm $(TESTDIR)/good.aln | diff -bB - $(TESTDIR)/good-c-m.res
	./$(EXE) -qb    $(TESTDIR)/ambig.aln | diff -bB - $(TESTDIR)/ambig.res

format:
	clang-format -i src/main.cpp