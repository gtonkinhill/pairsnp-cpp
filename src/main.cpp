#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <fcntl.h>
#include <vector>
#include <iterator>
#include <string>
#include "kseq.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <time.h>
#include <ctime>

#include "roaring.hh"
#include "roaring.c"
using namespace roaring;


#define VERSION "0.3.0"
#define EXENAME "pairsnp"

void show_help(int retcode)
{
  FILE* out = (retcode == EXIT_SUCCESS ? stdout : stderr);
  fprintf(out, "SYNOPSIS\n  Pairwise SNP distance matrices using fast matrix algerbra libraries\n");
  fprintf(out, "USAGE\n  %s [options] alignment.fasta[.gz] > matrix.csv\n", EXENAME);
  fprintf(out, "OPTIONS\n");
  fprintf(out, "  -h\tShow this help\n");
  fprintf(out, "  -v\tPrint version and exit\n");
  fprintf(out, "  -s\tOutput in sparse matrix form (i,j,distance).\n");
  fprintf(out, "  -d\tDistance threshold for sparse output. Only distances <= d will be returned.\n");
  fprintf(out, "  -k\tWill on return the k nearest neighbours for each sample in sparse output.\n");
  fprintf(out, "  -c\tOutput CSV instead of TSV\n");
  fprintf(out, "  -n\tCount comparisons with Ns (off by default)\n");
  fprintf(out, "  -t\tNumber of threads to use (default=1)\n");
  exit(retcode);
}


int main(int argc, char *argv[])
{

  // parse command line parameters
  int opt, quiet=0, csv=0, sparse=0, dist=-1, knn=-1;
  size_t nthreads=1;
  while ((opt = getopt(argc, argv, "hqsncvd:t:k:")) != -1) {
    switch (opt) {
      case 'h': show_help(EXIT_SUCCESS); break;
      case 'q': quiet=1; break;
      case 'c': csv=1; break;
      case 's': sparse=1; break;
      case 'd': dist=atoi(optarg); break;
      case 'k': knn=atoi(optarg); break;
      case 't': nthreads=atoi(optarg); break;
      case 'v': printf("%s %s\n", EXENAME, VERSION); exit(EXIT_SUCCESS);
      default : show_help(EXIT_FAILURE);
    }
  }

  char sep = csv ? ',' : '\t';

  // require a filename argument
  if (optind >= argc) {
    show_help(EXIT_FAILURE);
    return 0;
  }
  const char* fasta = argv[optind];

  // say hello
  if (!quiet) fprintf(stderr, "This is %s %s\n", EXENAME, VERSION);

  // open filename via libz
  int fp = open(fasta, O_RDONLY);
  if (! fp) {
    fprintf(stderr, "ERROR: Could not open filename '%s'\n", fasta);
    exit(EXIT_FAILURE);
  }

  //Initially run through fasta to get SNPS not matching reference
  size_t n_seqs = 0;
  int l=0;
  kseq seq;
  FunctorRead r;
  kstream<int, FunctorRead> ks(fp, r);

  // Use first sequence to set size of arrays
  l = ks.read(seq);
  size_t seq_length = seq.seq.length();

  // initialise roaring bitmaps
  std::vector<std::string> seq_names;
  std::vector<Roaring> A_snps;
  std::vector<Roaring> C_snps;
  std::vector<Roaring> G_snps;
  std::vector<Roaring> T_snps;

  do {
    seq_names.push_back(seq.name);
    Roaring As;
    Roaring Cs;
    Roaring Gs;
    Roaring Ts;

    for(size_t j=0; j<seq_length; j++){

      switch(std::toupper(seq.seq[j])){
        case 'A': As.add(j); break;
        case 'C': Cs.add(j); break;
        case 'G': Gs.add(j); break;
        case 'T': Ts.add(j); break;

        // M = A or C
        case 'M': As.add(j); 
          Cs.add(j);
          break;

        // R = A or G
        case 'R': As.add(j); 
          Gs.add(j);
          break;

        // W = A or T
        case 'W': As.add(j); 
          Ts.add(j);
          break;

        // S = C or G
        case 'S': Cs.add(j); 
          Gs.add(j);
          break;

        // Y = C or T
        case 'Y': Cs.add(j); 
          Ts.add(j);
          break;

        // K = G or T
        case 'K': Gs.add(j); 
          Ts.add(j);
          break;

        // V = A,C or G
        case 'V': As.add(j); 
          Cs.add(j); 
          Gs.add(j);
          break;

        // H = A,C or T
        case 'H': As.add(j); 
          Cs.add(j); 
          Ts.add(j);
          break;
        
        // D = A,G or T
        case 'D': As.add(j);
          Gs.add(j); 
          Ts.add(j);
          break;

        // B = C,G or T
        case 'B': Cs.add(j);
          Gs.add(j); 
          Ts.add(j);
          break;

        // N = A,C,G or T
        default:  
          As.add(j);
          Cs.add(j);
          Gs.add(j); 
          Ts.add(j);
          break;

      }
    }
    As.runOptimize();
    A_snps.push_back(As);
    Cs.runOptimize();
    C_snps.push_back(Cs);
    Gs.runOptimize();
    G_snps.push_back(Gs);
    Ts.runOptimize();
    T_snps.push_back(Ts);

    n_seqs++;
  } while((l = ks.read(seq)) >= 0);
  close(fp);

  if((knn!=-1) && (knn>=int(n_seqs))){
      fprintf(stderr, "kNN > number of samples. Running in dense mode!\n" );
      knn = -1;
  }

  #pragma omp parallel for ordered shared(A_snps, C_snps \
    , G_snps, T_snps, seq_length \
    , n_seqs, seq_names, dist, sep, sparse \
    , knn) default(none) schedule(static,1) num_threads(nthreads)
  for (size_t i = 0; i < n_seqs; i++) {

    std::vector<int> comp_snps(n_seqs);
    
    for(size_t j=0; j<n_seqs; j++){

      Roaring res = A_snps[i] & A_snps[j];
      Roaring intersect = C_snps[i] & C_snps[j];
      res = res | intersect;
      intersect = G_snps[i] & G_snps[j];
      res = res | intersect;
      intersect = T_snps[i] & T_snps[j];
      res = res | intersect;

      comp_snps[j] = seq_length - res.cardinality();

    }

    // if using knn find the distance needed
    int start;
    if (knn>=0){
      std::vector<int> s_comp = comp_snps;
      std::sort(s_comp.begin(), s_comp.end());
      dist = s_comp[knn+1];
      start=0;
    } else {
      start = i+1;
    }

    // Output the distance matrix to stdout
    #pragma omp ordered //#pragma omp critical
    if (sparse){
      for (size_t j=start; j < n_seqs; j++) {
        if ((dist==-1) || (comp_snps[j] <= dist)){
          printf("%d%c%d%c%d\n", int(i), sep, int(j), sep, comp_snps[j]);
        }
      }
    } else {
      printf("%s", seq_names[i].c_str());
      for (size_t j=0; j < n_seqs; j++) {
        printf("%c%d", sep, comp_snps[j]);
      }
      printf("\n");
    }
  }

  return 0;

  }


