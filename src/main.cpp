#include <stdio.h>
#include <getopt.h>
#include <zlib.h>
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
KSEQ_INIT(gzFile, gzread)


#define VERSION "0.3.0"
#define EXENAME "pairsnp"

void show_help(int retcode)
{
  FILE* out = (retcode == EXIT_SUCCESS ? stdout : stderr);
  fprintf(out, "SYNOPSIS\n  Fast pairwise SNP distance matrices\n");
  fprintf(out, "USAGE\n  %s [options] alignment.fasta[.gz] > matrix.csv\n", EXENAME);
  fprintf(out, "OPTIONS\n");
  fprintf(out, "  -h\tShow this help\n");
  fprintf(out, "  -v\tPrint version and exit\n");
  fprintf(out, "  -s\tOutput in sparse matrix form (i,j,distance).\n");
  fprintf(out, "  -d\tDistance threshold for sparse output. Only distances <= d will be returned.\n");
  fprintf(out, "  -k\tWill on return the k nearest neighbours for each sample in sparse output.\n");
  fprintf(out, "  -c\tOutput CSV instead of TSV\n");
  fprintf(out, "  -i\tOutput sequence index inplace of header (useful for downstream processing)\n");
  fprintf(out, "  -t\tNumber of threads to use (default=1)\n");
  exit(retcode);
}


int main(int argc, char *argv[])
{

  // parse command line parameters
  int opt, quiet=0, csv=0, sparse=0, dist=-1, knn=-1, index=0;
  size_t nthreads=1;
  while ((opt = getopt(argc, argv, "hqsncvid:t:k:")) != -1) {
    switch (opt) {
      case 'h': show_help(EXIT_SUCCESS); break;
      case 'q': quiet=1; break;
      case 'c': csv=1; break;
      case 's': sparse=1; break;
      case 'i': index=1; break;
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

  // open filename and initialise kseq
  int l;
  gzFile fp = gzopen(fasta, "r");
  kseq_t * seq = kseq_init(fp);

  //Initially run through fasta to get SNPS not matching reference
  size_t n_seqs = 0;
  size_t seq_length;

  // initialise roaring bitmaps
  std::vector<std::string> seq_names;
  std::vector<Roaring> A_snps;
  std::vector<Roaring> C_snps;
  std::vector<Roaring> G_snps;
  std::vector<Roaring> T_snps;

  while(true) {
    l = kseq_read(seq);

    if (l == -1)  // end of file
        break;
    if (l == -2) {
        std::cerr << "Error: incorrect FASTA format" << seq->name.s << "\n";
        return 1;
    }
    if (l == -3) {
        std::cerr << "Error reading " << fasta << "\n";
        return 1;
    }

    // check sequence length
    if ((n_seqs>0) && (seq->seq.l != seq_length)){
      std::cerr << "Error: incorrect FASTA format, variable sequence lengths! " << seq->name.s << "\n";
      return 1;
    }
    seq_length = seq->seq.l;

    seq_names.push_back(seq->name.s);
    Roaring As;
    Roaring Cs;
    Roaring Gs;
    Roaring Ts;

    for(size_t j=0; j<seq_length; j++){

      switch(std::toupper(seq->seq.s[j])){
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
  }
  kseq_destroy(seq);
  gzclose(fp);

  if((knn!=-1) && (knn>=int(n_seqs))){
      fprintf(stderr, "kNN > number of samples. Running in dense mode!\n" );
      knn = -1;
  }

  //If sparse output print sequence names in the header
  if (index){
    std::cout << '#' << sep;
    for (size_t j=0; j < n_seqs; j++) {
      std::cout << seq_names[j] << sep;
    }
    std::cout << std::endl;
  }

  #pragma omp parallel for ordered shared(A_snps, C_snps \
    , G_snps, T_snps, seq_length \
    , n_seqs, seq_names, dist, sep, sparse \
    , knn, index) default(none) schedule(static,1) num_threads(nthreads)
  for (size_t i = 0; i < n_seqs; i++) {

    std::vector<int> comp_snps(n_seqs);

    size_t start;
    if (sparse && (knn<0)){
      start = i+1;
    } else {
      start = 0;
    }
    
    for(size_t j=start; j<n_seqs; j++){

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
          if (index) {
            std::cout << i;
          } else {
            std::cout << seq_names[i];
          }
          std::cout << sep;
          if (index) {
            std::cout << j;
          } else {
            std::cout << seq_names[j];
          }
          std::cout << sep << comp_snps[j] << std::endl;
        }
      }
    } else {
      if (index) {
        std::cout << i;
      } else {
        std::cout << seq_names[i];
      }
      for (size_t j=0; j < n_seqs; j++) {
        std::cout << sep << comp_snps[j];
      }
      std::cout << std::endl;
    }
  }

  return 0;

  }


