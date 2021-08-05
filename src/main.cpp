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

#include <dynamic_bitset.hpp>


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
  std::vector<boost::dynamic_bitset<>> A_snps;
  std::vector<boost::dynamic_bitset<>> C_snps;
  std::vector<boost::dynamic_bitset<>> G_snps;
  std::vector<boost::dynamic_bitset<>> T_snps;
  std::string consensus;

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
    boost::dynamic_bitset<> As(seq_length);
    boost::dynamic_bitset<> Cs(seq_length);
    boost::dynamic_bitset<> Gs(seq_length);
    boost::dynamic_bitset<> Ts(seq_length);

    for(size_t j=0; j<seq_length; j++){

      seq->seq.s[j] = std::toupper(seq->seq.s[j]);

      switch(seq->seq.s[j]){
        case 'A': As[j] = 1; break;
        case 'C': Cs[j] = 1; break;
        case 'G': Gs[j] = 1; break;
        case 'T': Ts[j] = 1; break;

        // M = A or C
        case 'M': As[j] = 1; 
          Cs[j] = 1;
          break;

        // R = A or G
        case 'R': As[j] = 1; 
          Gs[j] = 1;
          break;

        // W = A or T
        case 'W': As[j] = 1; 
          Ts[j] = 1;
          break;

        // S = C or G
        case 'S': Cs[j] = 1; 
          Gs[j] = 1;
          break;

        // Y = C or T
        case 'Y': Cs[j] = 1; 
          Ts[j] = 1;
          break;

        // K = G or T
        case 'K': Gs[j] = 1; 
          Ts[j] = 1;
          break;

        // V = A,C or G
        case 'V': As[j] = 1; 
          Cs[j] = 1; 
          Gs[j] = 1;
          break;

        // H = A,C or T
        case 'H': As[j] = 1; 
          Cs[j] = 1; 
          Ts[j] = 1;
          break;
        
        // D = A,G or T
        case 'D': As[j] = 1;
          Gs[j] = 1; 
          Ts[j] = 1;
          break;

        // B = C,G or T
        case 'B': Cs[j] = 1;
          Gs[j] = 1; 
          Ts[j] = 1;
          break;

        // N = A,C,G or T
        default:  
          As[j] = 1;
          Cs[j] = 1;
          Gs[j] = 1; 
          Ts[j] = 1;
          break;

      }
    }
    // As.runOptimize();
    A_snps.push_back(As);
    // Cs.runOptimize();
    C_snps.push_back(Cs);
    // Gs.runOptimize();
    G_snps.push_back(Gs);
    // Ts.runOptimize();
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
    boost::dynamic_bitset<> res(seq_length);

    size_t start;
    if (sparse && (knn<0)){
      start = i+1;
    } else {
      start = 0;
    }
    
    for(size_t j=start; j<n_seqs; j++){

      res = A_snps[i] & A_snps[j];
      res |= C_snps[i] & C_snps[j];
      res |= G_snps[i] & G_snps[j];
      res |= T_snps[i] & T_snps[j];

      comp_snps[j] = seq_length - res.count();

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
            std::printf("%d%c%d%c%d\n", i, sep, j, sep, comp_snps[j]);
          } else {
            std::printf("%s%c%s%c%d\n", seq_names[i].c_str(), sep, seq_names[j].c_str(), sep, comp_snps[j]);
          }
        }
      }
    } else {
      if (index) {
        printf("%d", i);
      } else {
        printf("%s", seq_names[i].c_str());
      }
      for (size_t j=0; j < n_seqs; j++) {
        printf("%c%d", sep, comp_snps[j]);
      }
      printf("\n");
    }
  }

  return 0;

  }


