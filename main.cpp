#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <fcntl.h>
#include <vector>
#include "kseq.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <time.h>
#include <ctime>

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
  fprintf(out, "  -b\tBlank top left corner cell instead of '%s %s'\n", EXENAME, VERSION);
  exit(retcode);
}


int main(int argc, char *argv[])
{

  // parse command line parameters
  int opt, quiet=0, csv=0, corner=1, allchars=0, keepcase=0, sparse=0;
  int nthreads=1, dist=-1, knn=-1;
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
  int n_seqs = 0;
  int l=0;
  kseq seq;
  FunctorRead r;
  kstream<int, FunctorRead> ks(fp, r);

  // Use first sequence to set size of arrays
  l = ks.read(seq);
  int seq_length = seq.seq.length();
  int max_n_snps = 0, nsnps;
  std::string consensus = seq.seq;
  
  // initialise vectors to store SNPs and sequence names
  std::vector<std::string> seq_names;
  std::vector<std::vector<int> > A_snps;
  std::vector<std::vector<int> > C_snps;
  std::vector<std::vector<int> > G_snps;
  std::vector<std::vector<int> > T_snps;  

  do {
    seq_names.push_back(seq.name);
    std::vector<int> As;
    std::vector<int> Cs;
    std::vector<int> Gs;
    std::vector<int> Ts;

    nsnps=0;

    for(int j=0; j<seq_length; j++){

      //skip SNPs that match the reference
      if (seq.seq[j]==consensus[j]) continue;

      nsnps+=1; 
      
      switch(std::toupper(seq.seq[j])){
        case 'A': As.push_back(j); break;
        case 'C': Cs.push_back(j); break;
        case 'G': Gs.push_back(j); break;
        case 'T': Ts.push_back(j); break;

        // M = A or C
        case 'M': As.push_back(j); 
          Cs.push_back(j);
          break;

        // R = A or G
        case 'R': As.push_back(j); 
          Gs.push_back(j);
          break;

        // W = A or T
        case 'W': As.push_back(j); 
          Ts.push_back(j);
          break;

        // S = C or G
        case 'S': Cs.push_back(j); 
          Gs.push_back(j);
          break;

        // Y = C or T
        case 'Y': Cs.push_back(j); 
          Ts.push_back(j);
          break;

        // K = G or T
        case 'K': Gs.push_back(j); 
          Ts.push_back(j);
          break;

        // V = A,C or G
        case 'V': As.push_back(j); 
          Cs.push_back(j); 
          Gs.push_back(j);
          break;

        // H = A,C or T
        case 'H': As.push_back(j); 
          Cs.push_back(j); 
          Ts.push_back(j);
          break;
        
        // D = A,G or T
        case 'D': As.push_back(j);
          Gs.push_back(j); 
          Ts.push_back(j);
          break;

        // B = C,G or T
        case 'B': Cs.push_back(j);
          Gs.push_back(j); 
          Ts.push_back(j);
          break;

        // N = A,C,G or T
        default:  
          As.push_back(j);
          Cs.push_back(j);
          Gs.push_back(j); 
          Ts.push_back(j);
          break;

      }
    }

    if (nsnps>max_n_snps) max_n_snps=nsnps;

    A_snps.push_back(As);
    C_snps.push_back(Cs);
    G_snps.push_back(Gs);
    T_snps.push_back(Ts);

    n_seqs++;
  } while((l = ks.read(seq)) >= 0);
  close(fp);

  
  if((knn!=-1) && (knn>=n_seqs)){
      fprintf(stderr, "kNN > number of samples. Running in dense mode!\n" );
      knn = -1;
  }

  // fprintf(stderr, "  SizA: %d\n", A_snps.size());
  // fprintf(stderr, "  SizC: %d\n", C_snps.size());
  // fprintf(stderr, "  SizG: %d\n", G_snps.size());
  // fprintf(stderr, "  Size: %d\n", T_snps.size());
  // fprintf(stderr, "  seq_names: %d\n", seq_names.size());
  // fprintf(stderr, "  n_seqs: %d\n", n_seqs);

  #pragma omp parallel for ordered shared(A_snps, C_snps \
    , G_snps, T_snps \
    , n_seqs, max_n_snps \
    , seq_names, dist, sep, sparse \
    , knn) default(none) schedule(static,1) num_threads(nthreads)
  for (size_t i = 0; i < n_seqs; i++) {

    std::vector<int> comp_snps(n_seqs);
    std::vector<int> inter;
    inter.reserve(2*max_n_snps);
    std::vector<int> merged;
    merged.reserve(2*max_n_snps);
    
    for(int j=0; j<n_seqs; j++){

      inter.clear();
      merged.clear();
      int len=0;

      std::set_intersection(A_snps[i].begin(), A_snps[i].end(), 
        A_snps[j].begin(), A_snps[j].end(), std::back_inserter(inter) );
      
      merged.insert( merged.end(), A_snps[i].begin(), A_snps[i].end());
      len+=A_snps[i].size();
      merged.insert( merged.end(), A_snps[j].begin(), A_snps[j].end());
      std::inplace_merge (merged.begin(),merged.begin()+len,merged.end());
      len+=A_snps[j].size();

      std::set_intersection(C_snps[i].begin(), C_snps[i].end(), 
        C_snps[j].begin(), C_snps[j].end(), std::back_inserter(inter) );

      merged.insert( merged.end(), C_snps[i].begin(), C_snps[i].end());
      std::inplace_merge (merged.begin(),merged.begin()+len,merged.end());
      len+=C_snps[i].size();
      merged.insert( merged.end(), C_snps[j].begin(), C_snps[j].end());
      std::inplace_merge (merged.begin(),merged.begin()+len,merged.end());
      len+=C_snps[j].size();

      std::set_intersection(G_snps[i].begin(), G_snps[i].end(), 
        G_snps[j].begin(), G_snps[j].end(), std::back_inserter(inter));
      
      merged.insert( merged.end(), G_snps[i].begin(), G_snps[i].end());
      std::inplace_merge (merged.begin(),merged.begin()+len,merged.end());
      len+=G_snps[i].size();
      merged.insert( merged.end(), G_snps[j].begin(), G_snps[j].end());
      std::inplace_merge (merged.begin(),merged.begin()+len,merged.end());
      len+=G_snps[j].size();

      std::set_intersection(T_snps[i].begin(), T_snps[i].end(), 
        T_snps[j].begin(), T_snps[j].end(), std::back_inserter(inter));
      
      merged.insert( merged.end(), T_snps[i].begin(), T_snps[i].end());
      std::inplace_merge (merged.begin(),merged.begin()+len,merged.end());
      len+=T_snps[i].size();
      merged.insert( merged.end(), T_snps[j].begin(), T_snps[j].end());
      std::inplace_merge (merged.begin(),merged.begin()+len,merged.end());

      // find union
      std::sort(inter.begin(), inter.end());
      int d1 = std::unique(merged.begin(), merged.end()) - merged.begin();
      int d2 = std::unique(inter.begin(), inter.end()) - inter.begin();


      // std::cout << '\n';
      // for (std::vector<int>::const_iterator p = inter.begin(); p != inter.end(); ++p)
      //   std::cout << *p << ' ';
      // std::cout << '\n';


      // std::cout << 'i' << i << '\t';
      // std::cout << 'j' << j << '\t';
      // std::cout << d1<< '\t';
      // std::cout << d2<< '\t';
      // std::cout << d1-d2<< '\t';
      // std::cout << '\n';

      comp_snps[j] = d1 - d2;

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
      for (int j=start; j < n_seqs; j++) {
        if ((dist==-1) || (comp_snps[j] <= dist)){
          printf("%d%c%d%c%d\n", i, sep, j, sep, comp_snps[j]);
        }
      }
    } else {
      printf("%s", seq_names[i].c_str());
      for (int j=0; j < n_seqs; j++) {
        printf("%c%d", sep, comp_snps[j]);
      }
      printf("\n");
    }
  }

  return 0;

  }
