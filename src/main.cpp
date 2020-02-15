#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <fcntl.h>
#include <vector>
#include <armadillo>
#include "kseq.h"

#ifdef _OPENMP
  #include <omp.h>
  void omp_set_num_threads(int num_threads);
  int omp_get_num_threads();
#endif

#include <time.h>
#include <ctime>

#define VERSION "0.1.0"
#define EXENAME "pairsnp"
using namespace arma;

void show_help(int retcode)
{
  FILE* out = (retcode == EXIT_SUCCESS ? stdout : stderr);
  fprintf(out, "SYNOPSIS\n  Pairwise SNP similarity and distance matrices using fast matrix algerbra libraries\n");
  fprintf(out, "USAGE\n  %s [options] alignment.fasta[.gz] > matrix.csv\n", EXENAME);
  fprintf(out, "OPTIONS\n");
  fprintf(out, "  -h\tShow this help\n");
  fprintf(out, "  -v\tPrint version and exit\n");
  fprintf(out, "  -s\tOutput in sparse matrix form (i,j,distance).\n");
  fprintf(out, "  -d\tDistance threshold for sparse output. Only distances <= d will be returned.\n");
  fprintf(out, "  -c\tOutput CSV instead of TSV\n");
  fprintf(out, "  -n\tCount comparisons with Ns (off by default)\n");
  // fprintf(out, "  -t\tNumber of threads to use\n");
  fprintf(out, "  -b\tBlank top left corner cell instead of '%s %s'\n", EXENAME, VERSION);
  exit(retcode);
}

int main(int argc, char *argv[])
{
  // clock_t begin_time = clock();

  // parse command line parameters
  int opt, quiet=0, csv=0, corner=1, allchars=0, keepcase=0, sparse=0, count_n=0, dist;
  while ((opt = getopt(argc, argv, "hqsncvd:")) != -1) {
    switch (opt) {
      case 'h': show_help(EXIT_SUCCESS); break;
      case 'q': quiet=1; break;
      case 'c': csv=1; break;
      case 's': sparse=1; break;
      case 'n': count_n=1; break;
      case 'd': dist=atoi(optarg); break;
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


  // cout << "warmup. " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
  // begin_time = clock();

  //Initially run through fasta to get consensus sequence and dimensions of matrix
  int n = 0;
  int l = 0;
  kseq seq;
  FunctorRead r;
  kstream<int, FunctorRead> ks(fp, r);

  // Use first sequence to set size of arrays
  l = ks.read(seq);
  int seq_length = seq.seq.length();
  std::vector<std::vector<int>> allele_counts(5, std::vector<int>(seq_length));

  do {
    for(int j=0; j<seq_length; j++){
      if((seq.seq[j]=='a') || (seq.seq[j]=='A')){
        allele_counts[0][j] += 1;
      } else if((seq.seq[j]=='c') || (seq.seq[j]=='C')){
        allele_counts[1][j] += 1;
      } else if((seq.seq[j]=='g') || (seq.seq[j]=='G')){
        allele_counts[2][j] += 1;
      } else if((seq.seq[j]=='t') || (seq.seq[j]=='T')){
        allele_counts[3][j] += 1;
      } else {
        allele_counts[4][j] += 1;
      }
    }
    n++;
  } while((l = ks.read(seq)) >= 0);
  close(fp);


  // cout << "initial read.. " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
  // begin_time = clock();

  // Now calculate the consensus sequence
  uvec consensus(seq_length);
  int n_seqs = n;
  std::vector<std::string> seq_names;

  int max_allele = -1;

  for(int j=0; j<seq_length; j++){
    for(int i=0; i<5; i++){
      if(allele_counts[i][j]>max_allele){
        max_allele = allele_counts[i][j];
        consensus(j) = i;
      }
    }
    max_allele = -1;
  }

  // Find the total number of SNPs we add a buffer of 
  // 'seq_length' to deal with the skipped first line
  int n_total_snps = seq_length*n_seqs;
  int total_A_snps = 0;
  int total_C_snps = 0;
  int total_G_snps = 0;
  int total_T_snps = 0;
  int total_N_snps = 0;
  for(int j=0; j<seq_length; j++){
    n_total_snps =  n_total_snps - allele_counts[consensus(j)][j];
    if ((consensus(j)!=0)) total_A_snps += allele_counts[0][j];
    if ((consensus(j)!=1)) total_C_snps += allele_counts[1][j];
    if ((consensus(j)!=2)) total_G_snps += allele_counts[2][j];
    if ((consensus(j)!=3)) total_T_snps += allele_counts[3][j];
    if ((consensus(j)!=4)) total_N_snps += allele_counts[4][j];
  }


  // create matrices to store SNP locations prior to generating
  umat locationsA(2, total_A_snps, fill::zeros);
  umat locationsC(2, total_C_snps, fill::zeros);
  umat locationsG(2, total_G_snps, fill::zeros);
  umat locationsT(2, total_T_snps, fill::zeros);
  umat locationsN(2, total_N_snps, fill::zeros);
  
  // open a new stream to the fasta file TODO: I think there's a cleaner way of doing this.
  int fp2 = open(fasta, O_RDONLY);
  kstream<int, FunctorRead> ks2(fp2, r);
  char temp_char;
  int n_snps = 0;
  int nA_snps = 0;
  int nC_snps = 0;
  int nG_snps = 0;
  int nT_snps = 0;
  int nN_snps = 0;
  n=0;
  std::string nt ("AaCcGgTt");

  while((l = ks2.read(seq)) >= 0) {
    // Record sequence names
    // seq.name.replace(">","")
    auto it = std::find(seq.name.begin(), seq.name.end(), '>');
    if (it != seq.name.end()){
      seq.name.erase(it);
    }

    if(seq_length!=seq.seq.length()){
      std::cout << "Error: sequences are not all of the same length!" << std::endl;
      return 1;
    }

    seq_names.push_back(seq.name);

    for(int j=0; j<seq_length; j++){
      temp_char = seq.seq[j];
      if((((temp_char=='A') || (temp_char=='a')) && (consensus(j)!=0))){
        locationsA(1, nA_snps) = n;
        locationsA(0, nA_snps) = j; 
        nA_snps += 1;
        n_snps += 1;
      } else if((((temp_char=='C') || (temp_char=='c')) && (consensus(j)!=1))){
        locationsC(1, nC_snps) = n;
        locationsC(0, nC_snps) = j; 
        nC_snps += 1;
        n_snps += 1;
      } else if((((temp_char=='G') || (temp_char=='g')) && (consensus(j)!=2))){
        locationsG(1, nG_snps) = n;
        locationsG(0, nG_snps) = j; 
        nG_snps += 1;
        n_snps += 1;
      } else if((((temp_char=='T') || (temp_char=='t')) && (consensus(j)!=3))){
        locationsT(1, nT_snps) = n;
        locationsT(0, nT_snps) = j; 
        nT_snps += 1;
        n_snps += 1;
      } else if((nt.find(temp_char) == std::string::npos && (consensus(j)!=4))){
        locationsN(1, nN_snps) = n;
        locationsN(0, nN_snps) = j; 
        nN_snps += 1;
        n_snps += 1;
      }
    }
    n += 1;
  }

  sp_umat sparse_matrix_A(locationsA, ones<uvec>(nA_snps), seq_length, n_seqs);
  sp_umat sparse_matrix_C(locationsC, ones<uvec>(nC_snps), seq_length, n_seqs);
  sp_umat sparse_matrix_G(locationsG, ones<uvec>(nG_snps), seq_length, n_seqs);
  sp_umat sparse_matrix_T(locationsT, ones<uvec>(nT_snps), seq_length, n_seqs);
  sp_umat sparse_matrix_N(locationsN, ones<uvec>(nN_snps), seq_length, n_seqs);

  sp_umat binary_snps = spones(sparse_matrix_A + sparse_matrix_C + sparse_matrix_G + sparse_matrix_T + sparse_matrix_N);
  sp_umat t_binary_snps = binary_snps.t();

  // // build some matrices we'll need to deal with column where the consensus is 'N'
  umat total_n = umat(sum(sparse_matrix_N, 0));
  uvec cons_idsN = find(consensus == 4); // Find indices
  umat matrix_n_cols = zeros<umat>(n_seqs, cons_idsN.size());

  if(!count_n){
    for (int i=0; i<cons_idsN.size(); i++){
      sp_umat col(t_binary_snps.col(cons_idsN(i)));
      for (arma::sp_umat::iterator it = col.begin(); it != col.end(); ++it) {
        matrix_n_cols(it.row(), i) = 1;
      }
    }
  }

  uvec tot_cons_snps_N = uvec(sum(matrix_n_cols, 1));
  matrix_n_cols = matrix_n_cols.t();

  sp_umat t_sparse_matrix_A = sparse_matrix_A.t();
  sp_umat t_sparse_matrix_C = sparse_matrix_C.t();
  sp_umat t_sparse_matrix_G = sparse_matrix_G.t();
  sp_umat t_sparse_matrix_T = sparse_matrix_T.t();
  sp_umat t_sparse_matrix_N = sparse_matrix_N.t();

  umat total_snps = umat(sum(sparse_matrix_A + sparse_matrix_C + sparse_matrix_G + sparse_matrix_T + sparse_matrix_N, 0));


  #pragma omp parallel for ordered shared(t_sparse_matrix_A, sparse_matrix_A \
    , t_sparse_matrix_C, sparse_matrix_C \
    , t_sparse_matrix_G, sparse_matrix_G \
    , t_sparse_matrix_T, sparse_matrix_T \
    , t_sparse_matrix_N, sparse_matrix_N \
    , t_binary_snps, binary_snps \
    , matrix_n_cols, tot_cons_snps_N \
    , seq_names, dist, sparse, sep, total_n \
    , count_n, total_snps, n_seqs) default(none) schedule(static,1)
  for (size_t i = 0; i < n_seqs; i++) {

      umat comp_snps = umat(t_sparse_matrix_A * sparse_matrix_A.col(i));
      comp_snps = comp_snps + umat(t_sparse_matrix_C * sparse_matrix_C.col(i));
      comp_snps = comp_snps + umat(t_sparse_matrix_G * sparse_matrix_G.col(i));
      comp_snps = comp_snps + umat(t_sparse_matrix_T * sparse_matrix_T.col(i));

      umat diff_n = umat(t_sparse_matrix_N * sparse_matrix_N.col(i));
      comp_snps = comp_snps + diff_n;

      umat differing_snps = umat(n_seqs, 1);
      for(int j=0; j<n_seqs; j++){
        differing_snps(j,0) = total_snps(i) + total_snps(j);
      }
      
      if(count_n){
        comp_snps = differing_snps - umat(t_binary_snps * binary_snps.col(i)) - comp_snps;
      } else {

        umat  cons_snps_N = matrix_n_cols.col(i).t() * matrix_n_cols;
        
        for(int j=0; j<n_seqs; j++){
          diff_n(j,0) = total_n(i) + total_n(j) - 2*diff_n(j,0) + tot_cons_snps_N(i) + tot_cons_snps_N(j) - 2*cons_snps_N(j);
        }
        comp_snps = differing_snps - umat(t_binary_snps * binary_snps.col(i)) - comp_snps - diff_n;   
      }        

    // Output the distance matrix to stdout
    #pragma omp ordered //#pragma omp critical
    if (sparse){
      for (int j=0; j < n_seqs; j++) {
        if (int(comp_snps(j)) <= dist){
          printf("%d%c%d%c%d\n", i, sep, j, sep, int(comp_snps(j)));
        }
      }
    } else {
      printf("%s", seq_names[i].c_str());
      for (int j=0; j < n_seqs; j++) {
        printf("%c%d", sep, int(comp_snps(j)));
      }
      printf("\n");
    }
  }

  return 0;

  }
