#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <math.h>

#define ENOUGH 100   /* How many times has to be n_vals>tau_int to be considered reasonable */

/***************************
 *    Global Variables     *
 ***************************/

char progname[]="analyze"; // name of the program
double *data=NULL;         // array of data values


/**************************
*     Functions           *
***************************/

/* read-in the data */
int read_data (char *filename, int column, int discard);

/* binning analysis */
void analyze_bins (int n_blocks, int n_vals);

/********************************
 *        Implementation        *
 ********************************/

int read_data (char *filename, int column, int discard) { //{{{

  int i, j, // index variables
      n_vals = 0;
  FILE *f = fopen(filename, "r");

  if(f == NULL) {
    fprintf(stderr, "%s error: cannot open file %s", progname, filename);
    exit(1);
  }

  // first read-in the simulation parameters
  fscanf(f,"%d", &n_vals);
//printf("\nRead from input file: n_vals=%d\n", n_vals);

  n_vals -= discard; // discard what shoud not be used

//printf("\nDiscarding first %d samples; using %d samples for analysis.\n", discard, n_vals);

  // the allocate appropriate space
  data = (double*)malloc(n_vals*sizeof(double));

  // set all values to zero
  for(i = 0; i < n_vals; i++)
    data[i] = 0.0;

  // leave out first 'discard' lines - file pointer is now before the end of the commented line
  for (i = 0; i < (discard+1); i++) {
    while (getc(f) != '\n')
      ;
  }

  // read-in the actual data
  for (i = 0; i < n_vals; i++) {
    // skip previous columns
    for(j = 1; j < column; j++)
      fscanf(f, "%*f");

    // read data
    fscanf(f, "%lf", &data[i]);

    while (getc(f) != '\n')
      ;
  }

  return n_vals;
} //}}}

void analyze_bins (int n_blocks, int n_vals) { //{{{

  int i, // index variables
      block = 0, block_end = 0, block_size = 0, remainder = 0;
  double av = 0.0, av2 = 0.0, err = -1.0, tau_int = -1.0, // average, average_squared, error and autocorrelation time
         *av_block, // per-block averages
         block_stdev = 0.0; // and standard deviation of block averages

  // allocate space for individual block averages
  av_block = (double*)malloc(n_blocks*sizeof(double));

  for(i = 0; i < n_blocks; i++)
    av_block[i] = 0.0;

  remainder = n_vals % n_blocks;

  block_size = (n_vals - remainder) / n_blocks;
  block_end = remainder + block_size;

  // collect block averages
  for(i = remainder; i < n_vals; i++) {
    av += data[i];
    av2 += data[i]*data[i];
    av_block[block] += data[i];

    if(i == block_end) {
      block++;
      block_end += block_size;
    }
  }

  // now n_vals becomes the number of really considered values
  n_vals -= remainder;

  // divide by number of samples
  av /= n_vals;
  av2 /= n_vals;
  for(i = 0; i < n_blocks; i++) {
    av_block[i] /= block_size;
    block_stdev += (av_block[i] - av)*(av_block[i] - av);
  }

  block_stdev /= (n_blocks - 1);
  err = sqrt(block_stdev/n_blocks);
  tau_int = 0.5 * block_size * block_stdev / (av2 - av*av);

  printf("%4d %lf %lf %lf\n",n_blocks,av,err,tau_int);

  return;
} //}}}

int main ( int argc, char** argv ) {

  if (argc != 5) /* wrong number of parameters! */
    exit(fprintf(stderr,"\n Usage: %s <filename> <column> <discard> <n_blocks>\n\n\
    Analyze output of the ising model simulation. Data stored in file\n\
    <filename> should contain three columns with values: m, m^2 and e\n\
    and a header with simulation parameters: temperature, lattice and \n\
    number of samples.\n\n\
    Arguments:\n\
    <filename> [string]  name of data file\n\
    <column>   [int]  number of column in file containing the data to analyze\n\
    <discard>  [int]  number of data values considered as equilibration\n\
    <n_blocks> [int]  number of blocks for binning and jackknife\n\
    \nWARNING: no sanity checks on argument values are performed. If you mess them up, undefined behaviour occurs.\n\
    \n", progname));

  /* Assign all variables from given parameters: */
  char *filename = argv[1];
  int column = atoi(argv[2]);     // number of column to analyze
  int discard = atoi(argv[3]);    // number of samples discarded from the begining
  int n_blocks = atoi(argv[4]);   // number of bins for binning and jackknife

  /* For feedback, write down what we've read in */
//printf("Input parameters: filename: %s, column=%d, discard=%d, n_blocks=%d\n", filename, column, discard, n_blocks);

  // Read the data
  int n_vals = read_data(filename, column, discard);
  if(!n_vals) {
    fprintf(stderr,"\nCannot read data from file %s\n\n",filename);
    exit(1);
  }

  // And analyze them using bins
  analyze_bins(n_blocks, n_vals);

  return 0;
}
