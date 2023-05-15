// Copyright (C) 2014, 2023 Mitsubishi Electric Research Laboratories (MERL)
//
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define INF 1e20       //Pseudo Infitinte number for this code

#define STRIDE 2

int fv_size;
int fv_size_half;
int smoothing_width;


extern void compute_SST_feature(double *T, int indx, int w, double *fv, int smoothing_width);

// computes the distance between an SST exemplar and an SST feature
// The distance is the number of standard deviations between the mean
// of the exemplar and the SST feature minus 3.  Also any distance
// less than 0 is zeroed.
double distance(double *fv, double **E , const int j , const int fv_size_half, const double bsf)
{
    int i;
    double sum = 0;
    double c = ((double)(fv_size_half-7))/7.0;

    // distance due to the statistical components of the SST feature
    for ( i = 0; (i < fv_size_half-7) && (sum < bsf); i++ ) {
        sum += MAX(0.0, (E[j][i] - fv[i])/E[j][i + fv_size_half] - 3.0); // 4 works well
    }

    // distance due to the trajectory components of the SST feature
    for ( i = fv_size_half-7; (i < fv_size_half) && (sum < bsf); i++ ) {
        sum += (MAX(0.0, (E[j][i] - fv[i])/E[j][i + fv_size_half] - 3.0)*c);
    }

    return sum;
}


/// If serious error happens, terminate the program.
void error(int id)
{
    if(id==1)
        printf("ERROR : Memory can't be allocated!!!\n\n");
    else if ( id == 2 )
        printf("ERROR : File not Found!!!\n\n");
    else if ( id == 3 )
        printf("ERROR : Can't create Output File!!!\n\n");
    else if ( id == 4 )
        printf("ERROR: Invalid Number of Arguments!!!\n");

    exit(1);
}


int main(  int argc , char *argv[] )
{
  FILE *fp;              // the input file pointer
  FILE *qp;              // the query file pointer
  double *Q;             // test array
  double **E;            // matrix of exemplars
  double bsf;            // best-so-far
  int m_test;            // length of test series
  long long loc = 0;     // answer: location of the best-so-far match
  int w, w2;                 // length of testing subsequence
  int hw;
  double d;
  long long i, j, k;
  double ex, ex2;
  double *anomaly_score;
  int n_ex;
  double *fv1;
  double t1,t2;
  int len;
  int smoothing_width;

  // read in testing time series
  // read in the exemplars (learned using MERL_model_efficient.c)
  // for each subsequence in testing time series
  //    compute SST feature
  //    initialize bsf with distance to last nearest neighbor exemplar
  //    for each exemplar
  //       call ED distance function (with early stopping) to
  //       find distance between exemplar and SST test feature
  //       store index and distance to nearest exemplar

  t1 = clock();

  bsf = INF;
  i = 0;
  j = 0;
  ex = ex2 = 0;

  if (argc<=4) {
    printf("Usage: BF_AnomalyDet_MERL_SST model.txt test.txt test_length subsequence_length\n");
    error(4);
  }

  fp = fopen(argv[1],"r");  // exemplar file
  if (fp == NULL)
    error(2);

  qp = fopen(argv[2],"r");  // testing time series file
  if (qp == NULL)
    error(2);

  m_test = atol(argv[3]);
  w = atoi(argv[4]);
  hw = (int)w/2;
  smoothing_width = (int)((float)w/40.0 + 0.5);  // parameter used to compute
                                                 // SST features

  // Array for keeping the test data
  Q = (double *)malloc(sizeof(double)*m_test);
  if (Q == NULL)
    error(1);

  // Array for keeping the anomaly scores
  anomaly_score = (double *)calloc(m_test, sizeof(double));
  if (anomaly_score == NULL)
    error(1);

  // Read the test data from input file
  while (fscanf(qp,"%lf", &d) != EOF && i < m_test) {
    Q[i] = d;
    i++;
  }

  fclose(qp);

  // read in the exemplars from the model file
  fscanf(fp, "%d %d", &n_ex, &fv_size);
  w2 = STRIDE*((fv_size/2) - 7);
  fv_size_half = fv_size/2;

  if (w2 != w) {
    printf("The size of training window (in the model file) %d should be the same as the size of testing window %d!\n", w2, w);
    exit(-1);
  }

  E = (double **)malloc(n_ex * sizeof(double *));
  if( E == NULL )
    error(1);

  for (i=0; i<n_ex; i++)
    E[i] = (double *)malloc(fv_size * sizeof(double));

  for (i=0; i<n_ex; i++) {
    for (j=0; j<fv_size; j++) {
      fscanf(fp, "%lf", &d);
      E[i][j] = d;
    }
  }
  fclose(fp);

  double dist = 0;

  fv1 = (double *)calloc(fv_size, sizeof(double));

  // loop over all testing subsequences and find nearest exemplar to each
  for (k=smoothing_width; k<m_test-w-smoothing_width; k+=1) {

    if ((k % 1000) == 0)
      printf("k=%lld\n", k);

    compute_SST_feature(Q, k, w, fv1, smoothing_width);

    if (k > 0) {
      bsf = distance(fv1, E, loc, fv_size_half, INF);
    }

    for (i=0; i<n_ex; i++) {

      /// Calculate distance
      dist = distance(fv1, E, i, fv_size_half, bsf);

      if( dist < bsf ) {
	bsf = dist;
	loc = i;
      }
    }

    anomaly_score[k+hw] = bsf;
  }

  t2 = clock();

  printf("Total Execution Time : %lf sec\n", (t2-t1)/CLOCKS_PER_SEC);

  len = strlen(argv[2]);
  argv[2][len-3] = 'o';
  argv[2][len-2] = 'u';
  argv[2][len-1] = 't';
  printf("out file = %s\n", argv[2]);
  if ((fp = fopen(argv[2], "w")) == NULL) {
    printf("Cannot open file %s to write\n", argv[2]);
    exit(-1);
  }

  for (i=0; i<m_test; i+=1) {
    fprintf(fp, "%lf ", anomaly_score[i]);
  }
  fclose(fp);
}
