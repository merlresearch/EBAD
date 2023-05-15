// Copyright (C) 2014, 2023 Mitsubishi Electric Research Laboratories (MERL)
//
// SPDX-License-Identifier: AGPL-3.0-or-later


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define INF 1e20       //Pseudo Infitinte number for this code

//using namespace std;

/// Data structure for sorting the query.
typedef struct Index
    {   double value;
        int    index;
    } Index;


/// Comparison function for sorting the query.
/// The query will be sorted by absolute z-normalization value, |z_norm(Q[i])| from high to low.
int comp(const void *a, const void* b)
{   Index* x = (Index*)a;
    Index* y = (Index*)b;
    return abs(y->value) - abs(x->value);
}


/// Main function for calculating ED distance between the query, Q, and current data, T.
/// Note that Q is already sorted by absolute z-normalization value, |z_norm(Q[i])|
//double distance(const double * const Q, const double * const T ,
//const int& j , const int& m , const double& mean , const double& std
//, const int* const order, const double& bsf)
double distance(double *Q, double *T , const int j , const int m , const double mean , const double std , int *order, const double bsf)

{
    int i;
    double sum = 0;
    for ( i = 0 ; i < m && sum < bsf ; i++ )
    {
//        double x = (T[(order[i]+j)]-mean)/std;
        double x = T[(order[i]+j)];
        sum += (x-Q[i])*(x-Q[i]);
    }
//    printf("number of iterations completed in distance() =%d, j=%d, bsf=%lf, sum=%lf\n", i, j, bsf, sum);

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
    {
        printf("ERROR: Invalid Number of Arguments!!!\n");
    }
    exit(1);
}


int main(  int argc , char *argv[] )
{
  FILE *fp;              // the input file pointer
  FILE *qp;              // the query file pointer
  double *Q;             // test array
  double *Qs;            // query array
  double *T;             // array of current data
  int *order;            // ordering of query by |z(q_i)|
  double bsf;            // best-so-far
  int m_test;            // length of test series
  int m_train;           // length of train series
  long long loc = 0;     // answer: location of the best-so-far match
  int w;                 // length of testing subsequence
  int hw;
  double d;
  double dum1, dum2;
  long long i, j, k;
  double ex, ex2, mean, std;
  double *anomaly_score;
  double t1,t2;
  int len;

  // read in training time series
  // read in testing time series
  // for each subsequence in testing time series
  //    initialize bsf with distance to last NN
  //    sort it (or update previous sort with binary search)
  //    call ED distance function (with early stopping)

  t1 = clock();

  bsf = INF;
  i = 0;
  j = 0;
  ex = ex2 = 0;

  if (argc<=5) {
    printf("Usage: BF_AnomalyDet_UCR_ED train.txt test.txt train_length test_length subsequence_length\n");
    error(4);
  }

  fp = fopen(argv[1],"r");
  if( fp == NULL )
    exit(2);

  qp = fopen(argv[2],"r");
  if( qp == NULL )
    exit(2);

  m_train = atol(argv[3]);
  m_test = atol(argv[4]);
  w = atoi(argv[5]);
  hw = (int)w/2;

  // Array for keeping the test data
  Q = (double *)malloc(sizeof(double)*m_test);
  if( Q == NULL )
    error(1);

  // Array for keeping the query data
  Qs = (double *)malloc(sizeof(double)*w);
  if( Qs == NULL )
    error(1);

  // Array for keeping the query data
  anomaly_score = (double *)calloc(m_test, sizeof(double));
  if( anomaly_score == NULL )
    error(1);

  // Read the test data from input file and calculate its statistic such as mean, std
  while(fscanf(qp,"%lf", &d) != EOF && i < m_test)
  {
    Q[i] = d;
    i++;
  }
  fclose(qp);

  T = (double *)malloc(sizeof(double)*m_train);
  if( T == NULL )
    error(1);

  i = 0;
  // Read training data file, one value at a time
  while(fscanf(fp,"%lf", &d) != EOF )
  {
    T[i] = d;
    i++;
  }
  fclose(fp);

  double dist = 0;

  order = (int *)malloc(sizeof(int)*w);
  if( order == NULL )
    error(1);

  Index *Q_tmp = (Index *)malloc(sizeof(Index)*w);
  if( Q_tmp == NULL )
    error(1);

  for(k=0; k<m_test-w; k++)
  {
    if ((k % 1000) == 0)
      printf("k=%lld\n", k);

    // Sort the query data
    for( j=0 ; j < w ; j++ )
    {
      Q_tmp[j].value = fabs(Q[k+j]);
      Q_tmp[j].index = j;
    }
    qsort(Q_tmp, w, sizeof(Index),comp);
    for( j=0; j<w; j++)
    {
      order[j] = Q_tmp[j].index;
      Qs[j] = Q[k+order[j]];
    }

    if (k > 0) {
      bsf = distance(Qs, T, loc, w, mean, std, order, INF);
    }

    for (i=0; i<m_train-w; i++)
    {
      /// Calculate ED distance
      dist = distance(Qs,T,i,w,mean,std,order,bsf);
      if( dist < bsf )
      {
	bsf = dist;
	loc = i;
      }

    }
    anomaly_score[k+hw] = bsf;
  }

  t2 = clock();

  printf("Total Execution Time : %lf sec\n", (t2-t1)/CLOCKS_PER_SEC);
  free(Q_tmp);

  len = strlen(argv[2]);
  argv[2][len-3] = 'o';
  argv[2][len-2] = 'u';
  argv[2][len-1] = 't';
  printf("out file = %s\n", argv[2]);
  if ((fp = fopen(argv[2], "w")) == NULL)
  {
    printf("Cannot open file %s to write\n", argv[2]);
    exit(-1);
  }

  for (i=0; i<m_test; i++)
  {
    fprintf(fp, "%lf ", anomaly_score[i]);
  }
  fclose(fp);
}
