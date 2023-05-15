// Copyright (C) 2014, 2023 Mitsubishi Electric Research Laboratories (MERL)
//
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define INF 1e20       //Pseudo Infitinte number for this code

#define MIN(a, b) (a < b) ? a : b
#define MAX(a, b) (a > b) ? a : b

int fv_size;
int fv_size_half;
int smoothing_width;

extern void compute_SST_feature(double *T, int indx, int w, double *fv, int m);
void exemplar_selection_hierarchical(double *T, int m, int w, double ***E, int *n_ex);

int main(int argc, char *argv[])
{
  FILE *fp;              // the training file pointer
  double *T;             // array of current data
  int m_train, w;
  double **E;
  int n_ex;
  int i, j;
  double d, dum1, dum2;
  double t1,t2;

  if (argc != 5) {
    printf("Usage: MERL_model_efficient training.txt length_of_training_series window_length out_fname\n");
    exit(-1);
  }

  t1 = clock();

  // read in training series
  if ((fp = fopen(argv[1],"r")) == NULL) {
    printf("Cannot open file %s to read.\n", argv[1]);
    exit(-1);
  }

  m_train = atoi(argv[2]);
  w = atoi(argv[3]);
  fv_size = 2*((int)(w/2 + 0.5) + 7);
  fv_size_half = fv_size/2;

  smoothing_width = (int)((float)w/40.0 + 0.5);

  printf("m_train=%d, w=%d, fv_size = %d\n", m_train, w, fv_size);

  T = (double *)malloc(sizeof(double)*m_train);
  if( T == NULL )
    exit(-1);

  i=0;
  // Read the training time series from input file
  while(fscanf(fp,"%lf", &d) != EOF && i < m_train) {
    T[i] = d;
    i++;
  }

  fclose(fp);

  // call exemplar selection
  exemplar_selection_hierarchical(T, m_train, w, &E, &n_ex);
  printf("done selecting exemplars\n");

  // save exemplars
  if ((fp = fopen(argv[4],"w")) == NULL) {
    printf("Cannot open file %s to write.\n", argv[4]);
    exit(-1);
  }

  fprintf(fp, "%d %d\n", n_ex, fv_size);

  for (i=0; i<n_ex; i++) {
    for (j=0; j<fv_size; j++) {
      fprintf(fp, "%lf ", E[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  t2 = clock();

  printf("Total Execution Time : %lf sec\n", (t2-t1)/CLOCKS_PER_SEC);
}


double eucl_dist_weighted(double *fv1, double *fv2, int len)
{
  int i;
  double d=0.0;
  int w = len - 7;

  // compute distance between statistical components of SST feature
  for (i=w; i<len; i++)
    d += (fv1[i]-fv2[i])*(fv1[i]-fv2[i]);
  d *= (double)w/7.0;

  // compute distance between trajectory components of SST feature
  for (i=0; i<w; i++)
    d += (fv1[i]-fv2[i])*(fv1[i]-fv2[i]);

  return d;
}


double eucl_dist_weighted_earlystop(double *fv1, double *fv2, int len, double thresh)
{
  int i;
  double d=0.0;
  int w = len - 7;
  double mult_factor = (double)w/7.0;

  // compute distance between statistical components of SST feature
  for (i=w; ((i<len) && (d < thresh)); i++)
    d += (fv1[i]-fv2[i])*(fv1[i]-fv2[i])*mult_factor;

  // compute distance between trajectory components of SST feature
  for (i=0; ((i<w) && (d < thresh)); i++)
    d += (fv1[i]-fv2[i])*(fv1[i]-fv2[i]);

  return d;
}


void min_val(double *NN_dist, int n_exemplars, double *min_d, int *min_i)
{
  int i;

  *min_d = INF;

  for (i=0; i<n_exemplars; i++) {
    if (NN_dist[i] < *min_d) {
      *min_d = NN_dist[i];
      *min_i = i;
    }
  }
}


void nearest_neighbor_chunk(double **E, int *count, int j, int start_ind, int end_ind, int *indx, double *dist)
{
  // find the nearest neighbor to indx in E
  // use early stopping

  int i;
  double d;
  *dist = INF;

  for (i=start_ind; i<end_ind; i++) {
    if ((count[i] > 0) && (i != j)) {
      d = eucl_dist_weighted_earlystop(E[j], E[i], fv_size_half, *dist);
      if (d < *dist) {
	*dist = d;
	*indx = i;
      }
    }
  }
}


double initial_merge(double *T, int m, int w, double ***E, int **count, int *n_ex)
{
  int maxOffset;
  int Nsamples = MIN(1000, m);

  double sum_dist = 0.0;
  double sumsqr_dist = 0.0;
  int cnt = 0;
  int n_exemplars;

  double mu, sigma, thresh;
  int indx, i, j, k, ii;
  double x;

  double *fv1, *fv2;
  double dist;

  printf("In initial merge\n");

  if (w <= 100)
    maxOffset = 2;
  else if (w <= 200)
    maxOffset = 3;
  else if (w <= 300)
    maxOffset = 4;
  else if (w <= 800)
    maxOffset = 5;
  else
    maxOffset = 6;

  // allocate space for 32 exemplars to begin
  (*E) = (double **)malloc(32*sizeof(double *));
  (*count) = (int *)malloc(32*sizeof(int));

  for (i=0; i<32; i++)
    (*E)[i] = (double *)calloc(fv_size, sizeof(double));
  n_exemplars = 32;

  fv1 = (double *)calloc(fv_size, sizeof(double));
  fv2 = (double *)calloc(fv_size, sizeof(double));

  // use a random sample of neighborhoods to compute mean and variance of
  // neighborhood distances to get a reasonable threshold value
  for (k=0; k<Nsamples; k++) {
    i = random() % (m-w-2*smoothing_width-2*maxOffset-1) + maxOffset+smoothing_width;

    compute_SST_feature(T, i, w, fv1, smoothing_width);

    compute_SST_feature(T, i+maxOffset, w, fv2, smoothing_width);
    dist = eucl_dist_weighted(fv1, fv2, fv_size_half);

    sum_dist += dist;
    sumsqr_dist += dist*dist;
    cnt++;

    compute_SST_feature(T, i-maxOffset, w, fv2, smoothing_width);
    dist = eucl_dist_weighted(fv1, fv2, fv_size_half);

    sum_dist += dist;
    sumsqr_dist += dist*dist;
    cnt++;
  }

  mu = sum_dist/(double)cnt;
  sigma = sqrt(sumsqr_dist/(double)cnt - mu*mu);
  thresh = mu + 3.0 * sigma;
  printf("mu = %lf, sigma = %lf, thresh = %lf\n", mu, sigma, thresh);

  indx = 0;
  i = smoothing_width;

  // Now we have a threshold that tells us the max distance between
  // feature vectors for merging.  Loop through and merge consecutive
  // feature vectors within the threshold distance

  while (i < m-w-smoothing_width) {
    j = 1;
    compute_SST_feature(T, i, w, fv1, smoothing_width);
    compute_SST_feature(T, i+j, w, fv2, smoothing_width);
    x = eucl_dist_weighted_earlystop(fv1, fv2, fv_size_half, thresh);

    while ((i+j < m-w+1) && (x < thresh)) {
      j = j + 1;
      compute_SST_feature(T, i+j, w, fv2, smoothing_width);
      x = eucl_dist_weighted_earlystop(fv1, fv2, fv_size_half, thresh);
    }

    if (x >= thresh)
      j = j - 1;

    k = j+1;

    if (i+k+w-1 < m) {
      compute_SST_feature(T, i+j, w, fv1, smoothing_width);
      compute_SST_feature(T, i+k, w, fv2, smoothing_width);
      x = eucl_dist_weighted_earlystop(fv1, fv2, fv_size_half, thresh);

      while ((i+k < m-w) && (x < thresh)) {
	k = k + 1;
	compute_SST_feature(T, i+k, w, fv2, smoothing_width);
	x = eucl_dist_weighted_earlystop(fv1, fv2, fv_size_half, thresh);
      }
    }
    else
      k=j;

    if ((k>j) && (x >= thresh))
      k = k - 1;

    // now merge all windows from i to i+k

    // check if indx>n_exemplars and realloc if needed
    if (indx >= n_exemplars) {
      n_exemplars *= 2;
      (*E) = (double **)realloc((*E), n_exemplars*sizeof(double *));
      (*count) = (int *)realloc((*count), n_exemplars*sizeof(int));

      for (ii=n_exemplars/2; ii<n_exemplars; ii++)
	(*E)[ii] = (double *)calloc(fv_size, sizeof(double));
    }

    for (j=i; j<=i+k; j++) {
      compute_SST_feature(T, j, w, fv1, smoothing_width);
      for (ii=0; ii<fv_size_half; ii++) {
	(*E)[indx][ii] += fv1[ii];
	(*E)[indx][ii+fv_size_half] += fv1[ii]*fv1[ii];
      }

      //indices = [indices m];
    }

    for (ii=0; ii<fv_size_half; ii++) {
      (*E)[indx][ii] /= (double)(k+1);  // mean
    }

    (*count)[indx] = k+1;

    indx = indx + 1;
    i = i+k+1;
  }

  n_exemplars = indx;
  (*n_ex) = n_exemplars;

  (*E) = (double **)realloc((*E), n_exemplars*sizeof(double *));
  (*count) = (int *)realloc((*count), n_exemplars*sizeof(int));
  printf("total exemplars after initial merge = %d\n", n_exemplars);
  printf("thresh = %lf\n", thresh);

  free(fv1);
  free(fv2);

  return thresh;
}


void exemplar_selection_chunk(double **E, int *count, int start_ind, int end_ind, double thresh)
{
  double *fv1, *fv2;
  double dist;

  int *NN;
  double *NN_dist;
  int min_i, min_i2;
  double min_d;
  int i, indx;
  int mode;

  int n_exemplars = end_ind - start_ind;

  // for each exemplar, find nearest neighbor and L2 dist
  NN = (int *)malloc(n_exemplars*sizeof(int));
  NN_dist = (double *)malloc(n_exemplars*sizeof(double));

  for (i=0; i<n_exemplars; i++) {
    nearest_neighbor_chunk(E, count, start_ind+i, start_ind, end_ind, &indx, &dist);
    NN[i] = indx;
    NN_dist[i] = dist;
  }

  // find indices of 2 exemplars that are closest to each other
  min_val(NN_dist, n_exemplars, &min_d, &min_i);
  min_i2 = NN[min_i];
  min_i = min_i + start_ind;

  while (min_d < thresh) {

    // merge these two exemplars and put into slot of first one.
    // Set slot of second one to all zeros and set count to 0

    // compute mean of merged exemplar
    for (i=0; i<fv_size_half; i++) {
      E[min_i][i] = (E[min_i][i] * (double)count[min_i] + E[min_i2][i] * (double)count[min_i2])/(double)(count[min_i] + count[min_i2]);
      E[min_i][i+fv_size_half] += E[min_i2][i+fv_size_half];
    }

    count[min_i] += count[min_i2];
    count[min_i2] = 0;

    // compute nearest neighbor of merged exemplar
    NN[min_i2-start_ind] = 0;
    NN_dist[min_i2-start_ind] = INF;

    nearest_neighbor_chunk(E, count, min_i, start_ind, end_ind, &indx, &dist);
    NN[min_i-start_ind] = indx;
    NN_dist[min_i-start_ind] = dist;

    // for all entries of NN[i] that equal one of these two indices,
    // recompute the nearest neighbor for that exemplar

    for (i=0; i<n_exemplars; i++)  {
      if (((NN[i] == min_i) || (NN[i] == min_i2)) && (count[start_ind+i] >0))
	nearest_neighbor_chunk(E, count, start_ind+i, start_ind, end_ind, &(NN[i]), &(NN_dist[i]));
    }

    // find indices of 2 exemplars that are closest to each other
    min_val(NN_dist, n_exemplars, &min_d, &min_i);
    min_i2 = NN[min_i];
    min_i += start_ind;
  }

  free(NN);
  free(NN_dist);
}

void exemplar_selection_hierarchical(double *T, int m, int w, double ***E, int *n_ex)
{
  double thresh;
  int n_chunks;
  int chunk_size = 150;
  int *count;
  int mode;
  int i, j, k, nc;
  int *cnt, total_ex, old_total_ex;
  int *start_ind, end_ind;

  // first merge nearby neighbors (same as before)
  thresh = initial_merge(T, m, w, E, &count, n_ex);
  printf("done with initial merge\n");

  // yields matrix E of initial exemplars

  total_ex = (*n_ex);

  // split into chunk_size exemplar chunks
  n_chunks = (int)ceil((double)(*n_ex)/(double)chunk_size);
  start_ind = (int *)malloc((n_chunks+1)*sizeof(int));
  cnt = (int *)malloc((n_chunks+1)*sizeof(int));

  for (i=0; i<n_chunks; i++)
    start_ind[i] = i*chunk_size;
  start_ind[n_chunks] = *n_ex;

  k = n_chunks;

  do {
    n_chunks = k;

    // for each chunk, call exemplar_selection_chunk on range of
    // exemplar indices.  Yields new set of exemplars in same space as
    // input exemplars (with empty slots at the end)

    printf("n_ex = %d, n_chunks = %d\n", (*n_ex), n_chunks);

    old_total_ex = total_ex;
    total_ex = 0;
    cnt[0] = 0;
    for (nc=0; nc<n_chunks; nc++) {

      exemplar_selection_chunk(*E, count, start_ind[nc], start_ind[nc+1], thresh);

      cnt[nc+1] = 0;

      for (i=start_ind[nc]; i<start_ind[nc+1]; i++) {
	if (count[i] > 0)
	  cnt[nc+1]++;
      }

      total_ex += cnt[nc+1];
    }

    // consolidate exemplars (gets rid of empty exemplar slots,
    // signified by count[i] == 0)

    // move all used entries of E to fill up beginning of E and then
    // realloc to get rid of extra entries

    mode = 0;
    for (i=0; i<(*n_ex); i++) {
      if (mode == 0) {
	if (count[i] == 0) {
	  j = i;
	  mode = 1;
	}
      }
      else if (mode == 1) {
	if (count[i] > 0) {
	  // copy exemplar i to j
	  for (k=0; k<fv_size; k++)
	    (*E)[j][k] = (*E)[i][k];
	  count[j] = count[i];
	  count[i] = 0;

	  // find next hole after j
	  while (count[j] > 0) {
	    j++;
	  }
	}
      }
    }

    (*n_ex) = j;

    // copy cnt array to start_ind array, putting neighboring chunks together
    k=0;
    for (i=1; i<n_chunks; i+=2) {
      start_ind[k+1] = start_ind[k] + cnt[i] + cnt[i+1];
      k++;
    }
    if ((n_chunks % 2) == 1) {
      start_ind[k+1] = start_ind[k] + cnt[n_chunks];
      k++;
    }

    (*E) = (double **)realloc((*E), (*n_ex)*sizeof(double *));
    count = (int *)realloc(count, (*n_ex)*sizeof(int));
    printf("total_ex = %d, (*n_ex) = %d\n", total_ex, (*n_ex));
    if (total_ex != (*n_ex)) {
      printf("error, total_ex should be same as n_ex!\n");
      exit(-1);
    }
  } while (n_chunks > 1);  // keep combining chunks until there is only one left


  // Now compute standard deviation components of exemplars -
  // currently the sum of squared values is stored.

  for (j=0; j<(*n_ex); j++) {
    for (i=fv_size_half; i<fv_size; i++) {

      (*E)[j][i] = (*E)[j][i]/(double)count[j] - (*E)[j][i-fv_size_half] * (*E)[j][i-fv_size_half];
      // don't let stddev be too small - avoids divide by zero problem later
      if ((*E)[j][i] < 0.001)
	(*E)[j][i] = .001;
      (*E)[j][i] = sqrt((*E)[j][i]);
    }
  }

  free(start_ind);
  free(cnt);
}
