// Copyright (C) 2014, 2023 Mitsubishi Electric Research Laboratories (MERL)
//
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
//#include <iostream>

#define INF 1e20       //Pseudo Infitinte number for this code

#define MIN(a, b) (a < b) ? a : b
#define MAX(a, b) (a > b) ? a : b

//using namespace std;

cmpdoublep(const void *p1, const void *p2)
       {
           if ((*(double * const)p1) < (*(double * const)p2))
               return 1;
           else if ((*(double * const)p1) > (*(double * const)p2))
               return -1;
           else
               return 0;
       }

double find_threshold_yielding_desired_fprate(int *gt, double *anom_score, int w, int n, double fp_rate)
{
    int i,j;
    int w_half;
    int left, right;
    int tooclose;
    double thresh;
    int cnt=0;
    double *as;

    w_half = (int)((double)w/2.0+0.5);

    as = (double *)malloc(n*sizeof(double));

    // for each element of gt that is 0 and at least w/2 distance from
    // nearest 1 element, put anomaly score into new array as
    for (i=0; i<n; i++) {
        if (gt[i] == 0) {
            // see if any of values within w/2 of i are 1
            left = MAX(0, i-w_half);
            right = MIN(n-1, i+w_half);
            tooclose = 0;
            for (j=left; j<=right; j++)
                if (gt[j] == 1)
                    tooclose = 1;

            if (! tooclose) {
                as[cnt++] = anom_score[i];
            }
        }
    }

    printf("number of normal windows = %d\n", cnt);

    // now sort in descending order as and find value of (fp_rate*cnt)th element
    qsort(as, cnt, sizeof(double), cmpdoublep);

    i = (int)(fp_rate*(double)cnt + 0.5);
    printf("i = %d\n", i);

    thresh = as[i];

    free(as);

    printf("thresh = %lf\n", thresh);

    return thresh;
}


int main(  int argc , char *argv[] )
{
    FILE *fp;              // the input file pointer
    FILE *qp;              // the query file pointer

    int i, j;
    int n;
    int cnt;
    int d;
    int *gt;
    double *anom_score;
    int w;
    double thresh;
    double f;
    int mode;
    int n_detections;
    int n_anomalies;
    int hits;
    double fp_rate;

// read in ground truth series for the testing time series
// read in anomaly score time series (computed by BF_AnomalyDet_X)
// find threshold yielding desired fp rate
// compute detection rate given that threshold (detection rate is per block of gt 1's)

    if (argc != 5) {
        printf("Usage: compute_detection_rate gt.txt anomaly_scores.txt window_length fp_rate\n");
        exit(-1);
    }

    w = atoi(argv[3]);
    fp_rate = atof(argv[4]);

    if ((fp = fopen(argv[1],"r")) == NULL) {
        printf("Cannot open file %s to read\n", argv[1]);
        exit(-1);
    }

    // first count how long the gt vector is
    cnt = 0;
    n = fscanf(fp, "%d", &d);
    while (n != EOF) {
        n = fscanf(fp, "%d", &d);
        cnt++;
    }
    fclose(fp);

    printf("time series length = %d\n", cnt);

    // next store gt vector
    if ((fp = fopen(argv[1],"r")) == NULL) {
        printf("Cannot open file %s to read\n", argv[1]);
        exit(-1);
    }

    gt = (int *)malloc(cnt * sizeof(int));

    for (i=0; i<cnt; i++) {
        n = fscanf(fp, "%d", &d);
        gt[i] = d;
    }
    fclose(fp);

    anom_score = (double *)malloc(cnt * sizeof(double));

    if ((fp = fopen(argv[2],"r")) == NULL) {
        printf("Cannot open file %s to read\n", argv[2]);
        exit(-1);
    }

    for (i=0; i<cnt; i++) {
        n = fscanf(fp, "%lf", &f);
        anom_score[i] = f;
    }
    fclose(fp);

    // find threshold that yields 0 false positives
    thresh = find_threshold_yielding_desired_fprate(gt, anom_score, w, cnt, fp_rate);

    // for each block of 1's in gt, count the number of above-threshold
    // anomaly scores
    mode = 0;
    n_detections = 0;
    n_anomalies = 0;
    for (i=0; i<cnt; i++) {
        if (mode == 0) {
            if (gt[i] == 1) {
                mode = 1;
                n_anomalies++;
                hits = 0;
                if (anom_score[i] > thresh)
                    hits++;
            }
        }
        else {  // mode == 1
            if (gt[i] == 1) {
                if (anom_score[i] > thresh)
                    hits++;
            }
            else {
                mode = 0;
                printf("hits = %d\n", hits);
                if (hits > 0)
                    n_detections++;
            }
        }
    }
    if (mode == 1) {
        printf("hits = %d\n", hits);
        if (hits > 0)
            n_detections++;
    }
    printf("n_anomalies = %d, n_detections = %d\n", n_anomalies, n_detections);

}
