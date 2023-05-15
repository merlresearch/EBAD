// Copyright (C) 2014, 2023 Mitsubishi Electric Research Laboratories (MERL)
//
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIGN(a) (a < 0) ? -1 : (a > 0) ? 1 : 0
#define STRIDE 2

extern int fv_size;
extern int fv_size_half;

void compute_SST_feature(double *T, int indx, int w, double *fv, int m)
{

    int i, j, left, right;
    int lngth;
    double mu, stddev;
    double mu2;
    double mad;
    int npd, nzd, nmc;
    int num_runs, tot_run_len, run;
    int sign_d1;
    double d1, sm;

    //    int fv_size = 2*((int)(w/STRIDE + 0.5) + 7);
    //    int fv_size_half = fv_size/2;

    for (i=0; i<fv_size_half; i++) {
        fv[i] = 0.0;
    }

    // the first n components are the smoothed trajectory,
    //  the next 7 are the statistics of the window, and the remaining
    //  fv_size_half+7 components are the standard deviations of the
    //  first fv_size_half+7 components (initialized to 0 until
    //  clustering is done)

    // first compute smoothed trajectory (low pass filter of z minus mean)
//    m = (int)((w/25)+0.5);
    mu = 0.0;
    stddev = 0.0;

    for (i=indx; i<indx+w; i++) {
        mu += T[i];
        stddev += T[i]*T[i];
    }
    mu /= (double)w;
    stddev = sqrt(stddev/(double)w - mu*mu);

    if (stddev < .001) {
        stddev = .001;
    }

    // the following is much faster code for computing smoothed time series
    // that doesn't compute redundant parts of sums

//    mu2 = 0.0;
    left = indx-m;
    right = indx+m;
    sm = 0.0;
    for (i=left; i<=right; i++)
        sm += T[i];

    lngth = right-left+1;

    for (i=0; i<w-1; i++) {
        if ((i % STRIDE) == 1)
            fv[(i-1)/STRIDE] = sm/lngth - mu;

//        mu2 += fv[i];

        sm -= T[left];
        left++;

        right++;
        sm += T[right];
    }

    if (((w-1) % STRIDE) == 1)
        fv[(w-2)/STRIDE] = sm/lngth - mu;

//    mu2 += fv[w-1];
//    mu2 /= (double)w;

//    for (i=0; i<w; i++)
//        fv[i] -= mu2;

    // next compute the statistical components

    // mean
    fv[fv_size_half - 7] = mu;

    // std. dev.
    fv[fv_size_half - 6] = stddev;

    // mean abs difference
    mad = 0.0;
    nmc = 0;
    npd = 0;
    nzd = 0;
    run = 0;
    tot_run_len = 0;
    num_runs = 0;

    for (i=indx; i<indx+w; i++) {
        if (i < indx+w-1) {
            d1 = T[i+1] - T[i];
            mad += fabs(d1);

            if (d1 > 0.0)
                npd += 1;

            if (d1 == 0.0)
                nzd += 1;
//            printf("i=%d, d1=%lf, mad=%lf, npd=%d, nzd=%d\n", i, d1, mad, npd, nzd);

            if (((T[i+1] - mu) > 0.0) && ((T[i] - mu) < 0.0))
                nmc += 1;
            else if (((T[i+1] - mu) < 0.0) && ((T[i] - mu) > 0.0))
                nmc += 1;

            sign_d1 = SIGN(d1);

            if (run == 1)
                if (sign_d1 > 0)
                    tot_run_len++;
                else
                    run = 0;

            else if (sign_d1 > 0) {
                num_runs++;
                tot_run_len++;
                run = 1;
            }
        }
    }

    fv[fv_size_half - 5] = mad/(double)(w-1);
    fv[fv_size_half - 4] = (double)nmc/(double)(w-1);
    fv[fv_size_half - 3] = (double)npd/(double)(w-1);
    fv[fv_size_half - 2] = (double)nzd/(double)(w-1);
    if (num_runs == 0) {
        num_runs = 1;
        tot_run_len = 1;
    }

    fv[fv_size_half - 1] = (double)tot_run_len/(double)(num_runs*(w-1));
}
