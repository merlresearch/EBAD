<!--
Copyright (C) 2014, 2023 Mitsubishi Electric Research Laboratories (MERL)

SPDX-License-Identifier: AGPL-3.0-or-later
-->

# EBAD Exemplar-Based Anomaly Detection

## Features

This code implements the hierarchical exemplar learning, SST
exemplar-based anomaly detection and brute-force Euclidean distance
anomaly detection algorithms described in the paper "Exemplar Learning
for Extremely Efficient Anomaly Detection in Real-Valued Time Series"
by Michael Jones, Daniel Nikovski, Makoto Imamura and Takahisa Hirata
appearing in the Journal of Data Mining and Knowledge Discovery 2015.

## Installation

The code has been tested on Ubuntu 14.04 using the cc compiler.  It
has not been tested on Windows, although it should not be difficult to
get it to compile there.

There is a simple Makefile for each different program:

* `make_MERL_model_efficient`
* `make_BF_AnomalyDet_MERL_SST`
* `make_compute_detection_rate`
* `make_BF_AnomalyDet_UCR_ED`

To compile a particular program, simply run make on the approriate makefile.  For example:

```bash
make -f make_MERL_model_efficient
```

## Usage

To train an SST exemplar model on a training time series (not
containing any anomalies):

```bash
MERL_model_efficient <training time series> <number of time steps in training time series> <window length> <output model filename>
```

example:

```bash
MERL_model_efficient arma_train.txt 10000 100 arma_model.txt
```

To compute anomaly scores on a testing time series:

```bash
BF_AnomalyDet_MERL_SST <model filename> <testing time series> <number of time steps in testing time series> <window length>
```

example:

```bash
BF_AnomalyDet_MERL_SST arma_model.txt arma_test.txt 100000 100
```

This outputs a text file named arma_test.out containing the anomaly
scores for each window (windows overlap and have a step size of 1).
This output file can be read into Matlab, for example, and plotted.


To compute the detection rate given a ground truth file (of 0's and 1's) and an output file:

```bash
compute_detection_rate <ground truth testing time series> <anomaly score time series> <window length> <desired false positive rate from 0 to 1>
```

example:

```bash
compute_detection_rate arma_test_gt.txt arma_test.out 100 0
```

## Citation

If you use the software, please cite the following  ([TR2014-042](https://merl.com/publications/TR2014-042)):

```BibTeX
@inproceedings{Jones2014jun,
    author = {Jones, M. and Nikovski, D. and Imamura, M. and Hirata, T.},
    title = {Anomaly Detection in Real-valued Multidimensional Time Series},
    booktitle = {ASE Bigdata/Socialcom/Cyber Security Conference},
    year = 2014,
    month = jun,
    publisher = {Academy of Science and Engineering (ASE)},
    isbn = {978-1-62561-000-3},
    url = {https://www.merl.com/publications/TR2014-042}
}
```

## Contact

Mike Jones <mjones@merl.com>

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for our policy on contributions.

## License

Released under `AGPL-3.0-or-later` license, as found in the [LICENSE.md](LICENSE.md) file.

All files:

```
Copyright (C) 2014, 2023 Mitsubishi Electric Research Laboratories (MERL).

SPDX-License-Identifier: AGPL-3.0-or-later
```
