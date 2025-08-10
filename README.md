# Multiband-Multitaper Nonuniform Fast Fourier Transform (M<sup>2</sup>NuFFT)

> A computationally efficient suboptimal power spectrum estimator for fast exploration of nonuniformly sampled time series

## Introduction

This is the code for the paper pre-print [(Cui 2024)](https://arxiv.org/abs/2407.01943).

## Getting Started

1. Download and install [Chronux](https://github.com/jiecui/chronux) computational toolbox. Please use this fork of `Chronux` as some of the original codes need to be modified for compatibility.

1. Download and install [M<sup>2</sup>NuFFT](https://github.com/jiecui/mtnufft) package.

## Build and Test

1. Error analysis of MTNUFFT method

   ```mtnufft_error_analysis.m ```

1. Speed analysis of MTNUFFT method

   ```mtnufft_speed_analysis.m```

1. Analysis of example impedance signal

   ```imp_example_analysis.m```

## Contribute

 <img src="./images/mayo_logo.png" alt="[Mayo Clinic]" width="50"/> [Laboratory of Bioelectronics Neurophysiology and Engineering at Mayo Clinic](https://www.mayo.edu/research/labs/bioelectronics-neurophysiology-engineering/overview)

## References

* __J. Cui__, B. H. Brinkmann, G. A. Worrell, _A fast multitaper power spectrum estimation in nonuniformly sampled time series_, [arXiv, 5704101](https://arxiv.org/abs/2407.01943), 2024 [[PDF]](./docs/fast_mtnufft_arxiv_v1.pdf) (DSP under revision).
