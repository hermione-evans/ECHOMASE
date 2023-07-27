# Model-agnostic search for the quasinormal modes of gravitational wave echoes
**Di Wu <sup>1,2</sup> Pengyuan Gao <sup>1</sup>, Jing Ren <sup>1</sup>, Niayesh Afshordi <sup>3,4,5</sup>**

<sup>1 Institute of High Energy Physics, Chinese Academy of Sciences, Beijing 100049, China</sup><br />
<sup>2 School of Physics Sciences, University of Chinese Academy of Sciences, Beijing 100039, China</sup><br />
<sup>3 Waterloo Centre for Astrophysics, University of Waterloo, Waterloo, ON, N2L 3G1, Canada</sup><br />
<sup>4 Department of Physics and Astronomy, University of Waterloo, Waterloo, ON, N2L 3G1,
Canada</sup><br />
<sup>5 Perimeter Institute of Theoretical Physics, 31 Caroline St. N., Waterloo, ON, N2L 2Y5, Canada</sup><br />

<!-- This repository is a companion to [Wu et al., arXiv:23xx.xxxxx](https://arxiv.org/abs/).
We include the code to generate all the data presented in the paper and some best-fit value results to reproduce some of our results. We do not provide full posterior files because of the size limit. But we show our compressed posterior with mini-K-means method in `./echo_wave_injection`. -->

This repository is a supplement to our paper, [Wu et al., arXiv:23xx.xxxxx](https://arxiv.org/abs/). It includes the code used to generate all the data presented in the paper, along with the best-fit value results, which you can use to reproduce some of our results. Due to size constraints, we're not able to provide the full posterior files. However, we do present our compressed posterior using the mini-K-means method, which is located in `./echo_wave_injection`.

<!-- This repository struction is organized by our paper. The folder `./comparison_likelihood&SNR_newlikelihood_dependence` contains the code and data to replicate Fig 1 and Fig 2 in Sec 3.2, and Fig 16 in Appendix C. The folder `./uniew_injection` contain the notebook to replicate Fig 3,4,5,6 in Sec 3.3 and Fig 17 Top in Appexdix C, but do not include the posterior to reproduce Fig 5 and Fig 6 middle and bottom part. The folder `./echo_wave_injection` contains the code to replicate Fig 9,10,11,12,13 in Sec 4.2, Fig 17 bottom in Appexdix C, Table 3,4 and Fig 18 in Appendix D. It also include the compressed posterior to reproduce all the results. The folder `./superposition_QNMS` contains the data and code to replicate Fig 14 and Fig 15 in Appendix B.

We encourage use of these data in derivative works. If you use the material provided here, please cite the paper using the reference: -->

## Repository Structure
The repository is organized to mirror the structure of our paper:

1. `./comparison_likelihood&SNR_newlikelihood_dependence`: This directory contains the code and data used to make a simple comparison of the two phase-marginalized likelihoods in Sec. 3.2 (Figures 1 and 2), and to plot the new likelihood dependence of the echo SNR in Appendix C (Figure 16).

1. `./uniew_injection`: This directory contains the notebook that showcases the Bayesian search results for the UniEw injections in different realizations of Gaussian noise. Specifically, it can be used to reproduce Figures 3, 4, 5, 6 in Section 3.3 and Figure 17 (Top) in Appendix D. Note that the posterior files for reproducing the middle and bottom parts of Figures 4 and 5 are not included.

1. `./echo_wave_injection`: This directory contains the notebook that showcases the Bayesian search results for injections of the four benchmarks for echo waveforms. Specifically, it can be used to reproduce Figures 9, 10, 11, 12, 13 in Section 4.2, Figure 17 (Bottom) in Appendix C, Table 3, 4 and Figure 18 from Appendix D. It also includes the compressed posterior files needed to reproduce all the results.

1. `./superposition_QNMS`: This directory contains the data and code to reproduce Figures 14 and 15 in Appendix B that examine the influence of QNMs interference.

We invite you to explore these directories and make use of the resources within them to better understand our research and methods.

## Using this Repository
To use this repository, clone it to your local machine and follow the individual instructions within each directory to reproduce the respective results from the paper.

We encourage the use of these data in derivative works. If you make use of the material provided in this repository, we kindly ask you to cite our paper using the following reference:

```
@article{,
    author = "Di Wu, Pengyuan Gao , Jing Ren , Niayesh Afshordi",
    title = "{Model-agnostic search for the quasinormal modes of gravitational wave echoes}",
    eprint = "",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "",
    year = "2023"
}
```

## License
This project is licensed under the terms of the MIT license. For more information, please see the LICENSE file.
