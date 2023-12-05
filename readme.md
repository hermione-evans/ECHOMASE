# Model-independent search for the quasinormal modes of gravitational wave echoes
**Di Wu <sup>1,2</sup> Pengyuan Gao <sup>1</sup>, Jing Ren <sup>1</sup>, Niayesh Afshordi <sup>3,4,5</sup>**

<sup>1 Institute of High Energy Physics, Chinese Academy of Sciences, Beijing 100049, China</sup><br />
<sup>2 School of Physics Sciences, University of Chinese Academy of Sciences, Beijing 100039, China</sup><br />
<sup>3 Waterloo Centre for Astrophysics, University of Waterloo, Waterloo, ON, N2L 3G1, Canada</sup><br />
<sup>4 Department of Physics and Astronomy, University of Waterloo, Waterloo, ON, N2L 3G1,
Canada</sup><br />
<sup>5 Perimeter Institute of Theoretical Physics, 31 Caroline St. N., Waterloo, ON, N2L 2Y5, Canada</sup><br />

This repository is a supplement to our paper, [Wu et al., arXiv:2308.01017](https://arxiv.org/abs/2308.01017). The link of published version is https://journals.aps.org/prd/abstract/10.1103/PhysRevD.108.124006. It includes the code used to generate all the data presented in the paper, along with the best-fit value results, which you can use to reproduce some of our results. Due to size constraints, we're not able to provide the full posterior files. However, we do present our compressed posterior using the mini-K-means method, which is located in `./echo_wave_injection`.

## Repository Structure
The repository is organized to mirror the structure of our paper:

1. `./comparison_likelihood&SNR_newlikelihood_dependence`: This directory contains the code and data used to make a simple comparison of the two phase-marginalized likelihoods in Sec. III.B (Figures 1 and 2), and to plot the new likelihood dependence of the echo SNR in Appendix C (Figure 17).

2. `./uniew_injection`: This directory contains the notebook that showcases the Bayesian search results for the UniEw injections in different realizations of Gaussian noise. Specifically, it can be used to reproduce Figures 3, 4, 5, 6 in Section III.C and Figure 18 (Top) in Appendix D. Note that the posterior files for reproducing the middle and bottom parts of Figures 4 and 5 are not included.

3. `./echo_wave_generate`: This directory contains the notebook to generate the injections of the four benchmarks for echo waveforms. Specifically, it can be used to reproduce Figures 7, 8 in Section IV.A.

4. `./echo_wave_injection`: This directory contains the notebook that showcases the Bayesian search results for injections of the four benchmarks for echo waveforms. Specifically, it can be used to reproduce Figures 9, 10, 11, 12, 13 in Section IV.B, Figure 14 in Appendix B,Figure 18 (Bottom) in Appendix C, Table 3, 4 and Figure 19 from Appendix D. It also includes the compressed posterior files needed to reproduce all the results.

5. `./superposition_QNMS`: This directory contains the data and code to reproduce Figures 15 and 16 in Appendix B that examine the influence of QNMs interference.

We invite you to explore these directories and make use of the resources within them to better understand our research and methods.

## Using this Repository
To use this repository, clone it to your local machine and follow the individual instructions within each directory to reproduce the respective results from the paper.

We encourage the use of these data in derivative works. If you make use of the material provided in this repository, we kindly ask you to cite our paper using the following reference:

```
@article{Wu:2023wfv,
    author = "Wu, Di and Gao, Pengyuan and Ren, Jing and Afshordi, Niayesh",
    title = "{Model-independent search for the quasinormal modes of gravitational wave echoes}",
    eprint = "2308.01017",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.108.124006",
    journal = "Phys. Rev. D",
    volume = "108",
    number = "12",
    pages = "124006",
    year = "2023"
}
```

## License
This project is licensed under the terms of the MIT license. For more information, please see the LICENSE file.
