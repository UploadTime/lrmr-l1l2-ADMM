# README
____

## Describtion:

    - This is the Matlab code for my B.S. thesis "Low-rank Matrix Recovery Problem Minimizing a Ratio of Two Norms Approximating the Rank Function".

## Algorithms:

### SVT.m:
    - Function of SVT algorithm.
    - Based on: Cai, Jian-Feng and Cand¨¨s, Emmanuel J and Shen, Zuowei, "A singular value thresholding algorithm for matrix completion", SIAM Journal on optimization, 2010.
    - Author: Hao Liang. (https://github.com/HauLiang/Matrix-Completion-Methods.git)

### WNNM.m: 
    - Function of WNNM algorithm.
    - Based on: "Weighted nuclear norm minimization and its applications in low level vision". S. Gu, Q. Xie, D. Meng, W. Zuo, X. Feng, L. Zhang.
    - Code is written based on "The Augmented Lagrange Multiplier Method for Exact Recovery of Corrupted Low-Rank Matrices". Z. Lin, M. Chen, L. Wu. arXiv:1009.5055, 2010.

### Sp_lp_new.m:
    - Function of Schatten p norm combined with $\ell_p$ norm.
    - Based on: Nie, Feiping and Wang, Hua and Huang, Heng and Ding, Chris, "Joint Schatten p-norm and lp-norm robust matrix completion for missing value recovery", Knowledge and Information Systems, 2015.
    - Author: Hao Liang. (https://github.com/HauLiang/Matrix-Completion-Methods.git)

### NF.m:
    - Function of N/F-ADMM alogrithm.
    - Based on: Kaixin Gao, Zheng-Hai Huang, Lulu Guo, "Low-rank matrix recovery problem minimizing a new ratio of two norms approximating the rank function then using an ADMM-type solver with applications", Journal of Computational and Applied Mathematics, Volume 438, 2024.
    - Author: Me.

### l1l2.m:
    - Function of $\ell_1 \ell_2$-ADMM algorithm.
    - Algorithm is proposed by me.
    - Author: Me.

## Other Files:

    psnr.m:
        - Funcition for computing Peak signal-to-noise ratio.

    .mlx:
        - lrmr_re.mlx: Relative error for SVT, WNNM, Sp & lp and N/F.
        - lrmr_re_nf_l1l2.mlx: Relative error for N/F and l1l2.
        - lrmr_re_l1l2_nf_svt_wnnm.mlx: Relative error for SVT, WNNM, N/F and l1l2.
        - image_psnr_re.mlx: Image, psnr and re for two figures with different algorithms.

    figures:
        - house.png
        - tree.jpg

## Folder:

    PROPACK: Package for WNNM.