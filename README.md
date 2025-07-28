# Weighted least squares subdivision schemes for noisy data on triangular meshes

This repository includes MATLAB code that generates figures for the paper:

[1] Costanza Conti, Sergio López-Ureña, Dionisio F. Yáñez. Weighted least squares subdivision schemes for noisy data on triangular meshes.  	
https://doi.org/10.48550/arXiv.2507.18976

Please, cite the paper if you use this code in your work.

Abstract: This paper presents and analyses a new family of linear subdivision schemes to refine noisy data given on triangular meshes. The subdivision rules consist of locally fitting and evaluating a weighted least squares approximating first-degree polynomial. This type of rules, applicable to any type of triangular grid, including finite grids or grids containing extraordinary vertices, are geometry-dependent which may result in non-uniform schemes. For these new subdivision schemes, we are able to prove reproduction, approximation order, denoising capabilities and, for some special type of grids, convergence as well. Several numerical experiments demonstrate that their performance is similar to advanced local linear regression methods but their subdivision nature makes them suitable for use within a multiresolution context as well as to deal with noisy geometric data as shown with an example.

The dataset "sphere100" was generated using [spheretri](https://github.com/pgagarinov/spheretri).
The locfit and associated functions (`locfit`, `predict`, `mexlf`, `mexpp`) are from the [locfit](https://github.com/brian-lau/locfit) repository.

Copyright (c) 2025 Costanza Conti, Sergio López-Ureña, Dionisio F. Yáñez-Avedaño.
