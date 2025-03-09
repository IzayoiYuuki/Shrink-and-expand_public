# Welcome to shrink-and-expand technique project

## Introduction

The **shrink-and-expand (SE) technique** is a method designed to accelerate eigensolvers.
This technique can be widely applied to a broad class of eigensolvers, including well-known algorithms such as **subspace iteration**, the **steepest descent method**, and **LOBPCG**.

This Git project is based on our paper *'On a shrink-and-expand technique for block eigensolvers（https://arxiv.org/abs/2409.05572）'* and includes simple implementations of the **shrink-and-expand (SE) technique** in five algorithms: **subspace iteration**, the **steepest descent method**, **LOBPCG**, **trace minimization**, and the **Chebyshev-Davidson method**.

## How to test

The project includes tests of the **shrink-and-expand (SE) technique** with three strategies (**fix**, **slope**, **slopek**) across five algorithms. Note that the **Chebyshev-Davidson method** does not have implementations for the **slope** and **slopek** strategies.  

To quickly evaluate the performance of the three strategies on the five algorithms, you can run **`test_all.m`** in the root directory. By default, only the first three matrices will be tested.  

If you wish to test individual algorithms and strategies separately, you can run the corresponding script files. The implementations of the five algorithms and their test examples are organized into five folders under the **`./Eigensolver`** directory. You can test specific cases by running the **`test_<algorithm>_<strategy>.m`** files.  

To test all 12 examples, you can modify the **`fileNo`** variable in the scripts. However, please be aware that some matrices are quite large, which may result in longer testing times.