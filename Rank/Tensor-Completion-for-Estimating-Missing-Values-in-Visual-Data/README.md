# How to run MATLAB scripts?

## 1. Create a dirty image
```shell
init_image; parameters initialization 
```
Functions: 

```
Replacelnd.m;  Corrupting Image, Missing

FindInd.m; Missing index

Trimming.m; Matrix Decomposition

Mfold.m; Imputation

ten2mat.m; tensor operation

shrinkage.m; SVD
```

After that, make sure you're already initialized `Image` and `known` var


## 2. Run scripts


```shell
epoch = 50;
FaLRTC(Image, epoch, known); Imputation Algorithm
#Time Consumption prediction: 100s

HaLRTC(Image, epoch, known); Imputation Algorithm
#Time Consumption prediction: 200s

SiLRTC(Image, epoch*10, known); Imputation Algorithm
#Time Consumption prediction: 1000s
```
