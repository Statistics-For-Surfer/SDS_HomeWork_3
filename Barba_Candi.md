Statistics for Data Science - Homework 3
================
Barba Paolo, Candi Matteo

## Purpose and Statistical tools used

The goal of the project is to perform Friedman’s procedure as a proxy of
multiple studies and applyi it on fMRI (functional magnetic resonance
imaging) data in order to test wheter the data from td (tipacally
develpment) subjects and asd (authism spectrum disorder) subject came
from the same distribution or not. In order to reach this goal we used
statistical tools as machine learning algorithm and two-sample
Hypothesis testing.

## Exercise 3

Friedman’s procedure is used for solving two-sample hypothesis testing
when the each observation consists of many measured attributes
![x\_{i} = x\_{i_1}, x\_{i_2} . . . x\_{i\_{M}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_%7Bi%7D%20%3D%20x_%7Bi_1%7D%2C%20x_%7Bi_2%7D%20.%20.%20.%20x_%7Bi_%7BM%7D%7D "x_{i} = x_{i_1}, x_{i_2} . . . x_{i_{M}}")
and
![z\_{i} = z\_{i_1}, z\_{i_2} . . . z\_{i\_{M}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bi%7D%20%3D%20z_%7Bi_1%7D%2C%20z_%7Bi_2%7D%20.%20.%20.%20z_%7Bi_%7BM%7D%7D "z_{i} = z_{i_1}, z_{i_2} . . . z_{i_{M}}")
where M is the dimensionality of the dataset and we have
![n\_{0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n_%7B0%7D "n_{0}")
samples for
![\\underline{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cunderline%7Bx%7D "\underline{x}")
and
![n\_{1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n_%7B1%7D "n_{1}")
samples for
![\\underline{z}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cunderline%7Bz%7D "\underline{z}")
.

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for example:

![](Barba_Candi_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
