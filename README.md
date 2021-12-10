# Overview

Downloads the [ESA Earth System Model](https://isdc.gfz-potsdam.de/esmdata/esaesm/) data and converts it to matlab.

Each data file (containing Spherical Harmonic coefficients) is converted to a matlab data file. In these files, there's one single variable:

```
mod =

  struct with fields:

    mod: [16471x4 double]
     GM: 3.9860e+14
      R: 6378137
```
The `GM` and `R` field report Earth's gravitational constant and radius, respectively.
The `mod` field contains 4 columns: SH degree, SH order, cosine coefficient and sine coefficient, e.g.:

```
>> format long
>> mod.mod(1:10,:)

ans =

                   0                   0  -0.000000000161801                   0
   1.000000000000000                   0   0.000000000021131                   0
   1.000000000000000   1.000000000000000  -0.000000000027007   0.000000000093393
   2.000000000000000                   0  -0.000000000041471                   0
   2.000000000000000   1.000000000000000  -0.000000000039424  -0.000000000018933
   2.000000000000000   2.000000000000000  -0.000000000025483   0.000000000022824
   3.000000000000000                   0  -0.000000000120588                   0
   3.000000000000000   1.000000000000000  -0.000000000101073  -0.000000000046886
   3.000000000000000   2.000000000000000  -0.000000000072408  -0.000000000023873
   3.000000000000000   3.000000000000000   0.000000000036481   0.000000000001653
```

# How to rebuild this database

1. `download.sh`
2. `extract.sh`
3. `convert.m`