# Tridiagonal matrix with corners

system of equations:

| 4 1 0 0 0 0 1 |   |x1|   |1|   
| 1 4 1 0 0 0 0 |   |x2|   |2|
| 0 1 4 1 0 0 0 |   |x3|   |3|
| 0 0 1 4 1 0 0 | * |x4| = |4|
| 0 0 0 1 4 1 0 |   |x5|   |5|
| 0 0 0 0 1 4 1 |   |x6|   |6|
| 1 0 0 0 0 1 4 |   |x7|   |7|

I this solved problem for SIZE = 7 but you can change it

My solution uses Sherman - Morisson formula and Givens rotation.
It is very important, reduces complexity from O(N^3) to O(N)

Detailed descritpion of my program is in polish language in the comments.

This is one of the programs written by me from numerical methods.

## Program results
 
Vector Z
 z[1] = 0.2281
 z[2] = 0.315699
 z[3] = 0.509103
 z[4] = 0.647887
 z[5] = 0.899347
 z[6] = 0.754723
 z[7] = 2.08176
 
Vector Q
 q[1] = 0.366197
 q[2] = -0.0985915
 q[3] = 0.028169
 q[4] = -0.0140845
 q[5] = 0.028169
 q[6] = -0.0985915
 q[7] = 0.366197
 
 (v^T)*Z / 1+(v^T)*q = 1.33333
 
Vector X
 x[1] = -0.260163
 x[2] = 0.447154
 x[3] = 0.471545
 x[4] = 0.666667
 x[5] = 0.861789
 x[6] = 0.886179
 x[7] = 1.5935
 

## Authors

* **Bartłomiej Orawiec** - (https://github.com/OrawiecBartlomiej)
