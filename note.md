We aim to minimize the zero point of the following equation
```math
1 = \sum_{i = 1}^{s} \gamma^{ - n_i}\;,
```
where $s$ is the number of branches, and $n_i$ is the number of covered elements in the $i$-th branch.
Generally speaking, we want larger $n$ and smaller $s$.

The problem can be encoded as a clustering problem, where the evidences are treated as $l$ dim data points and the distance is the hamming distance between the evidences.
Then $s$ is the number of clusters, and $d_i = l - n_i$ is given by the maximum distance between the elements in the $i$-th cluster.