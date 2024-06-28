# Linear Preservers of Secant Varieties
This work is joint with [Benjamin Lovitz](https://www.benjaminlovitz.com/) and [Fulvio Gesmundo](https://fulges.github.io/). See the associated paper (link tbd!).

This repository holds an implementation of the lie algebra method of determining the preserver group of a secant variety to a Segre variety. 
We have computed cases where the lie algebra dimension derived from the minor flattenings matched the expected dimension.
The list of computations are given below.

| variety | dimension of monomial vector space |
|--|--|
| $` \hat{\sigma}_2(\mathbb{P}^{1} \times \mathbb{P}^1 \times \mathbb{P}^1 \times \mathbb{P}^1)`$ | 816 |
| $` \hat{\sigma}_3(\mathbb{P}^{1} \times \mathbb{P}^1 \times \mathbb{P}^1 \times \mathbb{P}^1)`$ | 3876 |
| $` \hat{\sigma}_2(\mathbb{P}^{1} \times \mathbb{P}^1 \times \mathbb{P}^1 \times \mathbb{P}^2)`$ | 2600 |
| $` \hat{\sigma}_2(\mathbb{P}^{1} \times \mathbb{P}^1 \times \mathbb{P}^2 \times \mathbb{P}^2)`$ | 8436 |
| $` \hat{\sigma}_3(\mathbb{P}^{1} \times \mathbb{P}^1 \times \mathbb{P}^1 \times \mathbb{P}^2)`$ | 17550 |
| $` \hat{\sigma}_2(\mathbb{P}^{1} \times \mathbb{P}^1 \times \mathbb{P}^1 \times \mathbb{P}^2) \times \mathbb{P}^2)`$ | 5984 |

The reason we give the dimension of the monomial vector space is that the blowup of the vector space size limits the computational feasability.
In the case of $` \hat{\sigma}_3(\mathbb{P}^{1} \times \mathbb{P}^1 \times \mathbb{P}^1 \times \mathbb{P}^1) `$, we must find the rank of a matrix with size $` \approx 2 \times 10^6`$.
