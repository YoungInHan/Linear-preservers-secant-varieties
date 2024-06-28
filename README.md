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
The runtime of the above ranges from 1 minute in the simplest case to 4 hours for the largest. These measurements were done via a home computer.

## Simple implementation

The simple implementation is a straightforward translation of the method to code. Much of the code is creating a translation 
between monomials in sage and sparse matrices.

## Sparse span trick

This is a trick to allow for computations with larger matrix dimensions by reducing the orthogonal complement computation.
Consider flattenings $`f_1,\ldots,f_k`$ and $`X`$ the action of $`\mathfrak{gl}_N`$ on $`\mathbb{C}[x_1,\ldots,x_N]^d`$.
For finding a basis for the orthogonal complement we can restrict the computation to the subspace $`F:=\text{span}\{f_1,\ldots,f_k,Xf_1,\ldots Xf_k \} \subseteq \mathbb{C}[x_1,\ldots,x_N]^d`$.
In particular, we can do the orthogonal complement computation in $`F`$ and complete the basis in $`\mathbb{C}[x_1,\ldots,x_N]^d`$ by simply adding
standard basis vectors in $`\mathbb{C}[x_1,\ldots,x_N]^d\setminus F`$.
