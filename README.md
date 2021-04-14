# Maximum pseudo-likelihood estimation for superpositions of Strauss-hard core and Poisson processes

This project includes an R package containing an algorithm to compute an maximum pseudo-likelihood estimator for the superposition of a Strauss process (possibly with hard core) and a Poisson process. 

This repository is supplementary to Wieditz, J., Pokern, Y., Schuhmacher, D., Huckemann, S. (2021+). Separating Point Patterns for Fingerprints. Submitted to GSI 2021.

# The SuperMPLE package

To use the SuperMPLE R-package follow the steps below:

1. Install the R-package SuperMPLE via

	`library(remotes)`

	`install_github('jwieditz/SuperMPLE/SuperMPLE')`.
	
2. Load the library via `library(SuperMPLE)`.

3. To compute an MPLE for a homogeneous example, run `example(SuperMPLE)`.

| <img src="/home/jwieditz/SuperMPLE/mple-example-pattern.png" alt="PL example pattern" style="zoom:50%;" /> | <img src="/home/jwieditz/SuperMPLE/mple-example-pl.png" alt="PL example pattern" style="zoom:150%;" /> |
| :----------------------------------------------------------: | ------------------------------------------------------------ |
| *One sampled point pattern of the superposition of a Strauss process with parameters <img src="https://render.githubusercontent.com/render/math?math=(\beta \mu, \gamma, r, R) = (1 \times 42, 0.4, 0.03, 0.1)"> and  process Poisson process of intensity 12 on the unit square.* | *The log pseudo-likelihood on a discrete grid of size 20x20x20 and one computed MPLE <img src="https://render.githubusercontent.com/render/math?math=\hat\theta = (0.29, 1, 32.78)">. The optimisation was started in the point <img src="https://render.githubusercontent.com/render/math?math=\theta_0 = (1, 0.4, 12)">.* |



# Licence

This package is released under the [GPL3.0 licence](https://github.com/jwieditz/SuperMPLE/blob/main/license).

