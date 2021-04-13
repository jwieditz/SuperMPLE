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

# Licence

This package is released under the [GPL3.0 licence](https://github.com/jwieditz/SuperMPLE/blob/main/license).

