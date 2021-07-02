#' Pre-computes graph structures of the point pattern and the integration grid to speed up the PL optimisation.
#'
#' @param pointPattern the considered point pattern.
#' @param hardCoreRadius the hard core distance $r>0$ of the Strauss process.
#' @param interactionRadius the interaction radius $R > r$ of the Strauss process.
#' @param trend the trend/ activation function of the Strauss process. Might be a real number or a pixel image.
#' @param gridSizeX the number of grid points to compute the normalisation constant of the pseudo-likelihood in x-direction.
#' @param gridSizeY the number of grid points to compute the normalisation constant of the pseudo-likelihood in y-direction.
#'
#' @return A list of containing the pre-computed graph structures of the point pattern and the integration grid.
#' @author Johannes Wieditz
#' @export
#'
compute_graphStructure <- function( pointPattern, interactionRadius, hardCoreRadius, trend, gridSizeX, gridSizeY, isHomogeneous ){

  graphStructure    <- list()

  message('Computing graph structure of the point pattern')
  pb = txtProgressBar(min = 0, max = pointPattern$n, initial = 0)
  for( i in 1:pointPattern$n ){

    setTxtProgressBar(pb, i)

    z               <- subset(pointPattern, i)
    hasNeighbours   <- FALSE
    informationList <- list()

    nodes           <- pointPattern
    distance        <- pairdist(nodes)
    noOfNeighbours  <- 0

    # If z does not interact with any point of the point pattern, the conditional intensity equals the intensity.
    if( any((distance < interactionRadius & distance > 0)[i, ]) ){

      hasNeighbours <- TRUE
      noOfNeighbours<- sum((distance < interactionRadius & distance > 0)[i, ])

      zetaWithoutZ  <- subset(pointPattern, (1:pointPattern$n)[-i])
      distance      <- pairdist(zetaWithoutZ)

      # Compute the interaction graph induced by the point pattern.
      adjacencyMatrix     <- (distance < interactionRadius & distance > 0)
      interactionGraph    <- graph_from_adjacency_matrix(adjacencyMatrix,
                                                         mode = "undirected")

      adjacencyMatrixHC   <- (distance < hardCoreRadius & distance > 0)
      hardcoreGraph       <- graph_from_adjacency_matrix(adjacencyMatrixHC,
                                                         mode = "undirected")

      # Label the vertices accordint to their order in the point pattern.
      vertex_attr(interactionGraph, "label") <- 1:vcount(interactionGraph)
      vertex_attr(hardcoreGraph, "label")    <- 1:vcount(hardcoreGraph)

      # Compute all the components and label the points according to their corresponding component.
      comps               <- components(interactionGraph)
      componentMembers    <- comps$membership
      numberOfComponents  <- comps$no

      for( k in 1:numberOfComponents ){

        # Compute the induced subgraph of the k-th component.
        hasNeighboursInComponent <- FALSE
        vCount        <- c()
        eCount        <- c()
        neighboursOfZ <- c()
        withinHardcore  <- c()
        trendProduct    <- c()
        indices       <- which(componentMembers == k)


        pointsWithZ   <- superimpose(subset(zetaWithoutZ, indices), z)
        distanceWithZ <- pairdist(pointsWithZ)

        # If z interacts with some point in component
        if( sum((distanceWithZ < interactionRadius & distanceWithZ > 0)[pointsWithZ$n,]) > 0 ){
          hasNeighboursInComponent              <- TRUE
          componentGraph                        <- induced_subgraph(interactionGraph, indices)
          vertex_attr(componentGraph, "label")  <- indices

          powSet <- powerSet(1:length(indices))

          # For all subgraphs of this component compute the sum of the numerator and the denominator of the conditional intensity. Note, that they only differ in the factor gamma^number of neighbours of z in the subgraph.
          for( j in powSet ){

            vCount <- c( vCount, vcount(induced_subgraph(componentGraph, j)))
            eCount <- c( eCount, ecount(induced_subgraph(componentGraph, j)))

            # Compute the number of neighbours of z in the considered subgraph.
            pointsWithZ   <- superimpose(subset(zetaWithoutZ, indices[j]), z)
            distanceWithZ <- pairdist(pointsWithZ)
            neighboursOfZ <- c(neighboursOfZ, sum((distanceWithZ < interactionRadius & distanceWithZ > 0)[pointsWithZ$n,]))
            withinHardcore  <- c(withinHardcore, sum((distanceWithZ < hardCoreRadius & distanceWithZ > 0)[pointsWithZ$n,]))

            if( isHomogeneous ){
              trendProduct <- c( trendProduct, 1)
            } else {
              trendProduct <- c( trendProduct, prod(trend[coords(subset(zetaWithoutZ, indices[j]))]))
            }
          }
        }

        informationList <- c(informationList, list( list(component = k, neighboursInComponent = hasNeighboursInComponent, noVertices = vCount, noEdges = eCount, noNeighbours = neighboursOfZ, trendProduct = trendProduct, withinHardcore = withinHardcore)))
      }
    }

    graphStructure <- c(graphStructure, list(list(point = z, interacting = hasNeighbours, noOfNeighbours = noOfNeighbours, info = informationList)))
  }

  integralGrid          <- gridcentres(Window(pointPattern), gridSizeX, gridSizeY)
  integralGridPattern   <- ppp(integralGrid$x, integralGrid$y, window = Window(pointPattern))
  integralGridStructure <- list()

  message('Computing graph structure for numerical integration')
  pb = txtProgressBar(min = 0, max = integralGridPattern$n, initial = 0)

  for( i in 1:(integralGridPattern$n) ){

    setTxtProgressBar(pb, i)

    u               <- integralGridPattern[i]
    isInWindow      <- inside.owin(u$x, u$y, Window(pointPattern))
    hasNeighbours   <- FALSE
    informationList <- list()

    zetaPlusU       <- superimpose(pointPattern, integralGridPattern[i])
    distance        <- pairdist(zetaPlusU)
    noOfNeighbours  <- 0

    if( isInWindow ){

      # If u does not interact with any point of the point pattern, the conditional intensity equals the intensity. Hence, no further computations are needed and the following procedure is skipped.
      if( any((distance < interactionRadius & distance > 0)[pointPattern$n + 1, ]) ){

        hasNeighbours <- TRUE
        noOfNeighbours<- sum((distance < interactionRadius & distance > 0)[pointPattern$n + 1, ])
        zeta          <- pointPattern
        distance      <- pairdist(zeta)

        # Compute the interaction graph induced by the point pattern.
        adjacencyMatrix     <- (distance < interactionRadius & distance > 0)
        interactionGraph    <- graph_from_adjacency_matrix(adjacencyMatrix,
                                                           mode = "undirected")

        adjacencyMatrixHC   <- (distance < hardCoreRadius & distance > 0)
        hardcoreGraph       <- graph_from_adjacency_matrix(adjacencyMatrixHC,
                                                           mode = "undirected")

        # Label the vertices accordint to their order in the point pattern.
        vertex_attr(interactionGraph, "label") <- 1:vcount(interactionGraph)
        vertex_attr(hardcoreGraph, "label")    <- 1:vcount(hardcoreGraph)

        # Compute all the components and label the points according to their corresponding component.
        comps               <- components(interactionGraph)
        componentMembers    <- comps$membership
        numberOfComponents  <- comps$no

        for( k in 1:numberOfComponents ){

          # Compute the induced subgraph of the k-th component.
          hasNeighboursInComponent <- FALSE
          vCount        <- c()
          eCount        <- c()
          neighboursOfU <- c()
          withinHardcore<- c()
          trendProduct  <- c()
          indices       <- which(componentMembers == k)

          pointsWithU   <- superimpose(subset(zeta, indices), u)
          distanceWithU <- pairdist(pointsWithU)

          # If z interacts with some point in component
          if( sum((distanceWithU < interactionRadius & distanceWithU > 0)[pointsWithU$n,]) > 0 ){
            hasNeighboursInComponent              <- TRUE
            componentGraph                        <- induced_subgraph(interactionGraph, indices)
            vertex_attr(componentGraph, "label")  <- indices

            powSet <- powerSet(1:length(indices))

            # For all subgraphs of this component compute the sum of the numerator and the denominator of the conditional intensity. Note, that they only differ in the factor gamma^number of neighbours of z in the subgraph.
            for( j in powSet ){

              vCount <- c( vCount, vcount(induced_subgraph(componentGraph, j)))
              eCount <- c( eCount, ecount(induced_subgraph(componentGraph, j)))

              # Compute the number of neighbours of z in the considered subgraph.
              pointsWithU   <- superimpose(subset(zetaPlusU, indices[j]), u)
              distanceWithU <- pairdist(pointsWithU)
              neighboursOfU <- c(neighboursOfU, sum((distanceWithU < interactionRadius & distanceWithU > 0)[pointsWithU$n,]))
              withinHardcore<- c(withinHardcore, sum((distanceWithU < hardCoreRadius & distanceWithU > 0)[pointsWithU$n,]))

              if( isHomogeneous ){
                trendProduct <- c(trendProduct, trend)
              } else {
                trendProduct  <- c(trendProduct, prod(trend[coords(subset(zetaPlusU, indices[j]))]))
              }
            }
          }

          informationList <- c(informationList, list( list(component = k, neighboursInComponent = hasNeighboursInComponent, noVertices = vCount, noEdges = eCount, noNeighbours = neighboursOfU, trendProduct = trendProduct, withinHardcore = withinHardcore)))
        }
      }

      integralGridStructure <- c(integralGridStructure, list(list(point = u, isInWindow = TRUE, interacting = hasNeighbours, noOfNeighbours = noOfNeighbours, info = informationList)))


    } else {

      integralGridStructure <- c(integralGridStructure, list(list(point = u, isInWindow = FALSE, interacting = hasNeighbours, noOfNeighbours = noOfNeighbours, info = informationList)))
    }
  }

  return( list(graphStructure, integralGridStructure) )
}

#' Computes the conditional intensity of a given point pattern stemming from a PoSt-hard core process.
#'
#' @param theta the model parameters.
#' @param graphStructure the interaction graph induced by the point pattern as computed by the function compute_graphStructure (1st component).
#' @param i the index of the point of the point pattern considered.
#' @param trend the trend/ activation function of the Strauss process. Might be a real number or a pixel image.
#' @param isHomogeneous a flag indicating whether the considered process is homogeneous (TRUE) or not (FALSE)
#'
#' @return The conditional intensity of the i-th point in the point pattern listed in the graph structure passed.
#' @author Johannes Wieditz
#' @export
#'

conditionalIntensity <- function( theta, graphStructure, i, trend, isHomogeneous ){

  theta   <- exp(theta)
  beta    <- theta[1]
  gamma   <- theta[2]
  lambda  <- theta[3]

  if( isHomogeneous ) {

    betaZ <- beta * trend
  } else {

    betaZ <- beta * interp.im(trend, graphStructure[[i]]$point$x, graphStructure[[i]]$point$y)
  }

  if( !graphStructure[[i]]$interacting ){

    return( lambda + betaZ )

  } else {

    components <- graphStructure[[i]]$info

    if( lambda == 0 ){

      return( betaZ * gamma^(graphStructure[[i]]$noOfNeighbours) )
    }

    if( isHomogeneous ){
      relativeTrend       <- beta * trend / lambda
    } else{
      relativeTrend       <- beta / lambda
    }
    interactionProduct  <- 1

    for( j in 1:length(components) ){

      if( components[[j]]$neighboursInComponent ){

        withinHardcore  <- (components[[j]]$withinHardcore == 0)
        denominator     <- relativeTrend^(components[[j]]$noVertices) * gamma^(components[[j]]$noEdges) * components[[j]]$trendProduct * withinHardcore
        numerator           <- denominator * gamma^(components[[j]]$noNeighbours) * withinHardcore
        interactionProduct  <- interactionProduct * sum(numerator) / sum(denominator)
      }
    }
    return( lambda + betaZ * interactionProduct )
  }
}

#' Computes the log pseudo-likelihood of a given point pattern stemming from a PoSt-hard core process.
#'
#' @param theta the model parameters.
#' @param graphStructure the interaction graph induced by the point pattern as computed by the function compute_graphStructure (1st component).
#' @param integralGridStructure the interaction graph induced by the point pattern and the integration grid as computed by the function compute_graphStructure (2nd component).
#' @param trend the trend/ activation function of the Strauss process. Might be a real number or a pixel image.
#' @param isHomogeneous a flag indicating whether the considered process is homogeneous (TRUE) or not (FALSE)
#'
#' @return The conditional intensity of the i-th point in the point pattern listed in the graph structure passed.
#' @author Johannes Wieditz
#' @export
#'

logPseudoLikelihood  <- function( theta, graphStructure, integralGridStructure, trend, isHomogeneous, windowVolume ){

  PL <- 0
  for( z in 1:length(graphStructure) ){

    CI <- conditionalIntensity(theta, graphStructure = graphStructure, z, trend, isHomogeneous)

    if( CI > 0 ){
      PL    <- PL + log(CI)
    }
  }
  ## Compute the integral over the window of the conditional intensity. Here, we use an equidistant rectangular grid to evaluate the conditional intensity at these points.
  normalisationConstant <- 0
  pointsInWindow <- 0

  for( i in 1:length(integralGridStructure) ){

    if( integralGridStructure[[i]]$isInWindow ){

      CIP <- conditionalIntensity(theta, graphStructure = integralGridStructure, i, trend, isHomogeneous)
      if(!is.na(CIP)){
        normalisationConstant <- normalisationConstant + CIP
        pointsInWindow <- pointsInWindow + 1
      }
    }
  }
  normalisationConstant   <- normalisationConstant / pointsInWindow * windowVolume
  return( - PL + normalisationConstant )
}

#' Computes a maximum pseudo-likelihood estimator for the passed point pattern in the PoSt-hard core model.
#'
#' @param pointPattern the considered point pattern.
#' @param hardCoreRadius the hard core distance $r>0$ of the Strauss process.
#' @param interactionRadius the interaction radius $R > r$ of the Strauss process.
#' @param trend the trend/ activation function of the Strauss process. Might be a real number or a pixel image.
#' @param gridSizeX the number of grid points to compute the normalisation constant of the pseudo-likelihood in x-direction.
#' @param gridSizeY the number of grid points to compute the normalisation constant of the pseudo-likelihood in y-direction.
#' @param theta_init intial value where the optimisation procedure starts from.
#'
#' @return A list of containing the pre-computed graph structures of the point pattern and the integration grid. Note, that the parameters are stored in the natural parameter space, i.e. on log-scale.
#' @author Johannes Wieditz
#' @export
#'
#' @examples
#'
#' \donttest{
#' interactionRadius <- 0.1
#' hardCoreRadius    <- 0.03
#' trend <- 42
#' set.seed(42)
#' straussPoints <- rStraussHard(beta = 1 * trend, gamma = 0.4, R = interactionRadius, H = hardCoreRadius)
#' poissonPoints <- rpoispp(lambda = 12)
#' pointPattern  <- superimpose(straussPoints, poissonPoints)
#' gridSizeX     <- gridSizeY <- 10
#' theta_init    <- c( beta = log(1), gamma = log(.4), lambda = log(12) )
#' MPLE          <- SuperMPLE(pointPattern, interactionRadius, hardCoreRadius, trend, gridSizeX, gridSizeY, theta_init, TRUE)
#' MPLE
#' }
#'
#' \dontrun{
#' data("ExampleFingerprint")
#' gridSizeX <- gridSizeY <- 25
#' theta_init    <- c( beta = log(1), gamma = log(.4), lambda = log(1E-4) )
#' MPLE          <- SuperMPLE(pointPattern, interactionRadius, hardCoreRadius, trend, gridSizeX, gridSizeY, theta_init, FALSE)
#' MPLE
#' }
#'


SuperMPLE <- function( pointPattern, interactionRadius, hardCoreRadius, trend, gridSizeX, gridSizeY, theta_init, isHomogeneous ){

  if( !isHomogeneous ){
    badTrend <- c()
    for( i in 1:pointPattern$n){
      if(is.na(interp.im(trend, pointPattern$x[i], pointPattern$y[i]))){
        badTrend <- c(badTrend, i)
      }
    }
    if( !is.empty(badTrend) )
    pointPattern <- pointPattern[-badTrend]
  }

  L <- compute_graphStructure(pointPattern, interactionRadius, hardCoreRadius, trend, gridSizeX, gridSizeY, isHomogeneous )
  graphStructure        <- L[[1]]
  integralGridStructure <- L[[2]]
  windowVolume          <- volume(Window(pointPattern))

  message('Searching local maximum of the pseudo-likelihood...')
  # Computes the MPLE in the natural parameter space
  MPLE <- optim( par = theta_init,
                 logPseudoLikelihood,
                 graphStructure = graphStructure,
                 integralGridStructure = integralGridStructure,
                 trend = trend,
                 isHomogeneous = isHomogeneous,
                 windowVolume = windowVolume,
                 method = "L-BFGS-B",
                 lower = c(-Inf, -Inf, -Inf),
                 upper = c(Inf, 0, Inf),
                 control = list(trace = 2))
  return(MPLE)
}
