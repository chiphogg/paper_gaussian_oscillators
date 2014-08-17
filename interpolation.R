require(ggplot2)
require(reshape2)

TimeSequence <- function(n.keyframes, n.points, loop=FALSE) {
  # A sequence from 0 to n.keyframes, with n.points outputs.
  #
  # Args:
  #   n.keyframes:  The number of keyframes -- i.e., independent Gaussian draws
  #     -- to use.
  #   n.points:  The number of points to show in the output.  Increasing this
  #     number gives smoother functions but bigger output files and longer
  #     computation times.
  #   loop:  Indicates this is intended for a looping animation.  If FALSE, the
  #     final time will be equal to n.keyframes.  Otherwise, the final-plus-one
  #     time would be equal to n.keyframes, but this isn't necessary because
  #     it's equivalent to t=0.
  dt <- n.keyframes / (n.points - ifelse(loop, 0, 1))
  seq(from=0, by=dt, length.out=n.points)
}

InterpolationConceptPlot <- function(
    interp.matrix.functions,
    n.points,
    n.keyframes) {
  # A figure to demonstrate the concept of interpolation: generate random
  # points; interpolate between them.
  #
  # Args:
  #   interp.matrix.functions:  A list of functions mapping
  #     (<number of keyframes>, <output times>) -> <interpolation matrix>.
  #   n.points:  The number of points to show in the output.  Increasing this
  #     number gives smoother functions but bigger PDF files and longer
  #     computation times.
  #   n.keyframes:  The number of keyframes -- i.e., independent Gaussian draws
  #     -- to use.
  #
  # Returns:
  #   A ggplot2 plot object showing random points, and different methods of
  #   interpolating between them.
  seeds <- rnorm(n=n.keyframes)
  seed.data <- data.frame(t=1:n.keyframes, value=seeds)
  t <- TimeSequence(n.keyframes, n.points)
  d.f <- cbind(t, data.frame(sapply(interp.matrix.functions,
                             function(FUN) FUN(n.keyframes, t) %*% seeds)))
  dfm <- melt(d.f, id.vars='t')
  (ggplot(data=dfm, aes(x=t, y=value))
    + geom_line(aes(colour=variable))
    + geom_point(data=seed.data)
    )
}

RowQuantileDataFrame <- function(mat, t, quantiles) {
  # Compute the quantiles for each row of the matrix.  Rows are assumed to
  # correspond to times, according to the vector t.
  #
  # Args:
  #   mat:  A matrix whose rows correspond to times, and whose columns
  #     correspond to independent Gaussian oscillators.
  #   t:  A vector giving the time corresponding to each row.
  #   quantiles:  Which quantiles to compute over the oscillators at each time.
  #
  # Returns:
  #   A data.frame giving the time, the index of the quantile, and the value of
  #   that quantile at that time.
  n.quantiles <- length(quantiles)
  computed.quantiles <- apply(mat, 1, quantile, probs=quantiles, type=4)
  data.frame(t=rep(t, each=n.quantiles),
             y=as.vector(computed.quantiles),
             group=rep(1:n.quantiles, nrow(mat)))
}

InterpolationQuantilePlot <- function(
    interp.matrix.functions, quantiles, n.draws, n.points, n.keyframes) {
  # Show selected quantiles for various methods of interpolation.
  #
  # Args:
  #   interp.matrix.functions:  A list of functions mapping
  #     (<number of keyframes>, <output times>) -> <interpolation matrix>.
  #   n.draws:  The number of functions to draw from the distribution for
  #     calculating quantiles.  Increasing this number gives better statistics
  #     but longer computation times.
  #   n.points:  The number of points to show in the output.  Increasing this
  #     number gives smoother functions but bigger PDF files and longer
  #     computation times.
  #   n.keyframes:  The number of keyframes -- i.e., independent Gaussian draws
  #     -- to use.
  #
  # Returns:
  #   A ggplot2 plot object showing the quantiles as a function of time.  (This
  #   visually conveys which interpolation methods are statistically correct.)
  t.out <- TimeSequence(n.keyframes, n.points, loop=FALSE)
  n.seeds <- n.keyframes * n.draws
  random.seeds <- matrix(rnorm(n=n.seeds), nrow=n.keyframes)
  d.frame <- data.frame()
  for (interp.method in names(interp.matrix.functions)) {
    f <- interp.matrix.functions[[interp.method]]
    new.frame <- RowQuantileDataFrame(f(n.keyframes, t.out) %*% random.seeds,
                                      t.out,
                                      quantiles)
    new.frame$label <- interp.method
    d.frame <- rbind(d.frame, new.frame)
  }
  (ggplot(data=d.frame, aes(x=t, y=y))
    + geom_line(aes(group=paste0(group, label), colour=label))
    )
}

# Functions for constructing interpolation matrices.
weight <- function(t.in, t.out, FUN) {
  ifelse(abs(t.in - t.out) > 1, 0, FUN(t.in, t.out))
}
periodic_weight <- function(t.in, t.out, period, FUN) {
  # Dirty hack to make the basis functions periodic (instead of always 0 at
  # t=0).
  pmax(weight(t.in, t.out, FUN), weight(t.in, t.out - period, FUN))
}
InterpMatrix <- function(n.frames, t.out, FUN) {
  t.frames <- 1:n.frames
  outer(t.out, t.frames, FUN)
}
InterpLinear <- function(n.frames, t.out) {
  InterpMatrix(n.frames, t.out, 
               function(x, y) {
                 periodic_weight(x, y, n.frames,
                                 function(a, b) 1 - abs(a - b))
               })
}
InterpTrig <- function(n.frames, t.out) {
  InterpMatrix(n.frames, t.out, 
               function(x, y) {
                 periodic_weight(x, y, n.frames,
                                 function(a, b) cos((a - b) * pi / 2))
               })
}
InterpSpline <- function(n.frames, t.out) {
  SplineBasis <- function(i) {
    y <- rep(0, n.frames)
    y[i] <- 1
    spline(method='periodic', xout=t.out, x=(0:n.frames) + 1, y=c(y, y[1]))$y
  }
  sapply(1:n.frames, SplineBasis)
}
DelocalizedMatrix <- function(n.frames, t.out) {
  # n.frames must be even (we need as many sines as cosines).
  stopifnot(n.frames %% 2 == 0)
  N <- n.frames / 2
  cbind(
    outer(t.out, 1:N, function(t, i) sin(pi * i * t / N)),
    outer(t.out, N:1, function(t, i) cos(pi * i * t / N))) / sqrt(N)
}
