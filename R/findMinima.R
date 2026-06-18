findLocalMinima <- function(mat, N = 5)
{
  tryCatch({

    shifts <- list(c(1,1), c(0,1), c(1,0), c(-1,-1),
                   c(-1,+1), c(-1,0), c(0,-1), c(1,-1))
    nr <- nrow(mat)
    nc <- ncol(mat)

    applyShifts <- function(mat, shifts, nr, nc) {
      if (shifts[1] == -1) mat <- rbind(mat[-1L, ],    Inf)
      if (shifts[1] ==  1) mat <- rbind(Inf,            mat[-nr, ])
      if (shifts[2] == -1) mat <- cbind(mat[, -1L],    Inf)
      if (shifts[2] ==  1) mat <- cbind(Inf,            mat[, -nc])
      mat
    }

    nmat <- vector("list", 8L)
    for (i in 1:8) {
      nmat[[i]] <- apply(applyShifts(mat, shifts[[i]], nr, nc) > mat, 2, as.numeric)
    }

    isLocalOptima <- Reduce("*", nmat)
    ao <- which(isLocalOptima == 1, arr.ind = TRUE)

    # which() returns a *named vector*, not a matrix, when only 1 cell qualifies.
    # All subsequent [, 1] / [, 2] indexing fails on a vector.
    if (!is.matrix(ao))
      ao <- matrix(ao, nrow = 1L, dimnames = list(NULL, c("row", "col")))

    ord1  <- order(mat[isLocalOptima == 1], decreasing = FALSE)
    ao1   <- ao <- ao[ord1, , drop = FALSE]
    errs1 <- mat[ao]

    aopp <- cbind(as.numeric(rownames(mat)[ao[, 1]]) * 5,
                  as.numeric(colnames(mat)[ao[, 2]]))

    # hclust requires >= 2 points; handle single-optimum case explicitly
    if (nrow(ao) < 2L) {
      clusts1 <- 1L
      bests1  <- 1L
    } else {
      clusts1 <- cutree(hclust(dist(aopp), method = "ward.D2"), h = 0.1)
      bests1  <- unlist(
        tapply(seq_along(clusts1), clusts1,
               function(x) x[which.min(errs1[x])]),
        use.names = FALSE)
    }

    ao   <- ao[bests1, , drop = FALSE]
    ord  <- order(mat[isLocalOptima == 1][ord1][bests1], decreasing = FALSE)
    ao   <- ao[ord,   , drop = FALSE]

    aopp <- cbind(as.numeric(rownames(mat)[ao[, 1]]) * 5,
                  as.numeric(colnames(mat)[ao[, 2]]))

    if (nrow(ao) < 2L) {
      clusts <- 1L
      errs   <- mat[ao]
      bests  <- 1L
    } else {
      clusts <- cutree(hclust(dist(aopp), method = "ward.D2"), h = 0.15)
      errs   <- mat[ao]
      bests  <- unlist(
        tapply(seq_along(clusts), clusts,
               function(x) x[which.min(errs[x])]),
        use.names = FALSE)
      bests  <- bests[order(mat[ao[bests, , drop = FALSE]], decreasing = FALSE)]
    }

    if (length(bests) < N)
      bests <- c(rep(bests[1], N - length(bests)), bests)

    ao <- ao[bests, , drop = FALSE]

    list(bao    = ao[1:N, , drop = FALSE],
         ao     = ao1,
         clusts = clusts1[ord1])

  },
  error   = function(e) { print(e) },
  warning = function(w) { print(w) })
}
