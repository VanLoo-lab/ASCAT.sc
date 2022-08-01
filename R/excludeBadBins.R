excludeBadBins <- function(lSe, chr, dna, max.N.freq=.05)
{
	if(is.null(dna)) return(lSe)
    require(Biostrings)
    which.bad <- which(sapply(1:length(lSe$starts),function(x)
    {
        prop.table(alphabetFrequency(dna[[chr]][lSe$starts[x]:lSe$ends[x]]))["N"]>max.N.freq
    }))
    nlSe <- list(starts=lSe$starts[-c(which.bad)],
                 ends=lSe$ends[-c(which.bad)])
    nlSe
}

# fix from: Haixi Yan
excludeBadBins <- function(lSe, chr, dna, max.N.freq = 0.05) {
      if (is.null(dna))
        return(lSe)
      require(Biostrings)
      which.bad <- which(sapply(1:length(lSe$starts), function(x) {
        prop.table(alphabetFrequency(dna[[chr]][lSe$starts[x]:lSe$ends[x]]))["N"] >
          max.N.freq
      }))
      if(length(which.bad) != 0) {
        nlSe <- list(starts = lSe$starts[-c(which.bad)], ends = lSe$ends[-c(which.bad)])
      } else {
        nlSe <- list(starts = lSe$starts, ends = lSe$ends)
      }
      return(nlSe)
    }
