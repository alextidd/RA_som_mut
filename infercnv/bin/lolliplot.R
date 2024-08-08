# jitterLabels and reAdjustLabels adapted from lolliplots package
jitterLabels <- function(coor, xscale, lineW, weight = 1.2) {
  if (weight == 1.2) {
    stopifnot(
      "Please sort your inputs by start position" =
        order(coor) == 1:length(coor)
    )
  }
  if (weight < 0.5) {
    return(coor)
  }
  stopifnot(length(xscale) == 2)
  pos <- (coor - 1) / (xscale[2] - 1)
  # pos <- convertX(unit(coor01, "native"), "npc", valueOnly=TRUE)
  pos.diff <- diff(c(0, pos, 1))
  idx <- which(pos.diff < weight * lineW)
  if (length(idx) < 1) {
    return(coor)
  }
  if (all(idx %in% c(1, length(pos) + 1))) {
    return(coor)
  }
  idx.diff <- diff(c(-1, idx))
  idx.grp <- rle(idx.diff)
  idx.grp$values[idx.grp$values == 1] <- length(pos) + 1:sum(idx.grp$values == 1)
  idx.grp <- inverse.rle(idx.grp)
  idx.grp.w <- which(idx.grp > length(pos)) - 1
  idx.grp.w <- idx.grp.w[idx.grp.w > 0]
  idx.grp[idx.grp.w] <- idx.grp[idx.grp.w + 1]
  idx.grp <- split(idx, idx.grp)
  flag <- as.numeric(names(idx.grp)) > length(pos)
  idx.grp.mul <- lapply(idx.grp[flag], function(.ele) {
    c(.ele[1] - 1, .ele)
  })
  idx.grp.sin <- lapply(idx.grp[!flag], function(.ele) {
    lapply(as.list(.ele), function(.ele) {
      c(.ele - 1, .ele)
    })
  })
  idx.grp.sin <- unlist(idx.grp.sin, recursive = FALSE)
  idx.grp <- c(idx.grp.mul, idx.grp.sin)

  adj.pos <- lapply(idx.grp, function(.ele) {
    .ele <- .ele[.ele > 0 & .ele <= length(pos)]
    this.pos <- pos[.ele]
    names(this.pos) <- .ele
    if (length(this.pos) %% 2 == 1) {
      center <- ceiling(length(this.pos) / 2)
    } else {
      center <- length(this.pos) / 2 + .5
    }
    if (length(this.pos) > 5) { ## too much, how to jitter?
      this.pos <- this.pos +
        ((1:length(this.pos)) - center) * (weight - .1) *
          lineW / ceiling(log(length(this.pos), 5))
    } else {
      this.pos <- this.pos +
        ((1:length(this.pos)) - center) * (weight - .1) * lineW
    }
    this.pos
  })
  names(adj.pos) <- NULL
  adj.pos <- unlist(adj.pos)
  coor[as.numeric(names(adj.pos))] <- adj.pos * diff(xscale) + xscale[1]

  Recall(coor, xscale = xscale, lineW = lineW, weight = weight - 0.2)
}

reAdjustLabels <- function(coor, lineW, xmax) {
  # resort
  coor <- sort(coor)
  bins <- ceiling(1 / lineW)
  pos <- (coor - 1) / (xmax - 1)
  # pos <- convertX(unit(coor, "native"), "npc", valueOnly=TRUE)
  pos.bin <- cut(pos, c(-Inf, (0:bins) * lineW, Inf), labels = 0:(bins + 1), right = FALSE)

  ## split the coors by into clusters
  ## give the clusters with more idx more spaces if there are spaces between clusters
  tbl <- table(pos.bin)
  if (all(tbl < 2)) {
    return(coor)
  }
  tbl.len <- length(tbl)
  if (tbl.len < 3) {
    return(coor)
  }
  loops <- 1000
  loop <- 1
  while (any(tbl == 0) && any(tbl > 1) && loop < loops) {
    tbl.bk <- tbl
    for (i in order(tbl.bk, decreasing = TRUE)) {
      if (tbl[i] > 1 && tbl.bk[i] == tbl[i]) {
        if (i == 1) {
          if (tbl[2] < tbl[1]) {
            half <- sum(tbl[1:2]) / 2
            tbl[2] <- ceiling(half)
            tbl[1] <- floor(half)
          }
        } else {
          if (i == tbl.len) {
            if (tbl[tbl.len] > tbl[tbl.len - 1]) {
              half <- sum(tbl[(tbl.len - 1):tbl.len]) / 2
              tbl[tbl.len - 1] <- ceiling(half)
              tbl[tbl.len] <- floor(half)
            }
          } else {
            if (tbl[i - 1] < tbl[i + 1]) {
              ## i-1 and i should be balanced
              half <- sum(tbl[(i - 1):i]) / 2
              tbl[i - 1] <- floor(half)
              tbl[i] <- ceiling(half)
            } else {
              half <- sum(tbl[i:(i + 1)]) / 2
              tbl[i] <- floor(half)
              tbl[i + 1] <- ceiling(half)
            }
          }
        }
      }
    }
    loop <- loop + 1
  }
  coef <- unlist(lapply(tbl, function(.ele) {
    if (.ele == 0) {
      return(0)
    }
    .ele <- seq(from = 0, to = 1, length.out = .ele + 1)
    (.ele[-length(.ele)] + .ele[-1]) / 2
  }))
  coef <- coef[coef != 0]
  coor <- (rep(as.numeric(names(tbl)), tbl) - 1 + coef) * lineW
  # coor <- convertX(unit(coor, "npc"), "native", valueOnly=TRUE)
  coor <- ((xmax - 1) * coor) + 1
  coor
}

# Covert genomic bp to coding sequence coordinate
getCDSindex <- function(intervals) {
  strand <- intervals$strand[1]
  allBP <- tibble::tibble()
  for (exon in 1:length(intervals$length)) {
    exonBP <- tibble::tibble(pos = seq(intervals$start[exon], intervals$end[exon]))
    allBP <- dplyr::bind_rows(allBP, exonBP)
  }
  # Forward strand standard
  if (strand == 1) {
    allBP <- allBP %>%
      dplyr::mutate(cdsPos = seq(1:sum(intervals$length)))
  } else {
    # Reverse strand flip cds coords so all genes show left to right
    allBP <- allBP %>%
      dplyr::mutate(cdsPos = seq(from = sum(intervals$length), to = 1))
  }
  return(allBP)
}

# For a mutation outside the CDS (splice or non-coding) set its CDS position as .5 off the nearest exon border
getCDSindexNonCDS <- function(cdsCoordsIndex, vars, strand) {
  nonCDSvars <- vars %>% dplyr::filter(impact %in% c("Essential_Splice", "Non-coding") | (impact == "Indel" & is.na(ntchange)))
  closestCDS <- cdsCoordsIndex %>%
    dplyr::slice(rep(1:dplyr::n(), each = dim(nonCDSvars)[1])) %>%
    dplyr::mutate(nonCDSpos = rep(nonCDSvars$POS, times = dim(cdsCoordsIndex)[1])) %>%
    dplyr::mutate(gap = nonCDSpos - pos) %>%
    dplyr::mutate(absGap = abs(gap)) %>%
    dplyr::group_by(nonCDSpos) %>%
    dplyr::filter(absGap == min(absGap)) %>%
    dplyr::ungroup() %>%
    # flip direction if reverse strand
    dplyr::mutate(gap = gap * strand) %>%
    # Set offset based on strand direction
    dplyr::mutate(cdsPos = dplyr::if_else(gap < 0, cdsPos - 0.5, cdsPos + 0.5)) %>%
    dplyr::select(POS = nonCDSpos, cdsPos)
  return(closestCDS)
}
