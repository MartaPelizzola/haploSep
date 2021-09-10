#' Plot the reconstruction of haplotypes
#'
#' Function for plotting objects from class haplo
#'
#' @param x an object from class haplo
#'
#' @param type character string giving the type of plot desired
#' \itemize{
#' \item
#' "\code{str}" plot the structure of haplotypes
#' \item
#' "\code{frq}" plot the frequency of haplotypes
#' \item
#' "\code{both}" plot both the structure and the frequency of haplotypes
#' }
#'
#' @param compare an object from class haplo to be compared with x
#'
#' @param ... Arguments to be passed to methods, such as \link{graphical parameters} (see \code{\link[graphics]{par}}).
#'
#'
#' @examples
#'  data(ExampleDataset)
#'  m <- haploSelect(data = Y)
#'  reconstruction <- haploSep(data = Y, nHaplo = m, stabEval = TRUE, bias = TRUE)
#'  
#'  plot(reconstruction, type = "both")
#'
#'
#' @export
#'
#'

plot.haplo <- function(x, type = c("both", "str", "frq"), compare = NULL, ...) {
  type = match.arg(type)
  values <- ind <- dfr <- dfrcomp <- NULL
  # frequency plot
  #     default graphic parameters
  regionCol = "green"
  trpR      = 0.3   # transparency ratio
  txtS      = 15    # text size
  #     build data frames
  nHaplo  = attr(x,"nHaplo")
  nT      = ncol(x$haploFrq)
  tm      = 0:(nT-1)
  dfc     = stack(as.data.frame(t(x$haploFrq)))
  dfc$t   = rep(tm, nHaplo)
  frqPlot = ggplot(dfc, aes(t,values,group=ind,color=ind, linetype = ind)) +
    geom_line() + labs(x="Time",y="Haplotype frequency") + ylim(0,1) +
    scale_x_continuous(breaks = tm, labels=paste0("T", tm)) +
    theme(text=element_text(size=txtS),legend.position="top",
          legend.title=element_blank())
  stabIntFrq = attr(x,"stabIntFrq")
  if (!is.null(stabIntFrq)) {
    dfr     = stack(as.data.frame(t(stabIntFrq[,c(1:nT,nT:1+nT)])))
    dfr$t   = rep(c(tm,rev(tm)), nHaplo)
    dfrF <- dfr
  }
  if (!is.null(compare)){
    if (ncol(compare$haploFrq) != ncol(x$haploFrq))
      stop("Input compare$haploFrq should have the same number of columns as x$haploFrq")
    rownames(compare$haploFrq) <- paste0(rownames(compare$haploFrq), ".comp")
    dfcomp     = stack(as.data.frame(t(compare$haploFrq)))
    dfcomp$t   = rep(tm, nrow(compare$haploFrq))
    frqPlot = frqPlot + geom_line(aes(t,values,group=ind,color=ind), data = dfcomp)
    stabIntFrqC = attr(compare,"stabIntFrq")
    if (!is.null(stabIntFrqC)) {
      rownames(stabIntFrqC) <- paste0(rownames(stabIntFrqC), ".comp")
      dfr_comp     = stack(as.data.frame(t(stabIntFrqC[,c(1:nT,nT:1+nT)])))
      dfr_comp$t   = rep(c(tm,rev(tm)), nrow(compare$haploFrq))
      dfrF <- rbind(dfrF, dfr_comp)
    }
  }
  if (!is.null(stabIntFrq) || !is.null(stabIntFrqC)){
    frqPlot = frqPlot + geom_polygon(aes(t,values,group=ind),
                                     dfrF, linetype=0,
                                     fill=regionCol, alpha=trpR)
  }
  # struture plot
  if (is.null(compare)){
    strPlot <- superheat(x$haploStr, grid.hline.size=0.1,
                         grid.hline.col="gray", pretty.order.rows=FALSE,
                         left.label="none", bottom.label.size=0.1,
                         print.plot=FALSE, legend=FALSE,
                         column.title="Haplotype Structure",
                         row.title="Genomic location")$plot
  } else {
    bm_s <- compare$haploStr[, ind <- haploMatch(x$haploStr, compare$haploStr)]
    od <- order(rowSums(compare$haploFrq[ind,]),decreasing=TRUE)
    bm_s <- bm_s[,od]

    id_hp_s <- c()
    for (j in 1:ncol(x$haploStr)){
      id_hp_s <- cbind(id_hp_s, ifelse(x$haploStr[,j]==bm_s[,j],j,0))
    }
    colnames(id_hp_s) <- paste0("HapID ", 1:ncol(id_hp_s))
    strPlot <- superheat(id_hp_s, pretty.order.cols = F, grid.hline.size=0.1,
                         grid.hline.col="gray", pretty.order.rows=FALSE,
                         left.label="none", bottom.label.size=0.1,
                         print.plot=FALSE, legend=FALSE,
                         column.title="Comparison of Haplotype Structure",
                         row.title="Genomic location")$plot
  }

  if (type == "frq") {
    grid.arrange(frqPlot)
  } else if (type == "str") {
    grid.arrange(strPlot)
  } else if (type == "both") {
    grid.arrange(frqPlot, strPlot, ncol=2)
  }
  return(invisible(list(frqPlot = frqPlot, strPlot = strPlot)))
}
