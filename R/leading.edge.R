leading.edge<-function(k, #TRUE if genes in genelist
                       forder,#order(s2n)
                       s2n,#scores
                       meanES,
                       Ngenes,#number of genes
                       doplot=F,
                       label="",
                       alless=NA,
                       col="black",
                       add=F,
                       pas=10,
                       ...)
{
    running.sumt <- vector(length = length(s2n))
    emplacement <- k[forder]
    rest <- !emplacement
    sum.k.signal <- sum(abs(s2n[k]))
    running.sumt <- s2n[forder]
    running.sumt[rest] <- -abs(1/(Ngenes - sum(k)))
    running.sumt[emplacement] <- abs(running.sumt[emplacement])/sum.k.signal
    running.sum <- cumsum(running.sumt)
    ylims <- c(-1, 1)
    sign <- if (max(running.sum) < abs(min(running.sum)))
                -1
            else 1
    res <- list(ess = max(abs(running.sum)), sign = sign, essmin = min(running.sum),
                essmax = max(running.sum))
    if (doplot) {
        if (add) {
            points(1:Ngenes, running.sum, type = "l", pch = ".",
                   cex = 1, col = col)
        }
        else {
            if (!is.na(alless))
                plot(1:Ngenes, running.sum, type = "l", pch = ".",
                     cex = 1, frame.plot = F,  xaxt = "n",
                     ylim = c(min(alless), max(alless)) * 1.5, col = col,
                     ...)
            else plot(1:Ngenes, running.sum, type = "l", lwd = 6,
                      frame.plot = F, cex.axis = 2,
                      cex.lab = 2, ylim = ylims, xaxt = "n", col = col,
                      ...)
            mtext(side=3, label, cex=.5)
            abline(h = 0, col = rgb(0.1, 0.1, 0.1), lty = 2,
                   lwd = 5)
        }
        ##mx <- quantile(s2n, probs = 0.9, na.rm=T)
        ##mn <- quantile(s2n, probs = 0.1, na.rm=T)
        s2nf <- s2n[forder]
        mtext(paste("ES = ", signif(res$ess * res$sign, 2), "  ",
                    "NES = ", signif(res$ess * res$sign/abs(meanES),
                                     2), sep = ""), side = 3, cex = 2)
    }
    return(list(res = res, rug = which(emplacement)))
}
