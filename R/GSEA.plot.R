GSEA.plot<-function(m,
                    geneset,
                    o=NULL,
                    nperms=1000,
                    metric=c("SNR","FC"),
                    col=rgb(.3,.73,.3))
  {
      scores <- get(paste("snr.", metric, sep = ""))(m, o)
      scoresorder <- order(scores, decreasing = T)
      k <- rownames(m) %in% geneset
      names(scores) <- rownames(m)
      if(is.na(meanES))
      {
          scrandom <- lapply(1:nperms, function(x) {
              s2n <- if (is.null(o))
                         get(paste("snr.", metric, sep = ""))(apply(m, 2,
                                                                    function(x) sign(rnorm(1)) * x), o)
              else get(paste("snr.", metric, sep = ""))(m, sample(o))
              forder <- order(s2n, decreasing = T)
              leading.edge(k = k, forder = forder, s2n = s2n, Ngenes = nrow(m),
                           doplot = FALSE, alless = NA, col = "black", add = F,
                           pas = 10)
          })
          meanES <- mean(unlist(lapply(scrandom, function(x) abs(c(x$res$essmax,x$res$essmin)))))
      }
      layout(mat = cbind(c(1, 2, 3)),
             widths = 1,
             heights = c(6,
                         1, 4))
      par(mar = c(0, 6, 5, 1))
      sc <- leading.edge(k = k, forder = scoresorder,
                         meanES = meanES,
                         s2n = scores, Ngenes = nrow(m),
                         doplot = TRUE, alless = NA,
                         col = col, add = F, pas = 10, xlab = "",
                         ylab = "Enrichment score (ES)", ...)
      par(mar = c(0, 6, 0, 1))
      plot(0, 0, col = rgb(0, 0, 0, 0), ylim = c(-1, 1),
           xlim = c(1,
                    nrow(m)),
           frame.plot = F, xaxt = "n", yaxt = "n", ylab = "",
           xlab = "")
      abline(v = sc$rug, col = rgb(0.2, 0.2, 0.2, 0.8))
      par(mar = c(4, 6, 0, 1))
      plot(1:nrow(m), scores[order(scores, decreasing = T)], type = "l",
           lwd = 2, xlab = "Rank",
           ylab = if (metric[1] == "SNR")
                                              "Signal to Noise"
                  else "Fold Change",
           cex.axis = 2, cex.lab = 2, main = "",
           frame.plot = F)
      return(meanES)
  }
