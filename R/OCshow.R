OCshow = function(x, ..., global=TRUE, percentage=TRUE, top=0.1, legend, 
                  lty, col, main, xlab, ylab)
#
# Name: OCshow
# Desc: plot the operating characteristics for the output of EOC, fdr1d, fdr2d,
#       i.e. plot the FDR (or fdr) as a function of the percentage of genes that
#       are declared DE based on their gene-wise fdr or FDR
# Auth: Alexander.Ploner@meb.ki.se 251005
#
# Chng: 021105 AP changed x-label to proportion
#
{
    # Collect the arguments, see boxplot.default
    args = list(x, ...)
    # Check the arguments for class
    cl = sapply(args, function(x) class(x)[1])
    if (length(setdiff(cl, c("fdr1d.result","fdr2d.result","FDR.result")))>0)
        stop("arguments must have class fdr.result or FDR.result")
    # Local fdr only available for all fdr cast
    if (!global & "FDR.result" %in% cl) {
        warning("argument of class FDR.result - global=FALSE will be ignored")
        global = TRUE
    }
    # Some graphic defaults
    k = length(args)
    if (missing(lty)) lty = 1:k
    if (missing(col)) col = rep(par("fg"), k)
    if (missing(xlab)) xlab = if (percentage) "Proportion DE" else "Number DE"
    if (missing(ylab)) ylab = if (global) "FDR" else "fdr"
    # The worker 
    xtract = function(x)
    {
        n = nrow(x)
        is.FDR = class(x)[1]=="FDR.result"
        Fdr = if (is.FDR)  x$FDR else x$fdr
        ord = order(Fdr)
        Fdr = Fdr[ord]
        cnt = 1:n
        prc = cnt/n
        ndx = prc <= top
        Fdr = Fdr[ndx]
        cnt = cnt[ndx]
        prc = prc[ndx]
        if (global & !is.FDR)
            Fdr = cumsum(Fdr)/(1:length(Fdr))
        list(count=cnt, perc=prc, Fdr=Fdr)
    }
    pp = lapply(args, xtract)
    # Prepare the range
    xr = range(unlist(sapply(pp, function(x) range(if (percentage) x$perc else x$count))))
    yr = range(unlist(sapply(pp, function(x) range(x$Fdr))))
    plot(xr, yr, type="n", xlab=xlab, ylab=ylab)
    for (i in 1:k) {
        ll = pp[[i]]
        lines(if (percentage) ll$perc else ll$count, ll$Fdr, lty=lty[i], col=col[i])
    }
    if (!missing(legend))
        legend("topleft", legend=legend, lty=lty, col=col)
    if (!missing(main))
        title(main)
}

