##
## EOC.R
##
##  Core functions for computing the global FDR for two-sample problems:
##      EOC
##      plot.FDR.result
##      
##
##  Alexander.Ploner@meb.ki.se 231005
##
###############################################################################

"EOC" = function(xdat, grp, p0, paired=FALSE, nperm=25, seed=NULL, plot = TRUE, ...) 
#
# Name: EOC (empirical operating characteristics)
# Desc: displays the FDR and sensitivity for a microarray experiment comparing
#       two groups as a function of the cutoff level
#       Uses reaadl data and requires in most cases package multtest, unless the
#       p-values turn up from somewhere else
# Auth: core by YP, wrapper by AP
# Date: 250405
#
# TODO: add alpha like in TOC?
#       TEST THE PAIRWISE OPTION
#
# Chng: 260405 AP added legend
#       270405 AP added some argument checks
#       170605 AP changed some defaults, output format
#       020905 AP changed pairwise to paired
#       231005 AP split plotting capability off as separate function,
#                 changed output (class etc.)
#       021105 AP changed arguments g to grp, x to xdat 
#       060206 AP corrected class of returned value
#
{
    # Simple check
    if (ncol(xdat) != length(grp)) {
        stop("number of columns in xdat and number of elements in grp do not agree")
    }
    # Process the grouping variable
    ndx = is.na(grp)
    nanum = length(which(ndx))
    if (nanum >0) {
        warning("missing values in grp: ", nanum," case(s) are removed")
        xdat = xdat[, !ndx]
        grp = grp[!ndx]
    }
    # Set the grouping levels to 0/1, as required
    gg = factor(grp) 
    if (length(levels(gg)) != 2) {
        stop("need two distinct groups in g")
    }
    grp = as.numeric(gg) - 1
    
    # Prepare the call of FDRp and check the data
    if (paired) {
        stat = "pairt"
     } else {
        stat = "t.equalvar"
    }
    
    # Do the hard work
    fdrperm = FDRp(xdat,grp,p0=p0, nperm=nperm, seed=seed, test=stat)

    # Process parameters
    p0.est = FALSE
    if (missing(p0)) {
        p0.est = TRUE
        p0 = fdrperm$p0
    }
    
    # Construct the object to be returned
    ret = data.frame(tstat=fdrperm$stat, pvalue=fdrperm$pvalue, FDR=fdrperm$FDR, 
                sens=fdrperm$sens)
    lg = levels(gg)    
    param = list(p0=p0, p0.est=p0.est, statistic = paste("t =",lg[2],"-",lg[1]),
                 paired=paired) 
    attr(ret,"param") = param 
    class(ret) = c("FDR.result","fdr.result","data.frame")

    # do it
    if (plot)
        plot(ret, ...)
    
    ret
}

plot.FDR.result = function(x, add=FALSE, sensitivity.show=TRUE, legend.show=FALSE, 
                           xlim, ylim=c(0,1), xlab, ylab, main, ...)
#
# Name: plot.FDR.result
# Desc: plots the output created by EOC
# Auth: Alexander.Ploner@meb.ki.se 231005
#
# Chng: 2006-02-06 AP added add-switch
#
{
    # Some necessry aux
    t = abs(x$tstat)
    ord = order(t)
    # Process the options
    if (missing(xlab))
        xlab = "Critical value of t statistic"
    if (missing(ylab))
        ylab = "Rates"
    if (missing(xlim)) {
        xlim = range(t)
    }
    p = attr(x, "param")   
    if (missing(main)) {
        main = paste("p0 =",format.pval(p$p0,2))
        main = paste(main, " (",if (p$p0.est) "estimated)" else "fixed)",sep="")
    }
    ## DIRTY: in case of legend, just extend the y-axis
    ## This will work okay provided that ylim=c(0,1) and the device big enough
    ## (As bad as the labels, really)
    if (legend.show) {
        ylim[2] = ylim[1] + (ylim[2]-ylim[1])*1.05
    }    

    # do it
    if (!add) {
        plot(t, x$FDR, type="n", xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
             main=main, ...)
        lines(t[ord], x$FDR[ord], lty=1)
    } else {
        lines(t[ord], x$FDR[ord], lty=1, ...)
    }        
    
    if (sensitivity.show) {
        lines(t[ord], x$sens[ord], lty=3) 
    }

    ## Draw the legend, if required
    if (legend.show) {
        xl = xlim[2]
        yl = ylim[2]
        ll = "FDR"; lty=1; ncol=1
        if (sensitivity.show) {
            ll = c(ll, "Sensitivity"); lty=c(lty,3); ncol=ncol+1
        }
        legend(xl,yl,ll, lty=lty, ncol=ncol, xjust=1, yjust=1)
    }
}


