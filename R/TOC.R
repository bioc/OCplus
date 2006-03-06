"TOC" = function(n=10, p0=0.95, sigma=1, D, F0, F1, n1=n, n2=n, paired=FALSE,
                 plot=TRUE, local.show=FALSE, alpha.show=TRUE, 
                 sensitivity.show=TRUE, nplot=100, xlim, ylim=c(0,1), main, 
                 legend.show=FALSE, ...)
#
# Name: TOC
# Desc: displays the FDR for a planned microarray experiment comparing
#       two groups as a function of the cutoff level, either for the 
#       t-statistic or the log fold change
#       By default sharp F0 and simple symmetric F1, but allows both F0 and F1 
#       to be any kind of mixtures of t-distributions
# Auth: core by YP, wrapper by AP
# Date: 291104
#
# TODO: adjust width of label space, legend space: how?
#       - probably by moving the whole plotting thing to lattice and isolating
#         the computational part in doTOC or similar (so TOC can return the 
#         graph as standard for lattice)
#       clean up the code for the x-axis settings
#        
#
# Chng: 161204 AP fixed bug for length(p0)==1
#       190405 AP changed name to TOC (was: Cutoff), cosmetics for plot
#       210405 AP many arguments, a legend of sorts 
#       240405 AP changed argument name to statistic
#       270405 AP removed the fold-change option
#       170605 AP changed return value to data frame
#       020905 AP changed pairwise to paired
#       301005 AP added local.show switch, changed output
#
{
    
    # Process the parameters
    equal.size = (n1==n) & (n2==n)
    if (!equal.size & paired)
        stop("Equal group sizes required for paired test")
    if (missing(F0)) {
        F0 = list(D=0, p=1)
    }
    if (missing(F1)) {
        if (missing(D)) {
            D = 1
        }
        D = unique(abs(D))
        D = c(-D, D)
        F1 = list(D=D, p=rep(1,length(D)))
    } else {
       if (!missing(D)) {
           warning("Both D and F1 specified: D ignored")
        }
    } 

    # Process the distributions
    if (is.null(F0$D))
        F0$D = 0
    if (is.null(F0$p))
        F0$p = rep(1,length(F0$D))
    if (length(F0$D) != length(F0$p))
        stop("F0 incorrectly specified - D and p must have equal length")
    F0$p = F0$p/sum(F0$p)
    if (is.null(F1$D))
        F1$D = c(-1,1)
    if (is.null(F1$p))
        F1$p = rep(1,length(F1$D))
    if (length(F1$D) != length(F1$p))
        stop("F1 incorrectly specified - D and p must have equal length")
    F1$p = F1$p/sum(F1$p)


    # Internally, we use the t-scale - this is the default range
    tlim = c(0,6)
    xlab = "Critical value of t-statistic"
    if (missing(xlim)){         
        xlim = c(0,max(abs(tlim)))
    } else {
        xlim = range(abs(xlim))
        tlim = xlim
    }
    # Set up the plot & internal axis
    xx = seq(xlim[1], xlim[2], length=nplot)    
    tt = seq(tlim[1], tlim[2], length=nplot)
    
    # Process other options
    if (missing(main)) {
        if (equal.size) {
            main = paste("n =",n,"arrays/group")
            if (paired) 
                main = paste(main, "(paired)")
        } else { 
            main = paste("n1 =",n1,"and","n2 =",n2)
        }
    }
    
    # This computes the component distributions for the mixture:
    # - Critical level via F0; note that we do not rely on F0 being symmetric 
    #   (more than slightly paranoid)
    # - Sensitivity level via F1; note that this might very well be asymmetric
    #   (though not by default)
    p0num = length(p0)
    myFDR = mylfdr= matrix(0, nrow=nplot, ncol=p0num)
    if (paired) {    
        alpha = 1 - CDF.paired(tt, n1, F0$D, F0$p, sigma) +
                    CDF.paired(-tt, n1, F0$D, F0$p, sigma)
        sens  = 1 - CDF.paired(tt, n1, F1$D, F1$p, sigma) +
                    CDF.paired(-tt, n1, F1$D, F1$p, sigma)
        for (i in 1:p0num) {
            myFDR[,i]  = FDR.paired(tt, n1, p0[i], F0$D, F0$p, F1$D, F1$p, sigma)
            mylfdr[,i] = lfdr.paired(tt, n1, p0[i], F0$D, F0$p, F1$D, F1$p, sigma)

        }
    } else {
        alpha = 1 - CDF(tt, n1, n2, F0$D, F0$p, sigma) +
                    CDF(-tt, n1, n2, F0$D, F0$p, sigma)
        sens  = 1 - CDF(tt, n1, n2, F1$D, F1$p, sigma) +
                    CDF(-tt, n1, n2, F1$D, F1$p, sigma)
        for (i in 1:p0num) {
            myFDR[,i]  = FDR(tt, n1, n2, p0[i], F0$D, F0$p, F1$D, F1$p, sigma)
            mylfdr[,i] = lfdr(tt, n1, n2, p0[i], F0$D, F0$p, F1$D, F1$p, sigma)
        }
    }        
    
    ## Okay, so this is where we get tricky: put the labels on the curves on 
    ## the side that has more space
    ## Fit in the labels
    lab.txt = format(p0)
    xlim.ext = diff(xlim)*0.05
    sd.left  = ifelse(p0num<2, 1, sd(myFDR[1,]))
    sd.right = ifelse(p0num<2, 2, sd(myFDR[nplot,]))
    if (sd.left > sd.right) {
        lab.ndx = 1
        lab.pos = 2
        xlim2 = c(xlim[1]-xlim.ext, xlim[2])
    } else {
        lab.ndx = nplot
        lab.pos = 4
        xlim2 = c(xlim[1], xlim[2]+xlim.ext)
    }        
    
    ## DIRTY: in case of legend, just extend the y-axis
    ## This will work okay provided that ylim=c(0,1) and the device big enough
    ## (As bad as the labels, really)
    if (legend.show) {
        ylim[2] = ylim[1] + (ylim[2]-ylim[1])*1.05
    }
    
    ## Set up the plot
    if (plot){
        plot(xlim,0:1,type="n", xlab=xlab, ylab="Rates", xlim=xlim2, ylim=ylim, 
         main=main, ...)
    
        ## Draw the curves
        if (alpha.show)
            lines(xx, alpha, lty=2)
        if (sensitivity.show)
            lines(xx, sens, lty=3)
        ff = if (local.show) mylfdr else myFDR
        for (i in 1:p0num) {
            lines(xx,ff[, i])
            text(xx[lab.ndx],ff[lab.ndx,i],lab.txt[i], cex=0.8, pos=lab.pos, offset=0.1)        
        }
    
        ## Draw the legend, if required
        if (legend.show) {
            xl = xlim2[2]
            yl = ylim[2]
            ll = if (local.show) "fdr" else "FDR"; lty=1; ncol=1
            if (alpha.show) {
                ll = c(ll, "Alpha"); lty=c(lty,2); ncol=ncol+1
            }
            if (sensitivity.show) {
                ll = c(ll, "Sensitivity"); lty=c(lty,3); ncol=ncol+1
            }
            legend(xl,yl,ll, lty=lty, ncol=ncol, xjust=1, yjust=1)
        }
    }
        
    ## A nice object to return invisibly
    ret1 = myFDR
    colnames(ret1) = paste("FDR",lab.txt,sep="")
    ret2 = mylfdr
    colnames(ret2) = paste("fdr",lab.txt,sep="")
    ret = cbind(stat=tt, ret1, ret2)
    if (alpha.show) {
        ret = cbind(ret, alpha=alpha)
    }
    if (sensitivity.show) {
        ret = cbind(ret, sensitivity=sens)
    }
    ret = as.data.frame(ret)
    param = list(p0=p0, n1=n1, n2=n2, sigma=1, F0=F0, F1=F1, paired=paired) 
    attr(ret,"param") = param
    invisible(ret)
}

