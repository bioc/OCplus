"samplesize" = function(n=seq(5,50, by=5), p0=0.99, sigma=1, D, F0, F1, 
                        paired=FALSE, crit, crit.style=c("top percentage", "cutoff"), 
                        plot=TRUE, local.show=FALSE, nplot=100, ylim=c(0,1), main, 
                        legend.show=FALSE, grid.show=FALSE, ...)
#
# Name: samplesize
# Desc: displays the FDR for a planned microarray experiment comparing
#       two groups as a function of sample size for different cutoffs on the
#       t-statistic AND returns a matrix of likely values
#       By default sharp F0 and simple symmetric F1, but allows both F0 and F1 
#       to be any kind of mixtures of t-distributions
# Auth: core by YP, wrapper by AP
# Date: 291104
#
# TODO: Think about it - how would unequal sample sizes look here?
#
# Chng: 190405 AP changed name to samplesize (was SampleSize)
#       210405 AP plot cosmetics, arguments
#       260405 AP added legend
#       270405 AP added parameter checks, removed fold-change option
#       170605 AP added pairwise option
#       020905 AP changed pairwise to paired
#       231005 AP corrected default for paired to FALSE
#       301005 AP added local.show, changed output
#
{
    # Process the parameters
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
    
    
    # Process the critical default values; for now, default values depend 
    # solely on the kind of selection, not the
    crit.style = match.arg(crit.style)    
    if (missing(crit)) {
        crit = switch(crit.style,"cutoff"=2, "top percentage"=0.01)
    }
    # Check the critical values
    if (crit.style=="cutoff") {
        if (any(crit <= 0)) 
            stop("Only positive critical values allowed for crit.style==cutoff")
    } else if (crit.style=="top percentage") {
        if (any((crit>=1) | (crit<=0))) 
            stop("Critical values must be between 0 and 1 for crit.style==top percentage")
    }
    # Check the p0
    if (length(p0)>1) {
        warning("Only one p0 per run possible - additional values ignored")
        p0 = p0[1]
    }
    
    # Process other options
    if (missing(main)) {
        if (paired)
            main = substitute(paste(p[0]==p0," (paired)"), list(p0=p0))
        else 
            main =  substitute(p[0]==p0, list(p0=p0))
    }
    ## DIRTY: in case of legend, just extend the y-axis
    ## This will work okay provided that ylim=c(0,1) and the device big enough
    ## (As bad as the labels, really)
    if (legend.show) {
        ylim[2] = ylim[1] + (ylim[2]-ylim[1])*1.05
    }    
    
    # Plotting more than tabulating: we have a different set of points for
    # return values and plotting
    nn = seq(min(n), max(n), length=nplot)    
    retnum  = length(n)
    critnum = length(crit)
 
    # The way we treat the fold changes is not really kosher; right now, we
    # just find the critical values from the t-statistic and blow them up 
    # afterwards
    # If we have cutoff, we can set up the matrix of critical values directly
    if (crit.style=="cutoff") {
        crit.ret   = matrix(crit, nrow=retnum, ncol=critnum, byrow=TRUE)
        crit.plot  = matrix(crit, nrow=nplot, ncol=critnum, byrow=TRUE) 
    # In case of the top percentage, we actually have to invert the distribution
    # function; note that we use global variables in qFmix
    } else if (crit.style=="top percentage") {
        if (paired) {
            qFmix = function(x) {
                lq = CDFmix.paired(-abs(x), yy, p0, F0$D, F0$p, F1$D, F1$p, sigma) 
                uq = 1 - CDFmix.paired(abs(x), yy, p0, F0$D, F0$p, F1$D, F1$p, sigma) 
                lq + uq - cc
            }
        } else { 
            qFmix = function(x) {
                lq = CDFmix(-abs(x), yy, yy, p0, F0$D, F0$p, F1$D, F1$p, sigma) 
                uq = 1 - CDFmix(abs(x), yy, yy, p0, F0$D, F0$p, F1$D, F1$p, sigma) 
                lq + uq - cc
            }
        }    
        # These are the default limits for the interval where to search for the root
        # These settings are entirely heuristical in that they work as-is
        # for the current set of examples(samplesize)
        # We rely on the "extendInt"-argument to uniroot() to adjust for
        # situations where the limit needs to be extended (upwards)
        from = 0
        to   = 15
        # Do it
        crit.ret  = matrix(0, nrow=retnum, ncol=critnum)
        crit.plot = matrix(0, nrow=nplot, ncol=critnum)
        for (i in 1:critnum) {
            cc = crit[i]
            for (j in 1:retnum) {
                yy = n[j]
				if (qFmix(from) <= 0) {
					stop("Weird distribution, can't run uniroot from 0")
				}
                crit.ret[j,i] = uniroot(qFmix, c(from, to), extendInt="downX")$root
            }
            for (j in 1:nplot) {
                yy = nn[j]
				if (qFmix(from) <= 0) {
					stop("Weird distribution, can't run uniroot from 0")
				}
                crit.plot[j,i] = uniroot(qFmix, c(from, to), extendInt="downX")$root
            }
        }
    }
        
    # Compute the return values    
    ret = retl = matrix(0, nrow=retnum, ncol=critnum)
    for (i in 1:critnum) {
        ret[, i] = if (paired) 
                        FDR.paired(crit.ret[, i], n, p0, F0$D, F0$p, F1$D, F1$p, sigma)
                   else
                        FDR(crit.ret[, i], n, n, p0, F0$D, F0$p, F1$D, F1$p, sigma)
        retl[, i] = if (paired) 
                        lfdr.paired(crit.ret[, i], n, p0, F0$D, F0$p, F1$D, F1$p, sigma)
                   else
                        lfdr(crit.ret[, i], n, n, p0, F0$D, F0$p, F1$D, F1$p, sigma)

    }
    
    ## Okay, so this is where we get tricky: put the labels on the curves on 
    ## the side that has more space; we can use the return values legitimately
    ## for this!
    xlim = range(n)
    xlim.ext = diff(xlim)*0.05
    sd.left  = ifelse(critnum>1, sd(ret[1,]), 0)
    sd.right = ifelse(critnum>1, sd(ret[retnum,]), 1)
    if (sd.left > sd.right) {
        lab.ndx = 1
        lab.pos = 2
        xlim2 = c(xlim[1]-xlim.ext, xlim[2])
    } else {
        lab.ndx = nplot
        lab.pos = 4
        xlim2 = c(xlim[1], xlim[2]+xlim.ext)
    }        
    
    ## Set up the plot
    if (plot) {
        plot(xlim, 0:1,type="n", xlab="Sample size (arrays/group)", 
             ylab=if (local.show) "fdr" else "FDR", xlim=xlim2, ylim=ylim, 
             main=main, ...)
    
        ## Draw the curves
        for (i in 1:critnum) {
            ff = if (paired & local.show) 
                    lfdr.paired(crit.plot[, i], nn, p0, F0$D, F0$p, F1$D, F1$p, sigma)
                 else if (!paired & local.show)
                    FDR(crit.plot[, i], nn, nn, p0, F0$D, F0$p, F1$D, F1$p, sigma)                    
                 else if (paired & !local.show) 
                    FDR.paired(crit.plot[, i], nn, p0, F0$D, F0$p, F1$D, F1$p, sigma)
                 else if (!paired & !local.show)
                    FDR(crit.plot[, i], nn, nn, p0, F0$D, F0$p, F1$D, F1$p, sigma)                    
            lines(nn,ff)
            text(nn[lab.ndx],ff[lab.ndx],crit[i], pos=lab.pos, offset=0.2)        
        }
    
        ## Draw the legend, if required
        if (legend.show) {
            xl = xlim2[2]
            yl = ylim[2]
            legend(xl,yl,paste(crit.style), xjust=1, yjust=1)
        }    
    
        # Fill in a grid for the desired values, if required
        if (grid.show) {
            minn = min(n)
            minf = min(ylim)
            for (i in 1:retnum) {
                for (j in 1:critnum) {
                    points(n[i], ret[i,j])
                    segments(minn, ret[i,j], n[i], ret[i,j], lty=3)
                    segments(n[i], minf, n[i], ret[i,j], lty=3)
                }
            }
        }
    }

    rownames(ret) = n
    colnames(ret) = paste("FDR",crit,sep="_")
    colnames(retl) = paste("fdr",crit,sep="_")
    ret = as.data.frame(cbind(ret, retl))    
    param = list(p0=p0, sigma=sigma, F0=F0, F1=F1, paired=paired, 
                 crit.style=crit.style)
    attr(ret,"param") = param    
    ret
}

