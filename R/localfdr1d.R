##
## localfdr1d.R
##
##  Routines for computing the one-dimensional local false discovery rate
##      fdr1d
##      plot.fdr1d.result
##      smooth1d
##      Rinv
##      
##
##  This is based on earlier code from package eBayes
##
##  Alexander.Ploner@meb.ki.se  120705
##
##      090805 prepared for OCplus
##      160805 changed 1d smoothing
##      310805 added to OCplus
## 
###############################################################################

##------------------ The main access function -----------------------------##

fdr1d = function (xdat, grp, test, p0, nperm=100, nr=50, seed=NULL, null=NULL, 
                  zlim=1, sv2=0.01, err=0.0001, verb=TRUE, ...)
#
# Name: fdr1d
# Desc: 
# Uses: 
# Auth: Alexander.Ploner@meb.ki.se 120705
#       based on earlier code by AP & Yudi.Pawitan@meb.ki.se 
#
# TODO:
#
# Chng: 090805 added null, p0 as parameters
#       160805 added automatic B-spline/Poisson smoothing
#       031105 changed function name to tstatistics
#
{
    vcat = function(...) if (verb) cat(...)
    
    if (missing(test))
        test = function(...) tstatistics(...,logse=FALSE)[,1]

    # some validity checks
    if (ncol(xdat)!=length(grp))
        stop("grouping variable has invalid length")
    # If null is explicitly specified, we don't need nperm or seed
    if (!is.null(null)) {
        if (!is.null(seed)) 
            warning("Parameter null overrides parameter seed - ignored\n")
    }
    # remove NAs from the expression data, grouping
    if (any(is.na(grp))) {
        Prefilt = !is.na(grp)
        xdat = xdat[, Prefilt]
        grp  = grp[Prefilt]
        vcat("Removed",length(which(!Prefilt)), "samples with missing grp\n")
    }
    ng = nrow(xdat)
    nc = ncol(xdat)

    # compute the original statistics
    Z = test(xdat, grp, ...)
    
    # set the bin range pragmatically - we don't want extreme values
    xbreaks = MidBreaks(Z, nr)
    xmids   = brk2mid(xbreaks)
    
    # Do the permutations: we are not efficient here in that we store them 
    # all instead of binning them immediately
    if (is.null(null)) {
        zstar = rep(0, ng*nperm)
        if (!is.null(seed)) set.seed(seed)
            vcat("Starting permutations...\n")        
        for (i in 1:nperm){
            vcat("\t", i,"\n")
            perm = sample(nc)
            zstar[((i-1)*ng+1):(i*ng)] = test(xdat, grp[perm])
        }
    } else { # Use pre-packaged null - not efficient either (double allocation)
        nperm = length(null)/ng
        if (nperm != floor(nperm)) 
            stop("Permutations not multiple of number of genes")
        zstar = null
        vcat("Using pre-computed permutations (nperm =",nperm,", seed =",
             attr(null, "seed"),")\n")
    }

    # Table it & smooth it
    vcat("Smoothing f\n")    
    count = hist(Z, xbreaks, plot=FALSE)$counts
    scount = smooth1d(count, err=err, sv2=sv2, verb=FALSE)$fit
    fz = scount/(length(Z)*(xbreaks[2]-xbreaks[1]))  
    vcat("Smoothing f0\n")    # Is this really necessary??
    zstar = zstar[xbreaks[1] <= zstar & zstar <= xbreaks[nr]]
    count = hist(zstar, xbreaks, plot=FALSE)$counts
    scount = smooth1d(count, err=err, sv2=sv2, verb=FALSE)$fit
    f0 = scount/(length(zstar)*(xbreaks[2]-xbreaks[1]))  

    # do all the auxiliary calculations based on the densities
    f0fz = f0/fz        # strictly aux
    # Only the simple approach is implemented
    if (missing(p0)) {
        p0.est = TRUE        
        p0 = 1/max(f0fz[abs(xbreaks)<zlim])
    } else {
        p0.est = FALSE
    }
    # what we really want
    fdr = p0*f0fz
    
    ifdr = approx(xmids, fdr, xout=Z, rule=2)$y

    # something nice to return
    # HACK: setting the first column name explicitly GAH!  
    ret = data.frame(tstat=Z, fdr.local=ifdr)
    rn = rownames(xdat)
    if (!is.null(rn)) rownames(ret) = rn
    param = list(p0=p0, p0.est=p0.est, fdr=fdr, xbreaks=xbreaks)
    attr(ret, "param") = param
    class(ret) = c("fdr1d.result","fdr.result","data.frame")
    ret
}

##------------------- Plotting --------------------------------------------##
## 130705

plot.fdr1d.result = function(x, add=FALSE, grid=FALSE, rug=TRUE, 
                             xlab="t-Statistic", ylab="fdr", lcol="black", ...) 
{
    p = attr(x, "param")
    xmid = brk2mid(p$xbreaks)
    if (add) {
        lines(xmid, p$fdr, col=lcol, ...)
    } else {
        plot(xmid, p$fdr, col=lcol, type="l", xlab=xlab, ylab=ylab, ...)
    }
    if (rug) {
        rug(x[,1])
    }
    if (grid) {
        segments(xmid, 0, xmid, p$fdr, lty=3, col=lcol)
    }    
}

##------------------- Poisson smoothing -----------------------------------##
## 130705

smooth1d = function(y, sv2=0.1, err=0.01, verb=TRUE)
#
# Name: smooth1d
# Desc: uses a mixed model to smooth a vector of counts (presumably histogram
#       counts)
# Auth: origical code by Yudi.Pawitan@meb.ki.se, see Y.P., In all Likelihood,
#       Oxford Press 2001, chapter 18.11 and
#       http://statistics.ucc.ie/staff/yudi/likelihood/index.htm
#       wrapping & beautifying by Alexander.Ploner@mep.ki.se
#
# Chng: 260903 AP added starting value for bw-parameter lambda
#       020204 AP added smoothing parameter df, removed some bugs
#       130705 AP re-write as smooth1d for frd1d()
#
{
    vcat = function(...) if (verb) cat(...)    
    # Set up auxiliaries: slightly wasteful in the loop
    n    = length(y)
    R1 = Rinv(n, 2)     # second differences
    beta = log(mean(y))
    b    = rep(0, n)
            
    # set up the iteration
    lambda = 1/sv2
    oldval = 2*lambda
    
    while (abs(oldval-lambda)/lambda > err & lambda< 1/err^2)
    {
        oldval = lambda
        eta = beta + b
        # Step 1: update working vector & weights, random parameters        
        p   = exp(eta)
        eps = (y-p)/(p+0.000001)
        Y   = eta + eps
        wt  = p        
        b = solve(diag(wt) + R1/sv2, wt*(Y-beta))
    
        # Step 2: update working vector & weights, fixed parameter
        eta = beta + b
        p   = exp(eta)
        eps = (y-p)/(p+0.000001)
        Y   = eta + eps
        wt  = p        
        beta = sum(wt*(Y - b))/sum(wt)        
    
        # Step 3: update sv2, df
        brb = drop(b%*%R1%*%b)
        # Sic! B symmetric implies trace(A%*%B) = sum(A*B)
        trc = sum(solve(diag(wt) + R1/sv2) * R1)
        sv2 = (brb + trc)/(n-2)
        df = sum(diag(solve((diag(wt) + R1/sv2))) * diag(wt) )
            
        lambda = 1/sv2
        vcat("sv2 =",sv2,"  df =",df,"\n")
    }

    list(fit=p, df=df, sv2=sv2)    
}

Rinv = function(n, d=2)
{
    # differencing matrix
    delta = matrix(0, ncol=n, nrow=n-1)
    for (i in 1:(n-1))
    {
        delta[i,i] = -1
        delta[i,i+1] = 1
    }
    cnt = 1
    while (cnt < d)
    {
        delta = delta[-1,-(1:cnt)]%*% delta
        cnt = cnt+1
    }
    Rinv  = t(delta) %*% delta
    Rinv
}

