##
## fdr2d.R
##
##  Core functions for computing the two-dimensional local false discovery rate:
##      fdr2d
##      approx2d
##      plot.fdr2d.result
##      Volcanoplot
##      Tornadoplot
##      DoConstrain
##      smooth2d.direct
##      smooth2d.basic
##      df2var
##      dgfree
##      rmatrix
##      gausei
##      
##
##  Alexander.Ploner@meb.ki.se 220605
##
##      090805 prepare for OCplus
##      120805 modified fdr2d
##      150805 modified gridding, added smoothing adjustment
##      010905 integarted into OCplus
##
###############################################################################

##----------------------- Computing the FDR2D ------------------------------##
## 090805

fdr2d = function(xdat, grp, test, p0, nperm=100, nr=15, seed=NULL, null=NULL, 
                 constrain = TRUE, smooth=0.2, verb=TRUE, ...)
#
# Name: fdr2d
# Desc: compute two-dimensional local fdr for t-statistics plus se
# Uses: 
# Auth: Alexander.Ploner@meb.ki.se 090805
#       based on earlier code by AP & Yudi.Pawitan@meb.ki.se 
#
# TODO:
#
# Chng: 120805 AP changed smoothing, scaling of fdr, added p0 warning
##      150805 modified gridding, added smoothing adjustment
#       031105 AP changed function name to tstatistics
#
{
    vcat = function(...) if (verb) cat(...)
    
    # The default function
    if (missing(test)) {
        test = function(...) tstatistics(..., logse=TRUE)
    }

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
    # Degree of smoothing all right?
    if ((smooth < 0.01) | (smooth > 0.99)) {
        smooth = min(max(smooth, 0.01), 0.99)
        warning("smooth must be between 0.01 and 0.99 - reset to ",smooth)
    }
    
    # compute the original statistics
    Z = test(xdat, grp, ...)   
    
    # Breaks and such 
    xbreaks = MidBreaks.t(Z[,1], nr)
    ybreaks = MidBreaks(Z[,2], nr)
    xmid = brk2mid(xbreaks)
    ymid = brk2mid(ybreaks)       
    
    # Do the permutations: we are not efficient here in that we store them 
    # all instead of binning them immediately
    if (is.null(null)) {
        zstar = matrix(0, nrow=ng*nperm, ncol=2)
        if (!is.null(seed)) set.seed(seed)
            vcat("Starting permutations...\n")        
        for (i in 1:nperm){
            vcat("\t", i,"\n")
            perm = sample(nc)
            zstar[((i-1)*ng+1):(i*ng), 1:2] = as.matrix(test(xdat, grp[perm]))
        }
    } else { # Use pre-packaged null - not efficient either (double allocation)
        nperm = nrow(null)/ng
        if (nperm != floor(nperm)) 
            stop("Permutations not multiple of number of genes")
        zstar = null
        vcat("Using pre-computed permutations (nperm =",nperm,", seed =",
             attr(null, "seed"),")\n")
    }    

    # Table it
    cut.stat.x = cut(c(Z[,1], zstar[,1]), xbreaks, include.lowest=TRUE)
    cut.stat.y = cut(c(Z[,2], zstar[,2]), ybreaks, include.lowest=TRUE)    
    cut.null.x = cut(zstar[,1], xbreaks, include.lowest=TRUE)
    cut.null.y = cut(zstar[,2], ybreaks, include.lowest=TRUE)    

    tab.stat = table(cut.stat.x, cut.stat.y)
    tab.null = table(cut.null.x, cut.null.y)

    # Smooth & scale
    df.smooth = (1-smooth)*(nr-1)^2
    b = smooth2d.direct(tab.null, tab.stat, df=df.smooth)    
    f0fz = b$fit/((1-b$fit)*nperm)
    
    # An explicit estimate for p0: taking the maximum for a central area, 
    # comparable to 1D
    if (missing(p0)) {
        zndx = findInterval(0,xbreaks)
        nn = length(ybreaks)
        mndx = floor(nn/4):ceiling(3*nn/4)
        p0 = min(1/f0fz[zndx, mndx])
        #p0 = 1/mean(f0fz[zndx, mndx])
        p0.est = TRUE
    } else {
        p0.est = FALSE
    }
    
    # Check how sensible p0 is 
    if (p0 > 1) {
        warning("p0 estimate > 1. This may be due to oversmoothing",
                " - consider reducing smooth")
    }

    # The fdr    
    fdr  = p0*f0fz
    dimnames(fdr) = dimnames(tab.stat)
    
    # Constrain to be monotonous, if required
    if (constrain) 
        fdr = DoConstrain(fdr, xmid)
    
    # Compute the individual approximations
    ifdr = approx2d(xmid, ymid, fdr, Z[,1], Z[,2])
    
    # Output
    ret = data.frame(cbind(Z, ifdr))
    colnames(ret)[3] = "fdr.local"
    param = list(p0=p0, p0.est=p0.est, fdr=fdr, xbreaks=xbreaks, ybreaks=ybreaks) 
    attr(ret, "param") = param
    class(ret) = c("fdr2d.result","fdr.result","data.frame")
    ret
}

##------------------- Helper: approximate linearly on a 2D grid --------------##
##

approx2d = function(x,y, z, xout, yout, ...)
#
# Name: approx2d
# Desc: Compute linear approximation on a grid, reasonably but not crazily 
#       efficient; points oustide the grid are asigned to the closest grid point
# Auth: Alexander.Ploner@meb.ki.se
#
# Chng:
#
{
	n = length(x)
	m = length(y)
	no = length(xout)
	zout = rep(NA, no)
	for (i in 2:n) {
		for (j in 2:m) {	
			ndx = x[i-1]<=xout & xout<=x[i] & y[j-1]<=yout & yout<=y[j]
			lower = z[i-1,j-1] + (xout[ndx]-x[i-1])*(z[i,j-1]-z[i-1,j-1])/(x[i]-x[i-1])
			upper = z[i-1,j]   + (xout[ndx]-x[i-1])*(z[i,j]  -z[i-1,j])/(x[i]-x[i-1])
            zout[ndx] = lower  + (yout[ndx]-y[j-1])*(upper-lower)/(y[j]-y[j-1])
        }
    }
    # For the undefined guys, we take the closest grid point
    undef = which(is.na(zout))
    if (length(undef>0)) {
        gg = as.matrix(expand.grid(x, y))
        for (i in undef) {
            d2 = (xout[i]-gg[,1])^2 + (yout[i]-gg[,2])
            nn = which.min(d2)
            zout[i] = z[nn]
        }
    }
    zout
}
    
##---------------------- Plot a 2D fdr result ------------------------------##
## 220605
## 180805 added akima


plot.fdr2d.result = function(x, levels, nr.plot=20, add=FALSE, grid=FALSE,  
                             pch=".", xlab, ylab, vfont=c("sans serif", "plain"), 
                             lcol="black",  ...)
#
# Name: plot.fdr2d.result
# Desc: a S3 plotting method for fdr2d results
# Uses: akima for pretty smoothing
# Auth: Alexander.Ploner@meb.ki.se 180805
#
# Chng:
#
{
    if (missing(levels)) 
        levels = c(0.05, 0.1, 0.2, 0.3)
    if (missing(xlab))
        xlab = colnames(x)[1]
    if (missing(ylab))
        ylab = colnames(x)[2]
    p = attr(x, "param")
    if (add) {
        points(x[,1], x[,2], pch=pch, ...)
    } else {
        plot(x[,1], x[,2], xlab=xlab, ylab=ylab, pch=pch, ...)
    }        

    cfdr = list() 
    # No smoothing
    if (nr.plot<1) {
        cfdr$x = brk2mid(p$xbreaks)
        cfdr$y = brk2mid(p$ybreaks)
        cfdr$z = p$fdr
    } else {
        # Actually, not using MidBreaks.t and MidBreaks seems to please more
        xo = seq(min(x[,1]), max(x[,1]), len=nr.plot)
        yo = seq(min(x[,2]), max(x[,2]), len=nr.plot)
        cfdr = interp(x[,1], x[,2], x[,3] , xo, yo)
    }
    contour(cfdr, add=TRUE, col=lcol, levels=levels, vfont=vfont)
    if (grid) {
        abline(v=p$xbreaks, lty=3)
        abline(h=p$ybreaks, lty=3)
    }
    invisible(x)
}

##------------------ Plotting on alternative scales ------------------------##
## 090805
## 190805

Tornadoplot = function(x, levels, nr.plot=20, label=FALSE, constrain=FALSE, pch=".", 
                       xlab, ylab, vfont=c("sans serif", "plain"), lcol="black",  ...)
{
    # We check whether we really have a suitable argument
    if (!("fdr2d.result" %in% class(x)))
        stop("requires object of class fdr2d.result")
    cc = colnames(x)
    if (cc[1] != "tstat" | cc[2] != "logse")
        stop("requires fdrd2d based on tstat and logse")
    
    # Some defaults
    if (missing(levels)) 
        levels = c(0.05, 0.1, 0.2, 0.3)
    if (missing(xlab))
        xlab = "Mean Difference"
    if (missing(ylab))
        ylab = "log(Standard Error)"    
    
    # Conversion routines
    conv = function(t, logse) 
    {
        md    = t*exp(logse)
        cbind(md, logse)
    }
    lconv = function(x) 
    {
        xx = conv(x$x, x$y)
        list(level=x$level, x=xx[,1], y=xx[,2])
    }
        
    # Now go
    cc = conv(x[,1], x[,2])
    md = cc[,1]
    logse  = cc[,2]
    
    # Compute the contour lines in the original scale
    cfdr = list() 
    if (nr.plot<1) {
        p = attr(x, "param")
        cfdr$x = brk2mid(p$xbreaks)
        cfdr$y = brk2mid(p$ybreaks)
        cfdr$z = p$fdr
    } else {
        # Actually, not using MidBreaks.t and MidBreaks seems to please more
        xo = seq(min(x[,1]), max(x[,1]), len=nr.plot)
        yo = seq(min(x[,2]), max(x[,2]), len=nr.plot)
        cfdr = interp(x[,1], x[,2], x[,3] , xo, yo)
    }
    if (constrain) {
        cfdr$z = DoConstrain(cfdr$z, cfdr$x)
    }
    cl = contourLines(cfdr, levels=levels)
    cl.conv = lapply(cl, lconv)

    # Do it    
    plot(md, logse, xlab=xlab, ylab=ylab, pch=pch, ...)
    DrawContourlines(cl.conv, label=label, col=lcol, vfont=vfont)
    invisible(x)    
}

Volcanoplot = function(x, df, levels, nr.plot=20, label=FALSE, constrain=FALSE, 
                       pch=".", xlab, ylab, vfont=c("sans serif", "plain"), 
                       lcol="black",  ...)
{
    # We check whether we really have a suitable argument
    if (!("fdr2d.result" %in% class(x)))
        stop("requires object of class fdr2d.result")
    cc = colnames(x)
    if (cc[1] != "tstat" | cc[2] != "logse")
        stop("requires fdrd2d based on tstat and logse")

    # Some defaults
    if (missing(levels)) 
        levels = c(0.05, 0.1, 0.2, 0.3)
    if (missing(xlab))
        xlab = "Mean Difference"
    if (missing(ylab))
        ylab = "-log10(p-Value)"
    
    # Conversion routine
    conv = function(t, logse) 
    {
        md    = t*exp(logse)
        logp  = -log10(2*pt(-abs(t), df=df))
        logp[logp==0] = min(logp[logp!=0])
        cbind(md, logp)
    }
    lconv = function(x) 
    {
        xx = conv(x$x, x$y)
        list(level=x$level, x=xx[,1], y=xx[,2])
    }
        
    # Now go
    cc = conv(x[,1], x[,2])
    md = cc[,1]
    logp  = cc[,2]
    
    # Compute the contour lines in the original scale
    cfdr = list() 
    if (nr.plot<1) {
        p = attr(x, "param")
        cfdr$x = brk2mid(p$xbreaks)
        cfdr$y = brk2mid(p$ybreaks)
        cfdr$z = p$fdr
    } else {
        # Actually, not using MidBreaks.t and MidBreaks seems to please more
        xo = seq(min(x[,1]), max(x[,1]), len=nr.plot)
        yo = seq(min(x[,2]), max(x[,2]), len=nr.plot)
        cfdr = interp(x[,1], x[,2], x[,3] , xo, yo)
    }
    if (constrain) {
        cfdr$z = DoConstrain(cfdr$z, cfdr$x)
    }
    cl = contourLines(cfdr, levels=levels)
    cl.conv = lapply(cl, lconv)
        
    plot(md, logp, xlab=xlab, ylab=ylab, pch=".", ...)
    DrawContourlines(cl.conv, label=label, col=lcol, vfont=vfont)
    invisible(x)    
}

# Another helper function: constrain
DoConstrain = function(fdr, xmid)
{
    rhs = xmid >= 0
    lhs = xmid < 0
    for (i in 1:ncol(fdr)){
        aa = fdr[rhs,i]
        fdr[rhs,i] = cummin(aa)
        bb = rev(fdr[lhs,i])
        fdr[lhs,i] = rev(cummin(bb))
    }
    fdr
}

###############################################################################
## De smoothing section: all writen by Yudi.Pawitan@meb.ki.se, but 
##   prettified by Alexander.Ploner@meb.ki.se
###############################################################################


smooth2d.direct = function(y, n, df = 100, err=0.01, edge.count=3, mid.start=6)
{
    # Technically, we can allow different nx/ny, but note that we use the 
    # same edge/middle definition for both!
    nx = nrow(y)
    ny = ncol(y)
    # Testing: we can tolerate pretty small nx, becausd the indices 
    # will run backward
    if (nx < mid.start | edge.count > ny)
        stop("Grid size too small")
    # Implicit imputation around the corners: all zero counts are set to 
    # fdr zero
    n = ifelse(n <= 1,1,n)
    # Modify this in the lower middle: here the fdr is set to zero 
    edge = 1:edge.count
    mid  = mid.start:(nx-mid.start+1)
    y[mid,edge]= ifelse(n[mid,edge]==1, 1, y[mid,edge])

    p = y/n
    sv2 = df2var(nx,df=df)
    fit = smooth2d.basic(p,sv2=sv2,err=err) 
    list(fit=fit, df=df, sv2=sv2)
}

# tranform df to sv2
df2var = function(nx,df=nx, lower=0.0001, upper=1000)
{
    ff = function(x) dgfree(nx,x)-df
    out = uniroot(ff, lower=lower, upper=upper)
    return(out$root)
}

# degrees of freedom of basic smooth
dgfree = function(nx, sv2=1, ny=nx)
{
    imp = matrix(0,nrow=nx, ncol=ny)             ## impulse function
    imp[nx/2, ny/2] = 1
    imp.resp = smooth2d.basic(imp,sv2=sv2)  ## impulse response
    df = nx*ny*imp.resp[nx/2, ny/2]
    return(df)
}


# y is data matrix on regular grid
smooth2d.basic = function(y, sv2=0.01,se2=1,err=0.01)
{
    W = 1/se2
    nx = nrow(y)
    ny = ncol(y)
    Y = c(y)    # vectorized
    R = rmatrix(nx,ny)
    d =  R$d/sv2 + W   # diagonal elements
    xx =  R$x  # neighbors index

    ybar = mean(Y)
    B =  W*(Y-ybar)    # Z'W (Y-Xb) 
    v =  gausei(d,xx,B,sv2, err=err, maxiter=100)
    est = ybar+v
    matrix(est,ncol=ny)
}

# Index of neighbors of regular nx.ny grids, maximum nx and ny is 1000
rmatrix= function(nx,ny)
{
    xymat = matrix(0,nrow=nx, ncol=ny)
    grid = c((row(xymat)-1)*1000 + col(xymat)-1)
    neighbor= c(grid+1000,grid-1000,grid+1,grid-1)
    loc = match(neighbor,grid,0)
    loc = matrix(loc,ncol=4)
    num = c((loc>0)%*% rep(1,4))
    list(x=loc,d=num)
}

#/***************  Gauss Seidel Function  **********************/
# solving A v =B
# d contains diagonal elements of A
# x contains index of neighbours
# sv2 = sigma2v 
gausei=function(d,xx,B,sv2, err=0.01, maxiter=10)
{
    ng = length(d)                      # number of grid points
    xx[xx<1] = ng+1                     # non-neighbors = last location 
    v = c(B/(d+0.000001),0)       # rough solution, use 0 for non-neighbors
                                  # at last location
    vold = 3*v
    iter = 1
    while ( (sd(vold-v, na.rm=TRUE)/sd(v, na.rm=TRUE)  > err) & (iter< maxiter) )
    {
        vold = v
        vmat =  matrix(v[xx],ncol=4)
        v[1:ng]  = (B + c(vmat%*%rep(1,4))/sv2) /d
        iter = iter+1
    }   # end while
    v[1:ng]
}

