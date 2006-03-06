##
##  Routines common to fdr1d and fdr2d:
##      tstatistics
##      PermNull
##      summary.fdr.result
##      p0
##      topDE
##      average.fdr
##      MidBreaks.t
##      Midbreaks
##      brk2mid
##      DrawContourlines
##
##  Alexander.Ploner@meb.ki.se  310805
##
###############################################################################

tstatistics = function(xdat, grp, logse=FALSE, paired=FALSE)
#
# Name: tstatistic
# Desc: computes multiple parallel t-statistics in pure R in reasonable time
# Auth: Alexander.Ploner@meb.ki.se
# Note: returns a data.frame for good reasons
#
# Chng: 020905 adapted for OCplus
#       031105 AP changed name to tstatistics
#
{
    glab = unique(grp)
    if (length(glab)!=2) {
        stop("Need exactly two different groups")
    }
    
    if (!paired) {
        ndx  = grp == glab[1]
        n1 = table(ndx)["TRUE"]
        n2 = table(ndx)["FALSE"]
        m1 = rowMeans(xdat[,ndx])
        m2 = rowMeans(xdat[,!ndx])
        xdat[ ,ndx]  = xdat[,ndx] - m1
        xdat[ ,!ndx] = xdat[,!ndx] - m2
        pvar = rowSums(xdat*xdat)/(n1+n2-2)        
        se = sqrt(pvar*(1/n1 + 1/n2))
        mn = m1-m2
        ret = mn/se
    } else {
        n2 = ncol(xdat)
        n  = floor(0.5*n2)
        if (n2/2 != n) {
            stop("Paired t test requires an even number of observations")
        }
        ondx = seq(1, n2, by=2)              # odd ndx
        endx = seq(2, n2, by=2)              # even ndx 
        sodd = (-1)^as.numeric(grp[ondx]!=glab[1])  # sign if t=glab[1]-glab[2]
        sevn = (-1)^as.numeric(grp[endx]!=glab[2])
        if (!all(sodd==sevn)) {
            stop("paired observations must be consecutive")
        }
        d = (xdat[, ondx] - xdat[, endx]) %*% diag(sodd)
        m = rowMeans(d)
        d = d - m
        pvar = rowSums(d*d)/(n-1)
        se = sqrt(pvar/n)
        mn = m
        ret = m/se
    }
    
    if (logse) {
        data.frame(tstat=ret, logse=log(se))
    } else {
        data.frame(tstat=ret)
    }
}

##------------------------------ Permutation null --------------------------##
## 220605

PermNull = function(xdat, grp, nperm=100, seed=NULL, logse=FALSE, paired=FALSE)
{
    ng = nrow(xdat)
    nc = ncol(xdat)
    if (logse) {
        ret = matrix(0, nrow=ng*nperm, ncol=2)
        colnames(ret) = c("tstat","logse")
    } else {
        ret = matrix(0, nrow=ng*nperm, ncol=1)
        colnames(ret) = c("tstat")        
    }        
    if (!is.null(seed)) set.seed(seed)
    for (i in 1:nperm){
        cat(i,"\n")
        perm = sample(nc)
        ret[((i-1)*ng+1):(i*ng),] = as.matrix(tstatistics(xdat, grp[perm], 
                                                    logse=logse, paired=paired))
    }
    attr(ret,"seed") = seed
    as.data.frame(ret)
}


##----------------- generic fdr.result -------------------##
## 090805
summary.fdr.result = function(object, ...) 
{
    ret = list()
    ret$Statistic = summary(object[,1])
    jj1 = cut(object$fdr.local, c(0,0.05, 0.1, 0.2, 1, Inf))
    jj2 = factor(as.numeric(object[,1] < 0), levels=0:1, labels=c("t<0","t>=0")) 
    ret$fdr = table(statistic=jj2, fdr=jj1)
    ret$p0  = p0(object, how=TRUE) 
    ret
}

p0 = function(x, how=FALSE)
{
    if (!("fdr.result" %in% class(x)))
        stop("x must be of class fdr.result")
    p = attr(x, "param")
    if (how)
        ret = list(Value=p$p0,Estimated=p$p0.est)
    else 
        ret = p$p0

    ret
}

topDE = function(x, co=0.1)
#
# Name: topDE
# Desc: extracts a table with the top-regulated genes from the 
#       output of EOC, fdr1d, fdr2d
# Auth: Alexander.Ploner@ki.se 180206
# 
# Chng:
#
{
	if (inherits(x, "FDR.result")){
		nn = 3
	} else if (inherits(x,"fdr1d.result")){
		nn = 2
	} else if (inherits(x,"fdr2d.result")){
		nn = 3
	} 
	
	ndx = x[,nn] <= co
	ret = x[ndx,]
	ret = ret[order(ret[,nn]),]
	ret
}


##--------------------- Average 2d fdr ------------------------------------##
## 150805

# Average the local fdr
average.fdr = function(x, breaks)
{
    # We check whether we really have a suitable argument
    if (!("fdr.result" %in% class(x)))
        stop("requires object of class fdr.result")
    cc = colnames(x)[1]
    if (cc[1] != "tstat")
        stop("requires fdrd based on tstat")
    t  = x$tstat
    f2 = x$fdr.local
    
    if (missing(breaks)) {
        breaks = attr(x, "param")$xbreaks
    }
    tcut = cut(t, breaks)
    avef = tapply(f2, tcut, mean, na.rm=TRUE)
    avet = tapply(t, tcut, mean, na.rm=TRUE)

    cbind(tstat=avet, fdr.local=avef)
}

##------------------------ Helpers: Breaks and Mids -------------------------##
## 100805

MidBreaks.t = function(x, nr, eps=0.01, iter=20)
#
# Name: Midbreaks.t
# Desc: finds a given number of equidistant breaks for a vector of 
#       t-statistics so that one category is centered on zero - 
#       fairly elaborate
# Auth: Alexander.Ploner@meb.ki.se 150805
# 
{
    # We need pos & neg values
    rr = range(x, na.rm=TRUE)
    if (prod(rr) >= 0)
        stop("Need both negative and positive t-statistics - what's going on?!")
    # Iterate to find any solution
    lower = eps
    upper = diff(rr)+2*eps
    i     = 1
    cur.cnt = nr-3
    while(i < iter & cur.cnt != nr-2) {
         cur.d = (lower+upper)/2
         left  = ceiling((-rr[1]+eps)/cur.d-0.5)
         right = ceiling((rr[2]+eps)/cur.d-0.5)
         cur.cnt = left + right 
         if (cur.cnt > nr-2) {
             lower = cur.d
         } else if (cur.cnt < nr-2) {
             upper = cur.d
         }
         i = i+1
    }
    # Now iterate to make the solution small
    cur.work = cur.d    
    sub.good = 0
    sub.work = eps
    i   = 1
    while(i < iter) {
         cur.d = cur.work - sub.work
         left  = ceiling((-rr[1]+eps)/cur.d-0.5)
         right = ceiling((rr[2]+eps)/cur.d-0.5)
         cur.cnt = left + right 
         if (cur.cnt == nr-2) {
             sub.good = sub.work
             sub.work = 2*sub.work
         } else {
             sub.work = (sub.work+sub.good)/2
         }
         i = i+1
    }
    # Do it
    cur.d = cur.work - sub.good
    left  = ceiling((-rr[1]+eps)/cur.d-0.5)
    right = ceiling((rr[2]+eps)/cur.d-0.5)
    # Finish
    ret.left  = -cur.d/2 + seq(from=0, by=-cur.d, len=left+1)
    ret.right =  cur.d/2 + seq(from=0, by=cur.d,  len=right+1)
    ret = c(ret.left, ret.right)
    sort(ret)
}    
    
MidBreaks = function(x, nr)
{
    mids   = seq(min(x), max(x), length=nr-1)
    d = mids[2]-mids[1]
    breaks = c(mids-d/2, mids[nr-1]+d/2)
    breaks
}

# Another popular helper
brk2mid = function(breaks) {breaks[-1]-diff(breaks)/2}

DrawContourlines = function(x, label=FALSE, cex=0.7, 
                            vfont=c("sans serif","bold"), ...)
#
# Name: DrawContourlines
# Desc: this takes output from contourLines and plots them; labelling is 
#       done if required, but note that is just slapdash hack to have SOME
#       kind of labels for the Volcano- and Tornadoplot functions
# Auth: Alexander.Ploner@meb.ki.se 190805 (and not proud of it)
#
{
    # Some prep steps, if required
    n  = length(x)
    # Testing only
    # Set up a plot, if required    
    #if (!add) {
    #    xrange = lapply(x, function(x) range(x$x, na.rm=TRUE))
    #    yrange = lapply(x, function(x) range(x$y, na.rm=TRUE))
    #    xrange = matrix(unlist(xrange),nrow=n, byrow=TRUE)
    #    yrange = matrix(unlist(yrange),nrow=n, byrow=TRUE)
    #    xmin = min(xrange[,1])
    #    xmax = max(xrange[,2])
    #    ymin = min(yrange[,1])
    #    ymax = max(yrange[,2])
    #    plot(c(xmin,xmax), c(ymin, ymax), type="n", ...)
    #}
    # Draw the lines (always, duh) 
    lapply(x, lines, ...)

    if (label) {
        # Tired of lpply, let's loop
        for (i in 1:n) {
            xc = x[[i]]$x
            yc = x[[i]]$y
            ll = x[[i]]$level
            if (max(xc) < 0) { # negative, we go leftmost
                j = which.min(xc)
                text(xc[j], yc[j], ll, pos=2, vfont=vfont, cex=cex)
            } else if (min(xc) > 0) { # positive, we go rightmost
                j = which.max(xc)
                text(xc[j], yc[j], ll, pos=4, vfont=vfont, cex=cex)
            } else { # hm, why not topmost?
                j = which.max(yc)
                text(xc[j], yc[j], ll, pos=3, vfont=vfont, cex=cex)
            }
        }
    }
                
}

