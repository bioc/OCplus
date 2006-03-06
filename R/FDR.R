FDRp = function(xdat, grp, test="t.equalvar", p0, nperm, seed=NULL)
#
# Name: FDRp
# Desc: computes FDR from data via permutations, using package multtest
# Auth: core by YP, wrapper by AP
# Date: 260405
#
# TODO: make F0 more memory efficient
#
# Chng: 260405 AP removed the pval-option from the original code - it is not 
#                 needed for calls from EOC, and it messes up the output object
#       170605 AP changed default seed to NULL for excellent resaons
#                 added the option of doing pairwise t-tests
#       231005 AP made the equal variance t test the default
#
{
    n = length(grp)
    
    require(multtest)
    tstat = mt.teststat(xdat, grp, test=test)
    n = length(tstat)
    ord.abs = order(abs(tstat))   # !!! take absolute value for symmetry
    tstat.o = tstat[ord.abs]
 
    # observed distribution of abs(tstat), folded at 0
    F = 0.5*(1- seq(1/n,1, len=n))
    
    # The permutation for paired t is restricted to flipping 0/1 in adjacent 
    # pairs
    if (test=="pairt") {
        PermGrp = function(g) {
                    n = length(g)
                    for (i in seq(1, n, by=2)) {
                        g[i] = sample(0:1, 1)
                        g[i+1] = 1 - g[i]
                    }
                    g
            }
    } else {
        PermGrp = function(g) sample(g)
    }
 
    # Null distribution of tstat using permutation
    if (!is.null(seed)) set.seed(seed)
    TSTAT = NULL
    for (i in 1:nperm){
        perm.grp = PermGrp(grp)
        TSTAT = c(TSTAT, mt.teststat(xdat, perm.grp, test=test))
    }
    ord = order(abs(TSTAT))
    TSTAT.o = TSTAT[ord]
    # Folded null of abs(tstat) interpolated to the observed values
    F0 = 0.5* approx(abs(TSTAT.o), 1-seq(1/(n*nperm),1, len= n*nperm), 
                     xout=abs(tstat.o), rule=2)$y   # extrapolate
 
    # estimating P(nonDE) = p0 using Efron's method, if required
    if (missing(p0)) {
        a = sum(tstat<0.5 & tstat>-0.5)/n
        b = sum(TSTAT<0.5 & TSTAT>-0.5)/n/nperm
        p0 = min(a/b,1)
    }
 
    # combine rank of the absolute tstat
    crank = rank(c(abs(tstat), abs(TSTAT)))[1:n]
    pval  = 1- (crank - rank(abs(tstat)))/n/nperm
    
    # Compute FDR based on the p-values
    Fp = rank(pval)/length(pval)   # observed distribution of pval
    FDRp = p0 * pmin(pval/Fp,1)

    # this FDR might not be monotone as a function of pvalue: make it so
    ord = order(pval)
    FDR.o = FDRp[ord]            # sort FDR as a function of pval
    b = rev(cummin(rev(FDR.o)))  # find cummin on the tail part
    FDR = rep(0, n)
    FDR[ord] = b  # return monotone function

    # estimating F1 as a function of ordered abs(tstat.o), make it monotone
    F1 = (F-p0*F0)/(1-p0)
    b = rev(cummax(rev(F1)))  # find max of the tail
    sens = pmin(2*b, 1)  
    sens= sens 
    
    # in the original order, not ordered by t.statistic!!!
    ord.orig= rank(abs(tstat))

    list(stat=tstat, pvalue=pval, FDR=FDR, 
                F=F[ord.orig], F0= F0[ord.orig], sens=sens[ord.orig], p0=p0)
}

