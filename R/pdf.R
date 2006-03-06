## pdf.R
##
## Helper functions for calculating densities and local false discoveries for 
##   mixtures of t-distributions, see also CDF.R
##
##      dmt
##      dmtmix
##      dmt.paired
##      dmtmix.paired
##      lfdr
##      lfdr.paired
##
## Alexander.Ploner@meb.ki.se 301005
##############################################################################


dmt = function (x, n1, n2, D, p, sigma) 
#
# Name: dmt
# Desc: density function for a mixture of t-distributions
#       useful for computation of local fdr and sample sizes
# Auth: AP, based on CDF 301005
#
# Note: * we assume clean parameters, ie sum(p)==1, length(p)==length(D) etc.
#       * wankerish recursion, so sue me
#
# Chng: 
#
{
    if (length(D) == 0) 
        return(0)
    ncp = -D[1]*sqrt((n1*n2)/(n1+n2))/sigma
    ret = p[1] * dt(x, df = n1+n2 - 2, ncp = ncp) + dmt(x, n1, n2, D[-1], p[-1], sigma)
    ret

}

dmtmix = function (x, n1, n2, pmix, D0, p0, D1, p1, sigma) 
#
# Name: dmtmix
# Desc: density of the mixture of f0 and f1, both of which are mixtures of 
#       t-distributions
#       useful for computation of local fdr and sample sizes
# Auth: AP, based on CDFmix 301005
#
# Note: we assume clean parameters, ie sum(p)==1, length(p)==length(D) etc.
#
# Chng: 
#
{
    pmix * dmt(x, n1, n2, D0, p0, sigma) + (1 - pmix) * dmt(x, n1, n2, D1, p1, sigma)
}

lfdr = function (x, n1, n2, pmix, D0, p0, D1, p1, sigma)
#
# Name: lfdr
# Desc: local false discovery rate - a ratio of densities, which are both 
#       mixtures of t-distributions
#       Note that this is for symmetric cutoff, ie essentially for |t|
# Auth: AP, based on FDR 301005
#
# Note: we assume clean parameters, ie sum(p)==1, length(p)==length(D) etc.
#
# Chng: 
#
{
    # It's all in the tail! See Efron, large-scale simultaneous hypothesis
    # testing, Remark B:
    x = -abs(x)
    pmix * (dmt(x, n1, n2, D0, p0, sigma) + dmt(-x, n1, n2, D0, p0, sigma))/
           (dmtmix(x, n1, n2, pmix, D0, p0, D1, p1, sigma) + 
            dmtmix(-x, n1, n2, pmix, D0, p0, D1, p1, sigma))
}

"dmt.paired" = function (x, n, D, p, sigma) 
#
# Name: dmt.paired
# Desc: density function for a mixture of t-distributions
#       useful for computation of local fdr and sample sizes for paired t-tests
# Auth: AP, based on CDF.paired 301005
#
# Chng: 
#
{
    if (length(D) == 0) 
        return(0)
    ncp = -D[1]*sqrt(n)/sigma
    ret = p[1] * dt(x, df = n - 1, ncp = ncp) + dmt.paired(x, n, D[-1], p[-1], sigma)
    ret

}

"dmtmix.paired" = function (x, n, pmix, D0, p0, D1, p1, sigma) 
#
# Name: dmtmix.paired
# Desc: density function for a mixture of f0 and f1, both
#       of which are mixtures of t-distributions
#       useful for computation of local fdr and sample sizes for paired t-tests
# Auth: AP, based on CDFmix.paired 301005
#
# Note: 
#
# Chng: 
#
{
    pmix * dmt.paired(x, n, D0, p0, sigma) + (1 - pmix) * dmt.paired(x, n, D1, p1, sigma)
}

"lfdr.paired" = function (x, n, pmix, D0, p0, D1, p1, sigma)
#
# Name: lfdr.paired
# Desc: local false discovery rate when both f0 and f1 are mixtures of 
#       t-distributions
#       Note that we are considering the lfdr of |t|, which looks a bit weird
# Auth: AP, based on FDR.paired 301005
#
# Note: 
#
# Chng: 
#
{
    # It's all in the tail! See Efron, large-scale simultaneous hypothesis
    # testing, Remark B:
    x = -abs(x)
    pmix * (dmt.paired(x, n, D0, p0, sigma) + dmt.paired(-x, n, D0, p0, sigma))/
           (dmtmix.paired(x, n, pmix, D0, p0, D1, p1, sigma) + 
            dmtmix.paired(-x, n, pmix, D0, p0, D1, p1, sigma))
}




