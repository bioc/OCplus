"CDF" = function (x, n1, n2, D, p, sigma) 
#
# Name: CDF
# Desc: a cumulative distribution function which is a mixture of t-distributions
#       useful for computation of FDR and sample sizes
# Auth: AP, based on code by YP 291104
#
# Note: * we assume clean parameters, ie sum(p)==1, length(p)==length(D) etc.
#       * wankerish recursion, so sue me
#
# Chng: 170605 AP adjusted for unequal group sizes
#
{
    if (length(D) == 0) 
        return(0)
    ncp = -D[1]*sqrt((n1*n2)/(n1+n2))/sigma
    ret = p[1] * pt(x, df = n1+n2 - 2, ncp = ncp) + CDF(x, n1, n2, D[-1], p[-1], sigma)
    ret

}

"CDFmix" = function (x, n1, n2, pmix, D0, p0, D1, p1, sigma) 
#
# Name: CDFmix
# Desc: a cumulative distribution function which is a mixture of F0 and F1, both
#       of which are mixtures of t-distributions
#       useful for computation of FDR and sample sizes
# Auth: AP, based on code by YP 291104
#
# Note: we assume clean parameters, ie sum(p)==1, length(p)==length(D) etc.
#
# Chng: 170605 AP adjusted for unequal group sizes
#
{
    pmix * CDF(x, n1, n2, D0, p0, sigma) + (1 - pmix) * CDF(x, n1, n2, D1, p1, sigma)
}

"FDR" = function (x, n1, n2, pmix, D0, p0, D1, p1, sigma)
#
# Name: FDR
# Desc: a cumulative distribution function which is a mixture of F0 and F1, both
#       of which are mixtures of t-distributions
#       useful for computation of FDR and sample sizes
# Auth: AP, based on code by YP 291104
#
# Note: we assume clean parameters, ie sum(p)==1, length(p)==length(D) etc.
#
# Chng: 190405 AP changed name to FDR (was: FDR2)
#       170605 AP adjusted for unequal group sizes
#
{
    # It's all in the tail! See Efron, large-scale simultaneous hypothesis
    # testing, Remark B:
    x = -abs(x)
    pmix * (CDF(x, n1, n2, D0, p0, sigma) + 1 - CDF(-x, n1, n2, D0, p0, sigma))/
           (CDFmix(x, n1, n2, pmix, D0, p0, D1, p1, sigma) + 
            1 - CDFmix(-x, n1, n2, pmix, D0, p0, D1, p1, sigma))
}

"CDF.paired" = function (x, n, D, p, sigma) 
#
# Name: CDF.paired
# Desc: a cumulative distribution function which is a mixture of t-distributions
#       useful for computation of FDR and sample sizes for paired t-tests
# Auth: AP, based on CDF 170605
#
# Chng: 
#
{
    if (length(D) == 0) 
        return(0)
    ncp = -D[1]*sqrt(n)/sigma
    ret = p[1] * pt(x, df = n - 1, ncp = ncp) + CDF.paired(x, n, D[-1], p[-1], sigma)
    ret

}

"CDFmix.paired" = function (x, n, pmix, D0, p0, D1, p1, sigma) 
#
# Name: CDFmix.paired
# Desc: a cumulative distribution function which is a mixture of F0 and F1, both
#       of which are mixtures of t-distributions
#       useful for computation of FDR and sample sizes for paired t-tests
# Auth: AP, based on CDFmix 170605
#
# Note: 
#
# Chng: 
#
{
    pmix * CDF.paired(x, n, D0, p0, sigma) + (1 - pmix) * CDF.paired(x, n, D1, p1, sigma)
}

"FDR.paired" = function (x, n, pmix, D0, p0, D1, p1, sigma)
#
# Name: FDR.paired
# Desc: a cumulative distribution function which is a mixture of F0 and F1, both
#       of which are mixtures of t-distributions
#       useful for computation of FDR and sample sizes for paired t-tests
# Auth: AP, based on FDR 170605
#
# Note: 
#
# Chng: 
#
{
    # It's all in the tail! See Efron, large-scale simultaneous hypothesis
    # testing, Remark B:
    x = -abs(x)
    pmix * (CDF.paired(x, n, D0, p0, sigma) + 1 - CDF.paired(-x, n, D0, p0, sigma))/
           (CDFmix.paired(x, n, pmix, D0, p0, D1, p1, sigma) + 
            1 - CDFmix.paired(-x, n, pmix, D0, p0, D1, p1, sigma))
}




