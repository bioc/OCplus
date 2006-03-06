##
## Simulations.R
##
##  A collection of routines useful for simulating microarray expression data 
##
##      MAsim: fixed DE, fixed variance
##      MAsim.var: fixed DE, random variance
##      MAsim.smyth: Smyth's model for random DE, random variance
##      MAsim.real: bootstrap from real data, fixed but scaled DE
##
## Alexander.Ploner@meb.ki.se 241005
###############################################################################

## Note: these guys use attr(x, "DE") to store the DE status, and the 
## colnames as group labels

#.................................................. Fixed effect, fixed variance

MAsim = function(ng=10000, n=10, n1=n, n2=n, D=1, p0=0.9, sigma=1) 
#
# Name: MAsim
# Desc: simulate a two-sample microarray data set with ng genes, n 
#       samples in each group, an effect size of D, and a proportion p0 of
#       non differentially expressed genes, using constant variance
# Auth: Alexander.Ploner@meb.ki.se
#
# Chng: 031105 AP added unequal group sizes
#
{
    nn = n1+n2
    group = rep(c(0,1),c(n1,n2))
    xdat = matrix(rnorm(nn*ng, mean=0, sd=sigma), nrow=ng, ncol=nn)
    fc = rep(0, ng)
    ran1 = runif(ng) > p0
    ran2 = runif(ng) > 0.5
    fc = ifelse(ran1 & ran2, -D, fc)
    fc = ifelse(ran1 & !ran2, D, fc)
    # Note: works because matrix is stored columnwise
    xdat[,group==1] = xdat[,group==1] + fc*sigma
    colnames(xdat) = as.character(group)
    attr(xdat, "DE") = ran1
    xdat
}

#............................................... Fixed effect, variable variance

MAsim.var = function(ng=10000, n=10, n1=n, n2=n, D=1, p0=0.9) 
#
# Name: MAsim.var
# Desc: simulate a two-sample microarray data set with ng genes, n 
#       samples in each group, an effect size of D, and a proportion p0 of
#       non differentially expressed genes, using variable variances
# Auth: Alexander.Ploner@meb.ki.se, following Yudi's suggestion
#
# Chng: 031105 AP added unequal group sizes
#
{
    sigma = sqrt(1/rexp(ng))
    nn = n1+n2
    group = rep(c(0,1),c(n1,n2))
    xdat = matrix(rnorm(nn*ng, mean=0, sd=sigma), nrow=ng, ncol=nn)
    fc = rep(0, ng)
    ran1 = runif(ng) > p0
    ran2 = runif(ng) > 0.5
    fc = ifelse(ran1 & ran2, -D, fc)
    fc = ifelse(ran1 & !ran2, D, fc)
    # Note: works because matrix is stored columnwise
    xdat[,group==1] = xdat[,group==1] + fc*sigma
    colnames(xdat) = as.character(group)
    attr(xdat, "DE") = ran1
    xdat
}

#......................................................................Smyth

MAsim.smyth = function(ng=10000, n=10, n1=n, n2=n, p0=0.9, d0=4, s2_0=4, v0=2) 
#
# Name: MAsim.smyth
# Desc: simulate a balanced two-sample microarray data set with ng genes, and n 
#       samples in each group; the distribution follows a mixture of normals,
#       see Smyth, Statistical Applications in Genetics & Molecular Biology, 2004
# Note: Smyth does not really simulate data - only estimators; we are more 
#       specific
# Auth: Alexander.Ploner@meb.ki.se 141005
#
# Chng: 031105 AP added unequal group sizes
#
{
    # Set up grouping and data matrix
    nn = n1+n2
    group = rep(c(0,1),c(n1,n2))
    xdat = matrix(0, nrow=ng, ncol=nn)
    # Generate the gene-wise variances
    s2_g = d0*s2_0/rchisq(ng, df=d0)
    # Generate the indices of the DE genes
    ndx = runif(ng) > p0
    nde = length(which(ndx))
    # Generate the effect sizes for the DE genes
    def = rnorm(nde, mean=0, sd=sqrt(v0*s2_g))
    # Fill in the effects 
    xdat[ndx, group==1] = matrix(def, nrow=nde, ncol=n2)
    # Now simulate the errors: This is somewhat fishy
    xdat = matrix(rnorm(nn*ng, mean=xdat,sd=sqrt(s2_g)), nrow=ng) 
    # Annotate the result
    colnames(xdat) = as.character(group)
    des = rep(FALSE, ng)
    des[ndx] = TRUE
    attr(xdat,"DE") = des
    xdat
}


#................................. simulating by sampling from real data,
#..................................imputing the effect size

MAsim.real = function(xdat, grp, n, n1, n2, D=1, p0=0.9, replace=TRUE)
#
# Name: MAsim.real
# Desc: simulate from a real micorarray data set, using either sub-sampling
#       or a bootstrap sample from the groupwise residuals
# Auth: Yudi.Pawitan$meb.ki.se
#
# Chng: 031105 AP 
#           added unequal group sizes
#           removed major inconsistency (sampling before/after residuals)
#
{
    # Check group parameters
    glab = unique(grp)
    if (length(glab)!=2) {
        stop("Need exactly two different groups")
    }
    ndx  = grp == glab[1]
    n1.dat = table(ndx)["TRUE"]
    n2.dat = table(ndx)["FALSE"]
    # Check/impute the size of the simulated data set
    if (missing(n1)) {
        n1 = if (missing(n)) n1.dat else n
    }
    if (missing(n2)) {
        n2 = if (missing(n)) n2.dat else n
    }
    if (!replace & (n1>n1.dat | n2>n2.dat))
        stop("Cannot sub-sample for values of n1, n2 - reduce or set replace=TRUE")

    # Sub-sample or bootstrap within group, as required    
    id1 = sample(which(ndx), n1, replace=replace)    
    id2 = sample(which(!ndx), n2, replace=replace)
    # Replace the original data by the resampled version
    xdat = cbind(xdat[ , id1], xdat[, id2])
    grp  = rep(glab, c(n1, n2))
    ndx  = grp==glab[1]
    # Now compute the groupwise residuals
    m1 = rowMeans(xdat[,ndx])
    m2 = rowMeans(xdat[,!ndx])
    xdat[ ,ndx]  = xdat[,ndx] - m1
    xdat[ ,!ndx] = xdat[,!ndx] - m2
    msd = sqrt(rowSums(xdat*xdat)/(n1+n2-2))        

    # Effect is scaled by observed standard deviation
    shift = D* msd      # DE
    ng = nrow(xdat)
    fc = rep(0, ng)   # unit fold change
    ran1 = runif(ng) >p0
    ran2 = runif(ng) >0.5
    fc = ifelse(ran1 & ran2, -1, fc)
    fc = ifelse(ran1 & !ran2, 1, fc)
    fc = fc*shift  

    dat1 = xdat[,ndx] + rt(ng,df=n1-1)*msd/sqrt(n1)
    dat2 = xdat[,!ndx] + rt(ng,df=n2-1)*msd/sqrt(n2) + fc 
    ret = cbind(dat1, dat2)
    colnames(ret) = rep(glab, c(n1, n2))
    attr(ret, "DE") = ran1
    ret
}




