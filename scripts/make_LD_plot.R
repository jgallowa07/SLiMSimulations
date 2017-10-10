# This is a script that makes LD plots from a certain set of files.
#
# It can be used as a starting point for a more general script
# that can be applied to other files.

basedir <- "SLiMSimulations/Output1/MyRecipe4_0_5"

base_strings <- c("oldlake100", "newlake100", "fresh1K", "marine1K")

eff <- scan(file.path(basedir, "effectSizes.txt"))

read_ped <- function (base) {
    # first six columns in PED are
    ped_cols <- c("family", "id", "X1", "X2", "X3", "phenotype")
    ped_data <- data.table::fread(file.path(basedir,paste0(base,".ped")))
    pos <- scan(file.path(basedir,paste0(base,"_pos.txt")))
    # matrix of number of Ts
    geno <- c(0,1,1,2)[match(as.matrix(ped_data[,-(1:6)]), c("AA", "AT", "TA", "TT"))]
    dim(geno) <- dim(ped_data) - c(0,6)
    info <- as.data.frame(cbind(ped_data[,1:6], pos))
    names(info) <- c(ped_cols, "geog_loc")
    info <- info[,!grepl("^X", names(info))]
    stopifnot(ncol(geno) == length(eff))
    return(list(geno=geno, info=info))
}

ld_cols <- colorRampPalette(c("blue", "white", "red"))(64)
ld_col_fn <- function (x) { ld_cols[as.numeric(cut(x,breaks=seq(-1,1,length.out=length(ld_cols))))] }

plot_mat <- function (segsites, LD, ...) {
    image(segsites, segsites, LD, zlim=c(-1,1), 
            col=ld_cols, xaxt='n', yaxt='n', ...)
    do_pts <- lower.tri(LD) & !is.na(LD)
    points(segsites[row(LD)[do_pts]], segsites[col(LD)[do_pts]], 
            pch=20, col=ld_col_fn(LD[do_pts]))
}

plot_ld <- function (base1, base2) {
    pop1 <- read_ped(base1)
    pop2 <- read_ped(base2)

    # segregating sites
    nindivs <- c(nrow(pop1$geno), nrow(pop2$geno))
    segsites <- which((colSums(pop1$geno) > 0 & colSums(pop1$geno) < 2*nindivs[1]) | 
                      (colSums(pop2$geno) > 0 & colSums(pop2$geno) < 2*nindivs[2]))
    nsnps <- length(segsites)

    is_fun <- (eff[segsites] != 0)

    # everything
    LD <- list(pop1=cor(pop1$geno[,segsites]),
               pop2=cor(pop2$geno[,segsites]),
               both=cor(rbind(pop1$geno[,segsites], pop2$geno[,segsites])))
    # zero out non-functional LD on the lower triangle
    for (k in seq_along(LD)) {
        LD[[k]][lower.tri(LD[[k]], diag=FALSE) & !outer(is_fun,is_fun,"&")] <- NA
    }

    lims <- c(-0.05,1.05)*ncol(pop1$geno)

    outbase <- paste0(base1, "_", base2)

    png(file=paste0(outbase,"_LD.png"), width=1200, height=1000, pointsize=10, res=1000/6, type='cairo')
    layout(matrix(c(2,1,4,3,5,5), nrow=2), widths=c(1,1,0.2))
    omar <- par("mar")
    par(mar=c(omar,0.1)[c(1,2,5,5)])
    # x <- LD[1:nsnps, 1:nsnps]
    # matplot(segsites[row(x)], segsites[col(x)], pch=20,
    #          xlim=lims, ylim=lims, col=ld_col_fn(x), 
    #          xaxt='n', yaxt='n', xlab=base1, ylab=base1)
    plot_mat(segsites, LD[["pop1"]], xlim=lims, ylim=lims, xlab=base1, ylab=base1)
    rug(segsites, side=1); rug(segsites, side=2)
    par(mar=c(omar,0.1)[c(5,2,3,5)])
    plot_mat(segsites, LD[["both"]], xlim=lims, ylim=lims, xlab="", ylab=base2)
    rug(segsites, side=2); rug(segsites, side=3)
    # par(mar=c(omar,0.1)[c(1,5,5,4)])
    plot(0, type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
    par(mar=c(omar,0.1)[c(5,5,3,4)])
    plot_mat(segsites, LD[["pop2"]], xlim=lims, ylim=lims, xlab=base2, ylab="")
    rug(segsites, side=1); rug(segsites, side=3); rug(segsites, side=4)
    par(mar=c(omar,1)[c(1,5,3,5)])
    plot(rep(0,length(ld_cols)), seq(-1,1,length.out=length(ld_cols)),
            main='cor', xaxt='n', xlab='', ylab='', pch=22, col=ld_cols, bg=ld_cols, cex=3)
    dev.off()

    invisible(return(list(base1=base1, base2=base2, LD=LD)))
}

for (j in seq_along(base_strings)[-length(base_strings)]) {
    for (k in seq(j+1,length(base_strings))) {
        plot_ld(base_strings[j], base_strings[k])
    }
}

