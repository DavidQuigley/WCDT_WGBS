# Plots Figure 1B

heatmap.3 = function(x,
                     showHeatmap=TRUE,
                     Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                     distfun = dist,
                     hclustfun = hclust,
                     dendrogram = c("both", "row", "column", "none"),
                     symm = FALSE,
                     scale = c("none","row", "column"),
                     na.rm = TRUE,
                     revC = identical(Colv,"Rowv"),
                     add.expr,
                     breaks,
                     symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                     col = "heat.colors",
                     colsep,
                     rowsep,
                     sepcolor = "white",
                     sepwidth = c(0.05, 0.05),
                     cellnote,
                     notecex = 1,
                     notecol = "cyan",
                     na.color = par("bg"),
                     trace = c("none", "column","row", "both"),
                     tracecol = "cyan",
                     hline = median(breaks),
                     vline = median(breaks),
                     linecol = tracecol,
                     margins = c(5,5),
                     ColSideColors,
                     RowSideColors,
                     side.height.fraction=0.3,
                     cexRow = 0.2 + 1/log10(nr),
                     cexCol = 0.2 + 1/log10(nc),
                     labRow = NULL,
                     labCol = NULL,
                     key = TRUE,
                     keysize = 1.5,
                     density.info = c("none", "histogram", "density"),
                     denscol = tracecol,
                     symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                     densadj = 0.25,
                     main = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     lmat = NULL,
                     lhei = NULL,
                     lwid = NULL,
                     ColSideColorsSize = 1,
                     RowSideColorsSize = 1,
                     KeyValueName="Value",...){
    
    invalid = function (x) {
        if (missing(x) || is.null(x) || length(x) == 0)
            return(TRUE)
        if (is.list(x))
            return(all(sapply(x, invalid)))
        else if (is.vector(x))
            return(all(is.na(x)))
        else return(FALSE)
    }
    
    x = as.matrix(x)
    scale01 = function(x, low = min(x), high = max(x)) {
        x = (x - low)/(high - low)
        x
    }
    retval = list()
    if (symm && missing(scale) ){
        scale = "none"
    }else{
        scale = match.arg(scale)
    }
    dendrogram = match.arg(dendrogram)
    trace = match.arg(trace)
    density.info = match.arg(density.info)
    if (length(col) == 1 && is.character(col) ){
        col = get(col, mode = "function")
    }
    if (!missing(breaks) && (scale != "none")){
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
                "specified can produce unpredictable results.", "Please consider using only one or the other.")
    }
    if (is.null(Rowv) || is.na(Rowv)){
        Rowv = FALSE
    }
    if (is.null(Colv) || is.na(Colv)){
        Colv = FALSE
    }
    else if (Colv == "Rowv" && !isTRUE(Rowv)){
        Colv = FALSE
    }
    di = dim(x)
    if (length(di) != 2 || !is.numeric(x)){
        stop("`x' must be a numeric matrix")
    }
    nr = di[1]
    nc = di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2){
        stop("`margins' must be a numeric vector of length 2")
    }
    if (missing(cellnote)){
        cellnote = matrix("", ncol = ncol(x), nrow = nrow(x))
    }
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                     c("both", "row"))) {
            if (is.logical(Colv) && (Colv)){
                dendrogram = "column"
            }else{ 
                dedrogram = "none"
            }
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                    dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                     c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram = "row"
            else dendrogram = "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                    dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr = Rowv
        rowInd = order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr = hclustfun(distfun(x))
        ddr = as.dendrogram(hcr)
        ddr = reorder(ddr, Rowv)
        rowInd = order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv = rowMeans(x, na.rm = na.rm)
        hcr = hclustfun(distfun(x))
        ddr = as.dendrogram(hcr)
        ddr = reorder(ddr, Rowv)
        rowInd = order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd = nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc = Colv
        colInd = order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc = ddr
            colInd = order.dendrogram(ddc)
        }
        else colInd = rowInd
    }
    else if (is.integer(Colv)) {
        hcc = hclustfun(distfun(if (symm)
            x
            else t(x)))
        ddc = as.dendrogram(hcc)
        ddc = reorder(ddc, Colv)
        colInd = order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv = colMeans(x, na.rm = na.rm)
        hcc = hclustfun(distfun(if (symm)
            x
            else t(x)))
        ddc = as.dendrogram(hcc)
        ddc = reorder(ddc, Colv)
        colInd = order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd = 1:nc
    }
    retval$rowInd = rowInd
    retval$colInd = colInd
    retval$call = match.call()
    x = x[rowInd, colInd]
    x.unscaled = x
    cellnote = cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow = if (is.null(rownames(x)))
            (1:nr)[rowInd]
    else rownames(x)
    else labRow = labRow[rowInd]
    if (is.null(labCol))
        labCol = if (is.null(colnames(x)))
            (1:nc)[colInd]
    else colnames(x)
    else labCol = labCol[colInd]
    if (scale == "row") {
        retval$rowMeans = rm = rowMeans(x, na.rm = na.rm)
        x = sweep(x, 1, rm)
        retval$rowSDs = sx = apply(x, 1, sd, na.rm = na.rm)
        x = sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans = rm = colMeans(x, na.rm = na.rm)
        x = sweep(x, 2, rm)
        retval$colSDs = sx = apply(x, 2, sd, na.rm = na.rm)
        x = sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks = 16
        else breaks = length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks = seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                         length = breaks)
        else {
            extreme = max(abs(x), na.rm = TRUE)
            breaks = seq(-extreme, extreme, length = breaks)
        }
    }
    nbr = length(breaks)
    ncol = length(breaks) - 1
    if (class(col) == "function")
        col = col(ncol)
    min.breaks = min(breaks)
    max.breaks = max(breaks)
    x[x < min.breaks] = min.breaks
    x[x > max.breaks] = max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei = c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid = c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat = rbind(4:3, 2:1)
        
        if (!missing(ColSideColors)) {
            #if (!is.matrix(ColSideColors))
            #stop("'ColSideColors' must be a matrix")
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat = rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei = c(lhei[1], 0.2, lhei[2])
            lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
        }
        
        if (!missing(RowSideColors)) {
            #if (!is.matrix(RowSideColors))
            #stop("'RowSideColors' must be a matrix")
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat = cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid = c(lwid[1], 0.2, lwid[2])
            lwid = c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
        }
        lmat[is.na(lmat)] = 0
    }
    
    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op = par(no.readonly = TRUE)
    on.exit(par(op))
    
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
            par(mar = c(margins[1], 0, 0, 0.5))
            image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(rownames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), cex.axis=2, las = 2, tick = FALSE)
            }
        }
    }
    
    if (!missing(ColSideColors)) {
        
        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), cex.axis=1, las = 2, tick = FALSE)
            }
        }
    }
    
    par(mar = c(margins[1], 0, 0, margins[2]))
    x = t(x)
    cellnote = t(cellnote)
    if (revC) {
        iy = nr:1
        if (exists("ddr"))
            ddr = rev(ddr)
        x = x[, iy]
        cellnote = cellnote[, iy]
    }
    else iy = 1:nr
    if( showHeatmap ){
        image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    }else{
        image(1, 1, matrix(0, nrow=1, ncol=1), xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    }
    retval$carpet = x
    if (exists("ddr"))
        retval$rowDendrogram = ddr
    if (exists("ddc"))
        retval$colDendrogram = ddc
    retval$breaks = breaks
    retval$col = col
    if( showHeatmap ){
        if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
            mmat = ifelse(is.na(x), 1, NA)
            image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
                  col = na.color, add = TRUE)
        }
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
         cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
         cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale = min(breaks)
    max.scale = max(breaks)
    x.scaled = scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline = vline
        vline.vals = scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                       lty = 2)
            }
            xv = rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv = c(xv[1], xv)
            yv = 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline = hline
        hline.vals = scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv = rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv = rev(c(yv[1], yv))
            xv = length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    
    # plot row dendrogram
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
             col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    
    # plot col dendrogram
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }else{ 
        plot.new()
    }
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks = breaks
        if (symkey) {
            max.raw = max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw = -max.raw
            tmpbreaks[1] = -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] = max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw = min(x, na.rm = TRUE)
            max.raw = max(x, na.rm = TRUE)
        }
        
        z = seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
              xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv = pretty(breaks)
        xv = scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens = density(x, adjust = densadj, na.rm = TRUE)
            omit = dens$x < min(breaks) | dens$x > max(breaks)
            dens$x = dens$x[-omit]
            dens$y = dens$y[-omit]
            dens$x = scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                  lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h = hist(x, plot = FALSE, breaks = breaks)
            hx = scale01(breaks, min.raw, max.raw)
            hy = c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                  col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable = data.frame(low = retval$breaks[-length(retval$breaks)],
                                   high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}


df_hmr_recurrent_full = read.delim(fn_hmr_recurrent, sep='\t', header=T, check.names=F)
df_hmr_recurrent = df_hmr_recurrent_full

samples_plot = intersect(sample_ids_wgbs,colnames(df_hmr_recurrent))
mat_hmr_recurrent_full = as.matrix(df_hmr_recurrent_full[,samples_plot])
mat_hmr_recurrent_full[is.nan(mat_hmr_recurrent_full)] = NA

m = match.idx( dimnames(mat_hmr_recurrent_full)[[2]], samples_cmp )
col=rep("grey", 100)
col[m$idx.A] = "blue"
#m = match.idx( dimnames(mat_hmr_recurrent_full)[[2]], samples_tscnc )
#col[m$idx.A] = "gold"

m = match.idx( dimnames(mat_hmr_recurrent_full)[[2]],  dimnames(tumorpurity)[[1]] )
df_meth = data.frame( 
    median = colMedians( mat_hmr_recurrent_full, na.rm=TRUE ),
    mean = colMeans( mat_hmr_recurrent_full, na.rm=TRUE ),
    purity=tumorpurity$Tumor.Purity.Comp[m$idx.B],
    color=col,
    stringsAsFactors=FALSE)
df_meth = df_meth[order(df_meth$median),]
    
barplot(df_meth$median, col=df_meth$color)
boxplot( df_meth$median~df_meth$color, names=c("CMP", "non-CMP"), main="median methylation level", ylim=c(0, 60), las=1)


# restrict to the top 10% of variability in HMR by standard deviation
sds = apply(df_hmr_recurrent[,samples_plot],1,sd)
rows_top_decile = sds >= quantile( sds, 0.9, na.rm=T )
rows_top_decile[ is.na(rows_top_decile) ] = FALSE
df_hmr_recurrent = df_hmr_recurrent[rows_top_decile,]
mat_hmr_recurrent = as.matrix(df_hmr_recurrent[,samples_plot])

gr = makeGRangesFromDataFrame(df_hmr_recurrent)
cgi = rep('white',dim(mat_hmr_recurrent)[1])
rsidebar = t(matrix(cgi))


toplot = data.frame(samples_plot)
toplot$ERG = allele_effect('ERG')$alleles[samples_plot,'activating_sv']
toplot$ETV1 = allele_effect('ETV1')$alleles[samples_plot,'activating_sv']
toplot$ETV4 = allele_effect('ETV4')$alleles[samples_plot,'activating_sv']
toplot$ETV5 = allele_effect('ETV5')$alleles[samples_plot,'activating_sv']
toplot$CHD1 = allele_effect('CHD1')$alleles[samples_plot,'n_alleles_inactivated']==2
toplot$SPOP = allele_effect('SPOP')$alleles[samples_plot,'inactivating_missense']
toplot$RB1 = allele_effect('RB1')$alleles[samples_plot,'n_alleles_inactivated']==2
toplot$PTEN = allele_effect('PTEN')$alleles[samples_plot,'n_alleles_inactivated']==2
toplot$MYC = allele_effect('MYC')$alleles[samples_plot,'CNA_amp']
toplot$TP53 = allele_effect('TP53')$alleles[samples_plot,'n_alleles_inactivated']==2

CMP = samples_plot %in% samples_cmp
ETS = toplot$ERG|toplot$ETV1|toplot$ETV4|toplot$ETV5
CHD1_SPOP = toplot$CHD1 | toplot$SPOP
PTEN = toplot$PTEN
RB1 = toplot$RB1
TP53 = toplot$TP53
MYC = toplot$MYC
Mutation = samples_plot %in% c(samples_braf,samples_idh1,samples_tet2,samples_dnmt3b)
Sites = metastasis_locations[samples_plot]
print(fisher.test(CMP,ETS))
print(fisher.test(CMP,TP53))
print(fisher.test(CMP,MYC))
print(fisher.test(CMP,Mutation))
print(fisher.test(CMP,Sites))

colside = rep('white',dim(mat_hmr_recurrent)[2])
names(colside) = colnames(mat_hmr_recurrent)
stopifnot(length(colside)==length(samples_plot))

tSCNC = colside
tSCNC[samples_tscnc] = 'gray30'
stopifnot(length(tSCNC)==length(samples_plot))

ETS = colside
ETS[toplot$ERG] = 'gray30'
ETS[toplot$ETV1] = 'gray30'
ETS[toplot$ETV4] = 'gray30'
ETS[toplot$ETV5] = 'gray30'
stopifnot(length(ETS)==length(samples_plot))

MYC = colside
MYC[toplot$MYC] = 'gray30'
stopifnot(length(MYC)==length(samples_plot))

SPOP = colside
SPOP[toplot$SPOP] = 'gray30'
stopifnot(length(SPOP)==length(samples_plot))

RB1 = colside
RB1[toplot$RB1] = 'gray30'
stopifnot(length(RB1)==length(samples_plot))

PTEN = colside
PTEN[toplot$PTEN] = 'gray30'
stopifnot(length(PTEN)==length(samples_plot))

TP53 = colside
TP53[toplot$TP53] = 'gray30'
stopifnot(length(TP53)==length(samples_plot))

CMP = colside
CMP[samples_cmp] = 'white'#'lightsteelblue'
CMP[samples_braf] = 'orange'
CMP[samples_idh1] = 'green4'
CMP[samples_tet2] = 'purple'
CMP[samples_dnmt3b] = 'deeppink'
#CMP['DTB-252-BL'] = 'black'
stopifnot(length(CMP)==length(samples_plot))

Site = colside
Site[intersect(names(metastasis_locations)[metastasis_locations=='Bone'],samples_plot)] = wes_palette('Chevalier1')[1] #darkgreen
Site[intersect(names(metastasis_locations)[metastasis_locations=='Lymph_node'],samples_plot)] = wes_palette('Chevalier1')[2] #gold
Site[intersect(names(metastasis_locations)[metastasis_locations=='Liver'],samples_plot)] = wes_palette('Chevalier1')[4] #brown
Site[intersect(names(metastasis_locations)[metastasis_locations=='Other'],samples_plot)] = wes_palette('Chevalier1')[3] #gray

purity = colside
purity[tumorpurity[colnames(mat_hmr_recurrent),'Tumor.Purity.Histo']<50] = 'black'
stopifnot(length(purity)==length(samples_plot))

csidebar = data.matrix(cbind(Site,CMP,tSCNC,TP53,RB1,PTEN,MYC,ETS) )
color = colorRampPalette(c('blue','white','red'))(n = 1000)

hclust2 = function(x, method="average") {
    hclust(x, method=method)
}
dist2 = function(x) {
    as.dist(1-cor(t(x), method="spearman"))
}

dimnames(csidebar)[[2]][2] = "Select muts."
main_title = ""
fn_figure1b_noheat = '/notebook/human_sequence_prostate_WGBS/drafts/NatureGenetics/revision/submitted/figure_1b_noheat.pdf'
pdf( fn_figure1b_noheat, height=12, width=9)
h=heatmap.3(mat_hmr_recurrent, 
          showHeatmap=FALSE,
            hclustfun=hclust, 
          distfun=dist,
          na.rm=TRUE, 
          scale='none', 
          dendrogram="column", margins=c(8,1),
          Rowv=TRUE, Colv=TRUE, 
          ColSideColors=csidebar, 
          RowSideColors=rsidebar, 
          symbreaks=FALSE, 
          key=FALSE, 
          symkey=FALSE,
          lhei=c(1,5), 
          lwid=c(0.75,4), 
          cexCol=0.5,
          cexRow=1,
          density.info="none", trace="none",
          labCol=colnames(mat_hmr_recurrent), 
          labRow=rep("",(dim(mat_hmr_recurrent)[1])), 
          main=main_title, col=color,
          ColSideColorsSize=8, RowSideColorsSize=1)
dev.off()

