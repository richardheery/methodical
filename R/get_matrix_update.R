get_matrix_updated = function (m, type = "M", add_loci = FALSE, in_granges = FALSE) 
{
    if (!is(m, "methrix")) {
        stop("A valid methrix object needs to be supplied.")
    }
    type <- match.arg(arg = type, choices = c("M", "C"))
    if (add_loci == FALSE & in_granges == TRUE) {
        warning("Without genomic locations (add_loci= FALSE), it is not possible to convert the results to GRanges, ", 
            "the output will be a data.table object. ")
    }
    if (type == "M") {
        d <- SummarizedExperiment::assay(x = m, i = 1)
    }
    else {
        d <- SummarizedExperiment::assay(x = m, i = 2)
    }
    if (add_loci) {
        if (methrix:::is_h5(m)) {
            d <- as.data.frame(cbind(SummarizedExperiment::rowData(x = m), 
              # Added as.matrix here
                as.data.frame(as.matrix(d))))
        }
        else {
            d <- as.data.frame(cbind(SummarizedExperiment::rowData(x = m), 
                d))
        }
        if (in_granges) {
            d$end <- d$start + 1
            d <- GenomicRanges::makeGRangesFromDataFrame(d, keep.extra.columns = TRUE)
        }
        else {
            data.table::setDT(x = d)
        }
    }
    d
}
