# Matrix exports the shared S4 generics for these base operations. GPUmatrix
# registers its methods when its namespace loads, without needing to export
# the generic names itself. nrow/ncol use base's dim() dispatch.
#' @importFrom GPUmatrix gpu.matrix
#' @importFrom Matrix as.matrix colSums mean rowSums t
NULL
