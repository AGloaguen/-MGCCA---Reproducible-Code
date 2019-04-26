#' The function derivation() allows computing the derivate of the g function. The g function must be convex
#' @param g  a convex function 
#' @return \item{dg}{the derivate of the g function}
#' @title compute the derivate of the g function
#' @export derivation

derivation = function(g){
arg <- formalArgs(g)
gchar <- gsub(sprintf("function\\(%s\\)", arg),"", capture.output(g))
dg <- deriv(as.formula( paste("~", gchar) ), arg) 
}