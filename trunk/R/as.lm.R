
# if object is an "nls" object then its often used like this:
# predict(as.lm(object), ...) where ... are any predict.lm args

# as.lm.nls effectively just does this:
# lm(lhs ~ gradient - 1, offset = fitted(object),
#   list(gradient = object$m$gradient(), lhs = object$m$lhs()))
# so most of the code is just to get the names right.

as.lm <- function(object, ...) UseMethod("as.lm")

as.lm.nls <- function(object, ...) {
    if (!inherits(object, "nls")) {
		w <- paste("expected object of class nls but got object of class:", 
			paste(class(object), collapse = " "))
		warning(w)
	}

	gradient <- object$m$gradient()
	if (is.null(colnames(gradient))) {
		colnames(gradient) <- names(object$m$getPars())
	}

	response.name <- if (length(formula(object)) == 2) "0" else 
		as.character(formula(object)[[2]])

	lhs <- object$m$lhs()
	L <- data.frame(lhs, gradient)
	names(L)[1] <- response.name

	fo <- sprintf("%s ~ %s - 1", response.name, 
		paste(colnames(gradient), collapse = "+"))
	fo <- as.formula(fo, env = as.proto.list(L))

	do.call("lm", list(fo, offset = substitute(fitted(object))))

}

