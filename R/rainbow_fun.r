rainbow_fun <- function(n, c=90, l=50, ...) {
if(requireNamespace("colorspace")) {
colorspace::rainbow_hcl(n, c = c, l = l, ...)		
	} else {
		rainbow(n, ...)
	}
}
