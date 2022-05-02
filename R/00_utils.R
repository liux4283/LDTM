.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to my package")
}


# .onLoad <- function(libname, pkgname) {
#   op <- options()
#   op.devtools <- list(
#     devtools.path = "~/R-dev",
#     devtools.install.args = "",
#     devtools.name = "Randy",
#     devtools.desc.author = "Randy Jin <xin.2.jin@cuanschutz.edu> [aut, cre]",
#     devtools.desc.license = "What license is it under?",
#     devtools.desc.suggests = NULL,
#     devtools.desc = list()
#   )
#   toset <- !(names(op.devtools) %in% names(op))
#   if(any(toset)) options(op.devtools[toset])
#
#   invisible()
# }


