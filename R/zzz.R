# zzz.R is loaded last (alphabetically) among all R files in the package,
# ensuring that all package code is in place before this runs.

# Problem: ChIPseeker registers a broken mutate.GRanges S3 method that
# overwrites plyranges' correct version in the S3 dispatch table. Load order
# determines which version wins, and ChIPseeker typically loads after plyranges,
# so its broken version takes precedence. The result is that dplyr::mutate()
# fails on GRanges objects when both packages are loaded.

# Why users can't work around it themselves: the problem is in S3 dispatch,
# not in which generic is called. Even dplyr::mutate(gr, ...) still dispatches
# to mutate.GRanges via S3, and finds ChIPseeker's version.

# Fix: when profileplyr attaches (i.e. when the user runs library(profileplyr)),
# re-register plyranges' correct mutate.Ranges method against dplyr's generic.
# This puts plyranges' method back on top in the dispatch table, regardless of
# what order packages were loaded in.

# we've made plyranges an Imported package, so mutate will always use the plyranges version
# if it was in suggests, then the fact its wrapped in requireNamespace() means it only fires if the user has
# plyranges installed — profileplyr would then not gain a hard dependency on plyranges.

.onAttach <- function(libname, pkgname) {
  if (requireNamespace("plyranges", quietly = TRUE)) {
    plyranges_mutate <- getFromNamespace("mutate.Ranges", "plyranges")
    registerS3method("mutate", "GRanges", plyranges_mutate,
                     envir = asNamespace("dplyr"))
    registerS3method("mutate", "Ranges", plyranges_mutate,
                     envir = asNamespace("dplyr"))

    plyranges_filter <- getFromNamespace("filter.Ranges", "plyranges")
    registerS3method("filter", "GRanges", plyranges_filter,
                     envir = asNamespace("dplyr"))
    registerS3method("filter", "Ranges", plyranges_filter,
                     envir = asNamespace("dplyr"))

    plyranges_arrange <- getFromNamespace("arrange.Ranges", "plyranges")
    registerS3method("arrange", "GRanges", plyranges_arrange,
                     envir = asNamespace("dplyr"))
    registerS3method("arrange", "Ranges", plyranges_arrange,
                     envir = asNamespace("dplyr"))
  }
}
