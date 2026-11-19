versions <- c(ggplot2="4.0.2", cowplot="1.2.0", jsonlite="2.0.0",
              ragg="1.5.0", svglite="2.2.2", scales="1.4.0", systemfonts="1.3.1")

has_version <- function(package) {
  path <- find.package(package, quiet=TRUE)
  length(path) == 1L &&
    identical(utils::packageDescription(package, lib.loc=dirname(path), fields="Version"),
              versions[[package]])
}

available <- vapply(names(versions), has_version, logical(1))
if (!all(available)) {
  if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes", repos="https://cloud.r-project.org")
  for (package in names(versions)[!available]) {
    remotes::install_version(package, versions[[package]], repos="https://cloud.r-project.org",
                             upgrade="never", dependencies=NA)
    if (!has_version(package))
      stop(sprintf("Could not install %s %s", package, versions[[package]]))
  }
}
