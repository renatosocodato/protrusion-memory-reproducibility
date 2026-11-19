versions <- c(ggplot2="4.0.2", cowplot="1.2.0", jsonlite="2.0.0",
              ragg="1.5.0", svglite="2.2.2", scales="1.4.0", systemfonts="1.3.1")
if (!requireNamespace("remotes", quietly=TRUE))
  install.packages("remotes", repos="https://cloud.r-project.org")
for (package in names(versions)) {
  installed <- requireNamespace(package, quietly=TRUE) &&
    as.character(utils::packageVersion(package)) == versions[[package]]
  if (!installed)
    remotes::install_version(package, versions[[package]], repos="https://cloud.r-project.org",
                             upgrade="never", dependencies=NA)
}
