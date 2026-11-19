args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2, args[[2]] %in% c("calibration", "heldout"))
source <- new.env(parent=emptyenv())
stopifnot(identical(load(args[[1]], envir=source), "d"))
d <- source$d
required <- c("plate_ID", "well_ID", "v", "treatment", "dose", "bad_well", "experiment")
stopifnot(all(required %in% names(d)), !anyNA(d[, required]))
eligible <- d$bad_well == 0
treated <- d$treatment == "CK666" & d$dose == "40"
control <- d$treatment == "DMSO" & d$dose == "DMSO"
d <- d[eligible & (treated | control), required]
stopifnot(is.numeric(d$v), all(is.finite(d$v)), all(d$v >= 0))
means <- aggregate(v ~ plate_ID + well_ID + treatment + experiment, d, mean)
stopifnot(all(means$v > 0))
plates <- sort(unique(as.character(means$plate_ID)))
rows <- lapply(plates, function(plate_id) {
    plate <- means[means$plate_ID == plate_id, ]
    vehicle <- plate$v[plate$treatment == "DMSO"]
    ck666 <- plate$v[plate$treatment == "CK666"]
    stopifnot(length(vehicle) > 0, length(ck666) > 0)
    if (args[[2]] == "calibration") {
        labels <- c("Arp2/3" = "d1-arp", "Statins" = "d1-statins")
        preparation <- unname(labels[as.character(unique(plate$experiment))])
        stopifnot(length(preparation) == 1, !is.na(preparation))
    } else {
        stopifnot(all(plate$experiment == "Pimozide"))
        preparation <- paste0("d2-", sub("\\.[12]$", "", plate_id))
    }
    data.frame(preparation=preparation, unit=plate_id, scenario="ck666_40",
               channel="log_speed_ratio", time_s=0,
               value=mean(log(ck666)) - mean(log(vehicle)), measurement_unit="log_ratio")
})
write.csv(do.call(rbind, rows), stdout(), row.names=FALSE)
