# =====================================================================
# Capture the R environment that runs a set of pipeline scripts.
#   Rscript env_capture.R <out_dir> <label> <script1.R> [script2.R ...]
# Writes to <out_dir>:
#   system.txt               OS, CPU, R, BLAS/LAPACK, RNG, modules
#   sessionInfo.txt          sessionInfo() with every script package attached
#   packages_used.csv        packages the scripts call + all their dependencies
#   packages_installed.csv   every installed package, version, library, build
#   renv.lock                renv lockfile for packages_used (if renv available)
# Read-only with respect to the scripts.
# =====================================================================
a <- commandArgs(trailingOnly = TRUE)
out <- a[1]; label <- a[2]; scripts <- a[-(1:2)]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# --- packages referenced by the scripts ------------------------------
src <- unlist(lapply(scripts, readLines, warn = FALSE))
src <- sub("#.*$", "", src)                                   # drop comments
hits <- c(regmatches(src, gregexpr("library\\(([A-Za-z][A-Za-z0-9.]*)\\)", src)),
          regmatches(src, gregexpr("requireNamespace\\(\"([A-Za-z][A-Za-z0-9.]*)\"", src)),
          regmatches(src, gregexpr("\\b([A-Za-z][A-Za-z0-9.]*)::", src)))
direct <- unique(gsub("library\\(|\\)|requireNamespace\\(|\"|::", "", unlist(hits)))
base_pk <- rownames(installed.packages(priority = "base"))
direct <- sort(setdiff(direct, base_pk))

# --- system ----------------------------------------------------------
sink(file.path(out, "system.txt"))
cat("Label:   ", label, "\nCaptured:", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), "\n")
si <- Sys.info(); cat("Host:    ", si[["nodename"]], "| OS:", si[["sysname"]], si[["release"]], "| machine:", si[["machine"]], "\n")
if (file.exists("/etc/os-release")) cat(grep("PRETTY_NAME", readLines("/etc/os-release"), value = TRUE), "\n")
if (nzchar(Sys.which("sw_vers"))) cat(system("sw_vers", intern = TRUE), sep = "\n")
cpu <- if (file.exists("/proc/cpuinfo")) unique(grep("model name", readLines("/proc/cpuinfo"), value = TRUE)) else
  tryCatch(system("sysctl -n machdep.cpu.brand_string", intern = TRUE), error = function(e) "?")
cat("CPU:     ", cpu, "\n")
cat("R:       ", R.version.string, "\nR_HOME:  ", R.home(), "\n")
cat("BLAS:    ", extSoftVersion()[["BLAS"]], "\nLAPACK:  ", La_library(), "\n")
for (v in c("FLEXIBLAS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "LOADEDMODULES"))
  if (nzchar(Sys.getenv(v))) cat(sprintf("%-8s: %s\n", v, Sys.getenv(v)))
cat("RNGkind: ", paste(RNGkind(), collapse = " / "), "\n")
cat("Library paths:\n"); cat(paste0("  ", .libPaths()), sep = "\n")
cat("Scripts scanned:\n"); cat(paste0("  ", basename(scripts)), sep = "\n")
cat("Packages called directly:", paste(direct, collapse = ", "), "\n")
missing <- direct[!vapply(direct, requireNamespace, logical(1), quietly = TRUE)]
cat("NOT INSTALLED here:", if (length(missing)) paste(missing, collapse = ", ") else "none", "\n")
sink()

# --- sessionInfo with everything attached ------------------------------
for (p in setdiff(direct, missing))
  suppressPackageStartupMessages(try(library(p, character.only = TRUE), silent = TRUE))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))

# --- installed + used packages ---------------------------------------
ip <- installed.packages(fields = c("Built", "Repository", "RemoteType", "RemoteSha", "Packaged"))
ip <- as.data.frame(ip[, intersect(c("Package", "Version", "LibPath", "Priority", "Built",
                                     "Repository", "RemoteType", "RemoteSha"), colnames(ip))],
                    stringsAsFactors = FALSE)
ip <- ip[!duplicated(ip$Package), ]                           # first on .libPaths() wins
write.csv(ip[order(tolower(ip$Package)), ], file.path(out, "packages_installed.csv"), row.names = FALSE)

deps <- tools::package_dependencies(setdiff(direct, missing), db = installed.packages(),
                                    which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE)
used <- sort(unique(c(setdiff(direct, missing), unlist(deps))))
used <- setdiff(used, base_pk)
u <- ip[ip$Package %in% used, ]
u$called_directly <- u$Package %in% direct
write.csv(u[order(!u$called_directly, tolower(u$Package)), ], file.path(out, "packages_used.csv"), row.names = FALSE)

# --- renv lockfile for exactly the used packages ----------------------
if (requireNamespace("renv", quietly = TRUE)) {
  proj <- file.path(tempdir(), "envproj"); dir.create(proj, showWarnings = FALSE)
  ok <- tryCatch({
    renv::snapshot(project = proj, lockfile = file.path(out, "renv.lock"),
                   packages = used, prompt = FALSE, force = TRUE)
    TRUE
  }, error = function(e) { message("renv snapshot failed: ", conditionMessage(e)); FALSE })
  cat("renv.lock written:", ok, "\n")
}
cat(sprintf("%s: %d packages called directly, %d incl. dependencies, %d installed in total\n",
            label, length(direct), nrow(u), nrow(ip)))
