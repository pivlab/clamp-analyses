CLAMP_REPO <- "chikinalab/CLAMP"
CLAMP_REF <- "818e13ba55d66840e0710c3f1ac15f6d97e1dd8b"

args <- commandArgs(trailingOnly = TRUE)
unknown_args <- setdiff(args, "--check")
if (length(unknown_args) > 0) {
    stop("Unknown argument(s): ", paste(unknown_args, collapse = ", "))
}

installed_sha <- function() {
    if (!requireNamespace("CLAMP", quietly = TRUE)) {
        return(NA_character_)
    }
    description <- utils::packageDescription("CLAMP")
    sha <- description$RemoteSha
    if (is.null(sha) || !nzchar(sha)) NA_character_ else tolower(sha)
}

is_pinned <- function() identical(installed_sha(), tolower(CLAMP_REF))

if ("--check" %in% args) {
    sha <- installed_sha()
    if (is_pinned()) {
        cat("CLAMP is installed at the required revision", CLAMP_REF, "\n")
        quit(status = 0)
    }
    cat(
        "CLAMP revision mismatch: expected", CLAMP_REF,
        "but found", if (is.na(sha)) "an unverified/local installation" else sha,
        "\n", file = stderr()
    )
    quit(status = 1)
}

message(
    "WARNING: clamp-analyses temporarily requires CLAMP revision ", CLAMP_REF,
    ". Do not update CLAMP independently; a follow-up PR will migrate these analyses ",
    "to the latest compatible CLAMP release."
)

if (is_pinned()) {
    message("CLAMP is already installed at the required revision; skipping installation.")
    quit(status = 0)
}

if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
}

remotes::install_github(
    repo = CLAMP_REPO,
    ref = CLAMP_REF,
    upgrade = "never",
    dependencies = FALSE,
    force = TRUE
)

if (!is_pinned()) {
    stop("CLAMP installation completed, but the installed RemoteSha does not match ", CLAMP_REF)
}

message("Installed CLAMP revision ", CLAMP_REF, ".")
