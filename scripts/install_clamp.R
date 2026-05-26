if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
}

remotes::install_github(
    repo = "chikinalab/CLAMP",
    ref = "4a6a32006624b942c847becd71f73baf7dedfed6",
    upgrade = "never",
    dependencies = FALSE
)
