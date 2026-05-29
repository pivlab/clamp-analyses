if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
}

remotes::install_github(
    repo = "chikinalab/CLAMP",
    ref = "818e13ba55d66840e0710c3f1ac15f6d97e1dd8b",
    upgrade = "never",
    dependencies = FALSE
)
