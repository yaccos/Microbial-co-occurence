# install.packages("renv")
# renv::activate()
# renv::restore()
# Uncomment the above lines to up the environment on R version 4.1.1
# If renv does not work, go through renv.lock and install the packages
# manually
source("read_selection_switch.R")
source("prepare_analysis.R")
source("run_HUNT_computations.R")
source("create_total_frame.R")
source("calculate_interaction_count.R")
source("export_interaction_count_ss.R")
source("generate_community_label_ss.R")
source("community_dynamics.R")
source("create_community_label_ss_robustness.R")
