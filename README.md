TODO: Add more to the title of your project here

# sex_fiber_type_proteome_training:

TODO: Give a brief description of what your project is about

This project involves mass spectrometry-based proteomics of human skeletal muscle fibers in males and females
before and after 8 weeks of resistance training. Single muscle fibers were dissected from all participants and pooled
as type I and type II fibers at pre and post time points based on dot blotting against MYH7 and MYH2. This enables
comparison of differences in the skeletal muscle proteome between sexes and fiber types and how the sex- and
fiber type-specific proteome adapts to resistance training.

# Brief description of folder and file contents

TODO: As project evolves, add brief description of what is inside the data, doc and R folders.

The following folders contain:

- `data/`:Data aquired through handling of raw data such as normalization, log2 transformation and adding keywords
- `data-raw/`: Raw data and metadata about the samples 
- `R/`: All code generated in R
- `figures/`: All figures created for the manuscript

# Installing project R package dependencies

If dependencies have been managed by using `usethis::use_package("packagename")`
through the `DESCRIPTION` file, installing dependencies is as easy as opening the
`sex_fiber_type_proteome_training.Rproj` file and running this command in the console:

    # install.packages("remotes")
    remotes::install_deps()

You'll need to have remotes installed for this to work.

# Resource

For more information on this folder and file workflow and setup, check
out the [prodigenr](https://rostools.github.io/prodigenr) online
documentation.
