#' @export
setup_aw <- function(run = TRUE){
  if(run == TRUE){
    # install_github("trinker/pacman")
    library(pacman) # package install/load
    pacman::p_load(
      ### Project and File Management ###
      here,          # file paths relative to R project root folder
      rio,           # import/export of many types of data; https://cran.r-project.org/web/packages/rio/vignettes/rio.html
      openxlsx,      # import/export of multi-sheet Excel workbooks
      remotes,

      ### Data Wrangling ###
      devtools,      # package development tools
      readxl,        # import excel files
      readr,         # friendly way to read rectangular data
      data.table,    # working with tabular data
      tidyverse,     # tidy data wrangling and presentation
      #dplyr,      # data management
      #tidyr,      # data management
      #ggplot2,    # data visualization
      #stringr,    # work with strings and characters
      #forcats,    # work with factors
      #lubridate,  # work with dates; Sys.Date()
      #purrr       # iteration and working with lists
      magrittr,      # pipes
      naniar,        # assessing missing data
      mice,          # missing data imputation
      zoo,           # ordered observations
      glue,          # format and interpolate a string
      dataCompareR,  # match and confirm differences between dataframes
      DataEditR,     # interactively view, enter, filter and edit data

      ### Statistics ###
      doBy,          # working with grouped data
      vtable,        # outputting automatic variable documentation
      janitor,       # tables and data cleaning
      rstatix,       # quickly run statistical tests and summaries
      broom,         # tidy up results from regressions
      lmtest,        # likelihood-ratio tests
      emmeans,       # estimated marginal means (EMMs) for many linear, generalized linear, and mixed models
      psych,         # multivariate analysis for psychometric theory
      multcomp,      # simultaneous tests and confidence intervals for general linear hypotheses
      multcompView,  # convert a logical vector or a vector of p-values or a correlation, difference, or distance matrix into a display
      mcr,           # Method Comparison Regression
      effectsize,    # indices of effect size and standardized parameters
      stats,
      #statcheck,     # automatically extract statistical null-hypothesis significant testing (NHST) results from articles
      pastecs,      # regularisation, decomposition and analysis of space-time series
      outliers,      #
      frequency,     #
      diffobj,       #
      mice,          #
      wrangle,       #
      sampling,      #
      icarus,        #
      survey,        #
      simPop,        #
      iterpc,        # can be used to obtain a listing of all possible samples (sample space) for w/ and w/o replacement sampling

      ### Data Visualization ###
      DataExplorer,  # first look of the data
      knitr,         # R Markdown report generation and html tables
      kableExtra,    # build common complex tables and manipulate table styles
      flextable,     # HTML tables
      scales,        # graphical scales map data to aesthetics
      lattice,       # visualize multivariate data
      parzer,        # parse geographic coordinates
      plotly,        # interactive graphs
      finalfit,      # elegant final results tables and plots
      visreg,        # visualize the fit of regression models
      # ggtree,        # phylogenetics

      ### Soup Up ggplot ###
      ggforce,       # collection of stats and geoms for ggplot2
      showtext,
      RColorBrewer,  # display color palettes
      ggfortify,     # define fortify and autoplot
      gridExtra,     # arrange grid plots
      ggh4x,         # change facet sizes
      plotly,        # interactive, publication-quality graphs
      ggthemes,      # extra themes, scales, and geoms
      ggpubr,        # publication quality figures
      ggsignif,      # add group-wise comparisons
      cowplot,       # publication quality figures, aligning graphs to grids
      ggbreak,       # scale functions for setting axis breaks
      patchwork,     # add ggplots together for multiplot layouts
      paletteer,     # streamline color palettes across R
      viridis,       # viridis color scale
      hrbrthemes,    # typography-centric themes and theme components
      esquisse,      # interactive drag and drop make ggplots
      #rayshader,     # 3D topography plots
      emojifont,     # using emoji font in R
      gganimate,     # include the description of animation
      ggstatsplot,   # creating graphics with details from statistical tests included
      EnvStats,      # graphical and statistical analyses of environmental data (stat_n_text())
      ggrepel,       # labeling points
      tidytext,      # reorder within
      ggalt,         # encircle, dumpbell, xspline. etc.
      wesanderson,   # palettes
      jcolors,       # palettes
      ggsci,         # palettes

      ### Modeling ###
      tidymodels,    # meta-package for modeling and machine learning
      mgcv,          # Mixed GAM Computation Vehicle
      modelr,        # elegant pipelines when modelling
      # lme4,          # fitting and analyzing mixed models
      gamm4,         # based on gamm from package mgcv, but uses lme4 rather than nlme
      leaps,         # all subsets regression
      car,           # fitted regression model
      lmerTest,      # provides p-values for anova and summary tables for linear mixed models
      drc,           # Dose Response Curve
      ape,           # Analyses of Phylogenetics and Evolution
      survival,      # core survival analysis routines
      survminer,     # facilitating survival analysis and visualization
      caret,         # Classification And REgression Training
      rpart,         # Recursive Partitioning and Regression Trees
      randomForest,  # Breiman and Cutler's Random Forests for Classification and Regression
      mlr3,          # object-oriented programming on the building blocks of machine learning
      xgboost,       # efficient linear model solver and tree learning algorithms
      performance,   # evaluate the quality of modelfit
      effects,       # interactions for various statistical models with linear predictors
      AICcmodavg,    # AIC comparison

      ### Miscellaneous ###
      randomizr,     # helpful for distributing treatments
      checkpoint,    # rewinds time for packages, good for reproducibility
      flow,          # browse code in flow diagrams
      weathermetrics,# convert temperatures
      tinytex,       # custom LaTeX distribution

      ### Mapping ###
      rnaturalearthdata,
      #rnaturalearthhires,
      rnaturalearth, # vector data commonly used in world mapping
      sf,            # encode spatial vector data

      ### Ross Cunning ###
      steponeR,      # devtools::install_github("jrcunning/steponeR"); QuantStudio StepOneR
      IPAM2R         # devtools::install_github("jrcunning/IPAM2R"); WALZ Imaging PAM
    )

    showtext_auto()
    showtext_opts(dpi = 300)

    theme_aw <- function(){
      theme_bw(base_size = 14,
               base_family = 'serif') +
        theme(legend.position = 'bottom')
    }

    ggsave_aw <- function(filename = default_name(plot), height= 8, width= 10, ...) {
      ggsave(filename=filename, height=height, width=width, ...)
    }
  }
}
