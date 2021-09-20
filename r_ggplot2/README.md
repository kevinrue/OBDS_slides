# Introduction to `ggplot2`

## Instructor(s)

- Kevin Rue-Albrecht (@kevinrue)
- Nicki Gray (@nickigray)

Past instructors:

- Lucy Garner (@lucygarner)

## Lesson goals and objectives

<!--
Refer to:
https://github.com/Bioconductor/BioC2019/blob/master/docs/workshop-syllabus.md#a-note-about-learning-goals-and-objectives-bloom
https://cft.vanderbilt.edu/guides-sub-pages/blooms-taxonomy/
-->

### Learning goals

<!--
High-level "big picture" objectives of the learning process.
-->

- Describe how to build and display `ggplot` figures.
- Identify the key differences between `ggplot` objects and base R graphics.
- Understand the concept of tidy data.

### Learning objectives

<!--
More concrete and measurable outputs.
-->

- Use the *[ggplot2](https://CRAN.R-project.org/package=ggplot2)* package to generate figures from various data sets.
- Overlay multiple layers for the same data on a single plot.
- Use faceting and plotting grids.
- Customise plot themes.

## Pre-requisites

- A clone of the shared GitHub repository for this course.
- A working installation of [R](https://www.r-project.org/) (4.0.3).
- A working installation of [git](https://git-scm.com/).
- A working installation of [RStudio](https://rstudio.com/).

## Data sets

- `iris` (built-in)
- `ChickWeight` (built-in)
- `ggplot2::diamonds` (built-in)
- `ggplot2::mtcars` (built-in)

## Time outline

| Activity                                      |  Time |
|-----------------------------------------------|-------|
| [Introduction to `renv`](../4_r_renv)         | 09:30 |
| **Morning Break**                             | 10:50 |
| Setup                                         | 11:00 |
| Lecture: Introduction to `ggplot2`            | 11:20 |
| Exercise: make a `ggplot` scatterplot         | 11:30 |
| Exercise: predict the difference              | 11:40 |
| Themes                                        | 11:50 |
| Plot grids                                    | 12:00 |
| Facets                                        | 12:10 |
| Save plot to file                             | 12:20 |
| **Lunch Break**                               | 12:30 |
| Exercise: smoothed trend lines                | 13:30 |
| Exercise: make a bar plot                     | 14:00 |
| Exercise: arrange a plot grid                 | 14:30 |
| **Afternoon Break**                           | 14:50 |
| Exercise: pair programming                    | 15:00 |
| **Day End**                                   | 16:00 |
