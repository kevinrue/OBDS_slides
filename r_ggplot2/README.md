# Introduction to `ggplot2`

<!--
This title should match exactly the link in the main README.
-->

## Instructor(s)

<!--
Instructors should be listed in the order:
- Speaker
- Helper

Past instructors:

- Name
-->

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
This should be no more than 3 bullet points.
-->

- Describe how to build and display `ggplot` figures.
- Identify the key differences between `ggplot` objects and base R graphics.
- Understand the concept of tidy data.

### Learning objectives

<!--
More concrete and measurable outputs.
This can range from 3 to 8 bullet points.
-->

- Use the *[ggplot2](https://CRAN.R-project.org/package=ggplot2)* package to generate figures from various data sets.
- Overlay multiple layers for the same data on a single plot.
- Use faceting and plotting grids.
- Customise plot themes.

## Pre-requisites

<!--
May be a combination of:
- Requirements asked from participants before the day.
- Links to other OBDS course days with goals or objectives feeding in this day.
-->


- An account on the CCB cluster or [RStudio Cloud](https://rstudio.cloud/).
- Alternatively, a working installation of [R](https://www.r-project.org/) (4.1.2) and [RStudio Desktop](https://www.rstudio.com/products/rstudio/download/).
- A clone of the shared GitHub repository for this course.

## Data sets

<!--
Ideally, links to data sets that participants must download.
Even better, we add a page to this repository, that lists all data sets used; and this section links to some of those data sets.
Realistically, a list describing data sets that we will make them download on the day.
-->

- `iris` (built-in)
- `ChickWeight` (built-in)
- `ggplot2::diamonds` (built-in)
- `ggplot2::mtcars` (built-in)

## Time outline

<!--
Breakdown of time segments for lecture and exercises addressing the objectives listed above.
These are example times; adapt time, and insert/remove rows as needed.
Requirements:
- The day starts at 9:30
- There is a 10+ min break in the morning
- There is a 1+ h lunch break
- There is a 10+ min break in the afternoon
- The day ends at 16:00
-->

| Activity                                                      |  Time |
|---------------------------------------------------------------|-------|
| [Object oriented programming](../object_oriented_programming) | 09:30 |
| **Morning Break**                                             | 10:50 |
| Setup                                                         | 11:00 |
| Lecture: Introduction to `ggplot2`                            | 11:20 |
| Exercise: make a `ggplot` scatterplot                         | 11:40 |
| Exercise: predict the difference                              | 12:00 |
| Themes                                                        | 12:10 |
| Facets                                                        | 12:20 |
| **Lunch Break**                                               | 12:30 |
| Exercise: trend lines                                         | 13:30 |
| Exercise: bar plots                                           | 14:00 |
| Plot grids                                                    | 14:30 |
| Save plot to file                                             | 14:40 |
| **Afternoon Break**                                           | 14:50 |
| Exercise: arrange a plot grid                                 | 15:00 |
| Exercise: pair programming                                    | 15:20 |
| End of day: commit and push                                   | 15:40 |
| **Day End**                                                   | 16:00 |
