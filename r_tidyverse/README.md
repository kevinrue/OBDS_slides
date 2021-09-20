# The `tidyverse`

<!--
This title should match exactly the link in the main README.
-->

## Instructor(s)

<!--
Instructors should be listed in the order:
- Speaker
- Helper
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

- Describe tidy data.
- Appreciate that scripts can combine base R and tidyverse.
- Assemble iteratively longer workflows using the tidyverse pipe.

### Learning objectives

<!--
More concrete and measurable outputs.
This can range from 3 to 8 bullet points.
-->

- Import and manipulate tidy data using packages from the tidyverse.
- Compute summary statistics on groups of observations in data frames.
- Filter data frames to select observations that meet certain criteria.
- Combine information from multiple data frames using matching information in certain columns.

## Pre-requisites

<!--
May be a combination of:
- Requirements asked from participants before the day.
- Links to other OBDS course days with goals or objectives feeding in this day.
-->

- A clone of the shared GitHub repository for this course.
- A working installation of [R](https://www.r-project.org/) (4.1).
- A working installation of [git](https://git-scm.com/).
- A working installation of [RStudio](https://rstudio.com/).

## Data sets

<!--
Ideally, links to data sets that participants must download.
Even better, we add a page to this repository, that lists all data sets used; and this section links to some of those data sets.
Realistically, a list describing data sets that we will make them download on the day.
-->

- [`iris.csv`](https://github.com/OBDS-Training/OBDS_Syllabus/blob/main/datasets/iris.csv)
- [EH2011.xlsx](https://github.com/OBDS-Training/OBDS_Syllabus/blob/main/datasets/EH2011.xlsx)

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

| Activity                                        |  Time |
|-------------------------------------------------|-------|
| Setup                                           |  9:30 |
| Introduction to tidy data                       | 10:00 |
| Import and export data                          | 10:15 |
| **Morning Break**                               | 10:50 |
| The `%>%` operator and the `tibble` class       | 11:00 |
| Manipulate tidy data using `dplyr`              | 11:15 |
| **Lunch Break**                                 | 12:30 |
| Tidy data using `tidyr`                         | 13:30 |
| Work with strings using `stringr`               | 14:00 |
| Integrated exercise                             | 14:15 |
| **Afternoon Break**                             | 14:50 |
| Integrated exercise                             | 15:00 |
| **Day End**                                     | 16:00 |
