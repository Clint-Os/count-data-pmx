# Modelling Count Data in Pharmacometrics

Companion code for a three-part blog series based on:

> Plan EL. Modeling and simulation of count data. *CPT Pharmacometrics Syst Pharmacol.* 2014;3:e129. doi:10.1038/psp.2014.27

The running example is a **simulated** two-arm trial of hydroxyurea in sickle cell
anaemia, with vaso-occlusive crises counted per patient per month. No patient data
are used.

## Layout

| Path | Contents |
|---|---|
| `posts/` | The posts, in Markdown with LaTeX equations |
| `R/00_setup.R` | Packages, colours, plot theme |
| `R/01_simulate_trial.R` | Simulates the trial and writes `data/voc_trial_sim.csv` |
| `R/fig0x_*.R` | One script per figure |
| `figures/` | The figures as PNG |
| `data/` | The simulated dataset |

## Reproducing the figures

Requires R (4.1 or later, for the `|>` pipe) and four packages:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork"))
```

Run the scripts from the repository root:

```sh
Rscript R/01_simulate_trial.R
for f in R/fig0*.R; do Rscript "$f"; done
```

## Converting a post

The posts are plain Markdown. To paste one into a platform that does not render
LaTeX, convert it with [pandoc](https://pandoc.org):

```sh
pandoc posts/01-from-events-to-the-poisson-model.md -o post01.docx                           # Word, native equations
pandoc posts/01-from-events-to-the-poisson-model.md -s --mathml --embed-resources -o post01.html   # web page
```

Run these from inside `posts/` so the figure paths resolve.
