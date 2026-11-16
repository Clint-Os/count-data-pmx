# Modeling Count Data in Pharmacometrics

A three-part learning series on how to model and simulate count data, such as
seizures per week or crises per month, with every figure reproducible in R.

The running example is sickle cell anemia: vaso-occlusive crises counted per
patient per month, with hydroxyurea as the treatment. **All data are simulated.
No patient data are used.**

## Roadmap

| Part | What it covers | Post |
|---|---|---|
| 1 | Data types, the hazard, the Poisson model, mixed effects, the mean count function | [From events to the Poisson model](posts/01-from-events-to-the-poisson-model.md) |
| 2 | Over- and underdispersion, truncation, dependence between events, non-constant rates | [When the Poisson assumptions fail](posts/02-when-the-poisson-assumptions-fail.md) |
| 3 | Model-development workflow, diagnostics, simulation, drug-development example | Coming |

## Find your way around

| Folder | Contents |
|---|---|
| `posts/` | The posts |
| `R/` | One script per figure, plus the simulation scripts |
| `figures/` | The figures |
| `data/` | The simulated datasets |

## Run the code

Install the packages once, then run the scripts from the repository root:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork"))
```

```sh
Rscript R/01_simulate_trial.R
Rscript R/03_simulate_violations.R
for f in R/fig*.R; do Rscript "$f"; done
```

## Source

Based on Plan EL. Modeling and simulation of count data.
*CPT Pharmacometrics Syst Pharmacol.* 2014;3:e129. doi:10.1038/psp.2014.27

Related: [Markov-Models](https://github.com/Clint-Os/Markov-Models), on
continuous-time and hidden Markov models.
