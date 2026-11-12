# Figure 2 -- the same events seen as exact times (RTTE) and as monthly counts
source("R/00_setup.R")
set.seed(22)   # chosen so that one month holds two events

h    <- 4.5 / 12                           # constant hazard, events per month
gaps <- rexp(40, rate = h)                 # waiting times between events
t_ev <- cumsum(gaps)
t_ev <- t_ev[t_ev <= 12]                   # keep the first 12 months

events <- tibble(time = t_ev)
counts <- tibble(month = 1:12,
                 count = tabulate(ceiling(t_ev), nbins = 12))

p_events <- ggplot(events, aes(time, 0)) +
  geom_vline(xintercept = 0:12, colour = "#e9e8e4", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = col_muted, linewidth = 0.3) +
  geom_point(shape = 124, size = 7, colour = col_hu) +
  scale_x_continuous(breaks = 0:12, limits = c(0, 12), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(title = "What happened: exact event times (repeated time-to-event data)",
       x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank())

p_counts <- ggplot(counts, aes(month - 0.5, count)) +
  geom_col(fill = col_placebo, width = 0.9) +
  geom_text(aes(label = count), vjust = -0.5, size = 3.4, colour = col_ink) +
  scale_x_continuous(breaks = 0:12, limits = c(0, 12), expand = c(0.01, 0)) +
  scale_y_continuous(breaks = 0:3, limits = c(0, max(counts$count) + 0.6)) +
  labs(title = "What we recorded: events per month (count data)",
       x = "Time (months)", y = "VOCs per month") +
  theme(panel.grid.major.x = element_blank())

save_fig(p_events / p_counts + plot_layout(heights = c(1, 2.2)),
         "fig02_events_to_counts.png", height = 4.8)
