library("dplyr")
library("ggplot2")

## pad trial with frames after response
pad_trial <- function(x, trials) {
    trial_id <- x[["TrialID"]][1]
    clicked <- filter(trials, TrialID == trial_id) %>%
        `[[`("Clicked")
    sgn <- 1
    if (clicked == "nontarget") {
        sgn <- -1
    } else {}
    padded <- seq(max(x[["ms"]]), 5049, 20)
    res <- x
    if (length(padded)) {
        res <- bind_rows(x,
                         data_frame(TrialID = trial_id,
                                    ms = padded,
                                    dist_from_targ = sgn))
    } else {}
    res
}
## load in the data
session <- read.csv("data/session.csv", fileEncoding = "ascii")
trials <- read.csv("data/trials.csv", fileEncoding = "ascii")
list_info <- read.csv("data/list_info.csv", fileEncoding = "ascii")
test_info <- read.csv("data/test_info.csv", fileEncoding = "ascii")
mouse <- read.csv("data/mouse.csv", fileEncoding = "ascii")
screens <- read.csv("data/screens.csv", fileEncoding = "ascii")

## pull out mouse data from test trials only
## and time-lock to speech onset
mouse_test <- trials %>%
    inner_join(session, "SessionID") %>%
    inner_join(test_info, c("ListID", "ItemID")) %>%
    inner_join(mouse, c("TrialID")) %>%
    inner_join(screens %>% filter(identity == "target"),
               c("ItemID")) %>%
    mutate(ms = Msec - FPonset,
           raw_dist = (X - 512)/100 * (2 * (x1 == 612) - 1),
           dist_from_targ = ifelse(abs(raw_dist) > 1,
                                   1 * sign(raw_dist),
                                   raw_dist)) %>%
    filter(ms >= 0, Msec < ClickMs, ms <= 5049) %>%
    select(TrialID, ms, dist_from_targ) %>%
    as_data_frame()

mouse2 <- mouse_test %>%
    group_by(TrialID) %>%
    do(pad_trial(., trials))

## create means for each bin
mouse_bin <- mouse2 %>%
    inner_join(trials, "TrialID") %>%
    inner_join(session, "SessionID") %>%
    inner_join(test_info, c("ListID", "ItemID")) %>%
    mutate(bin = floor((ms + 50) / 100) * 100) %>%
    group_by(FillPause, Spkr, bin) %>%
    summarize(DFT = mean(dist_from_targ)) %>%
    ungroup() %>%
    mutate(FillPause = factor(FillPause), Spkr = factor(Spkr)) %>%
    filter(bin <= 5000)

## plot it
ggplot(mouse_bin, aes(bin, DFT, color = Spkr,
                      linetype = FillPause)) +
                          geom_point() + geom_line() +
                          ylab("Distance From Target")
