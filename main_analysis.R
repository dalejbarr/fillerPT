## analysis done as ANOVA rather than as t-test as in original paper 
library("dplyr")
library("tidyr")

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

## distance travelled during filled interval
dist_fill_int <- mouse2 %>%
    mutate(ms2 = ms - 1871) %>%
    group_by(TrialID) %>%
    filter(ms == min(ms) | abs(ms2) == min(abs(ms2))) %>%
    slice(1:2) %>%
    mutate(win = paste0("W", row_number())) %>%
    ungroup() %>%
    select(-ms, -ms2) %>%
    spread(win, dist_from_targ) %>%
    mutate(win = "filled_int", diff = W2 - W1)

## distance travelled during referring expression
dist_ref <- mouse2 %>%
    mutate(ms2 = ms - 1871,
           ms3 = ms - 4471) %>%
    group_by(TrialID) %>%
    filter(abs(ms2) == min(abs(ms2)) |
           abs(ms3) == min(abs(ms3))) %>%
    group_by(TrialID) %>%    
    slice(1:2) %>%
    mutate(win = paste0("W", row_number())) %>%
    ungroup() %>%
    select(-ms, -ms2, -ms3) %>%    
    spread(win, dist_from_targ) %>%
    mutate(win = "ref", diff = W2 - W1)

dist <- bind_rows(dist_fill_int, dist_ref)

## link to trial information
dist2 <- dist %>%
    inner_join(trials, "TrialID") %>%
    inner_join(session, "SessionID") %>%
    inner_join(test_info, c("ListID", "ItemID")) %>%
    select(SessionID, ItemID, win, diff, FillPause, Spkr)

## recreate table 1
dist2 %>%
    group_by(win, FillPause, Spkr) %>%
    summarize(mDist = mean(diff, na.rm = TRUE))

## analysis
dist_aov <- dist2 %>%
    mutate(Spkr = factor(Spkr),
           SessionID = factor(SessionID),
           FillPause = factor(FillPause)) %>%
    group_by(win, SessionID, Spkr, FillPause) %>%
    summarize(mDist = mean(diff)) %>%
    ungroup()

aov_filled <- aov(mDist ~ Spkr * FillPause + Error(SessionID),
                  data = subset(dist_aov, win == "filled_int"))

summary(aov_filled)

aov_ref <- aov(mDist ~ Spkr * FillPause + Error(SessionID),
               data = subset(dist_aov, win == "ref"))

summary(aov_ref)
