# Follow installation instructions here https://github.com/rmcelreath/rethinking

library(tidyverse)
library(rstan)

mydarktheme <-
    theme_minimal() +
    theme(
        panel.background = element_rect(fill = rgb(.2, .21, .27)),
        plot.background = element_rect(fill = rgb(.2, .21, .27)),
        text = element_text(colour = "grey"),
        axis.text = element_text(colour = "grey"),
        panel.grid = element_line(colour = "grey"),
        strip.text = element_text(color = "grey")
    )

d <- read_csv("data/DoseResponseMaino2018.csv", col_types = "ccniii") %>%
    mutate(dead = dead + incapacitated) %>%
    dplyr::select(-incapacitated) %>%
    group_by(treatment, pop, dose) %>%
    summarise_all(sum) %>%
    mutate(total = dead + alive) %>%
    filter(dose != 0)
d

m1 <- glm(cbind(dead, alive) ~ pop / log(dose) - 1,
    family = binomial("logit"),
    data = d
)

model.matrix(m1)

summary(m1)

m1pred <- expand.grid(
    pop = c("control", "resistant"),
    dose = 10^seq(-2, 4, length = 100)
)
m1pred$logdose <- log(m1pred$dose)
pred <- predict(m1, newdata = m1pred, type = "response", se.fit = TRUE)
m1pred$p <- pred$fit
m1pred$se <- pred$se.fit

ggplot() +
    geom_point(data = d, aes(dose, dead / total, colour = pop)) +
    geom_line(data = m1pred, aes(dose, p, colour = pop)) +
    geom_ribbon(
        data = m1pred,
        aes(dose, ymin = p - 1.96 * se, max = p + 1.96 * se, fill = pop),
        alpha = 0.5
    ) +
    coord_cartesian(expand=FALSE, ylim=c(0,1), xlim=c(1e-2, 1e4)) +
    mydarktheme +
    xlab("omethoate (mg/L)") +
    ggtitle("Frequentist") +
    scale_x_log10()
ggsave("plots/plot1.png", height = 6, width = 10)

db <- d %>%
    rowwise() %>%
    mutate(outcome = list(rep(c(0, 1), c(alive, dead)))) %>%
    dplyr::select(-alive, -dead, -total) %>%
    unnest(cols = outcome)
head(db)

m1 <- glm(outcome ~ pop / log(dose) - 1,
    family = binomial("logit"),
    data = db
)
summary(m1)

############ FIT BAYESIAN MODEL WITH STAN ###############

stan_data <- list(
    N = nrow(db),
    y = db$outcome,
    pop = ifelse(db$pop == "control", 1, 2),
    dose = db$dose
)

fit <- stan(
    file = "model.stan",
    data = stan_data,
    warmup = 500,
    iter = 3000,
    chains = 4,
    cores = 4,
    thin = 1,
    seed = 123
)

print(fit, probs = c(0.025, 0.975))

e <- rstan::extract(fit)
a <- tibble(control = e$a[, 1], resistant = e$a[, 2]) %>%
    mutate(par = "a", id = 1:n())
b <- tibble(control = e$b[, 1], resistant = e$b[, 2]) %>%
    mutate(par = "b", id = 1:n())

dl <- bind_rows(a, b) %>%
    pivot_longer(-c(par, id))

ggplot(dl, aes(x = value, color = name, fill = name)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~par, scales = "free") +
    mydarktheme
ggsave("plots/plot2.png", height = 6, width = 10)

dl %>%
    pivot_wider(names_from = par, values_from = value) %>%
    ggplot(aes(a, b, color = name)) +
    geom_point(alpha = 0.1) +
    facet_wrap(~name, scales = "free") +
    mydarktheme +
    guides(color = "none")
ggsave("plots/plot3.png", height = 5, width = 10)

pred <-
    expand.grid(dose = 10^seq(-2, 4, by = 0.1), pop = c(1, 2)) %>%
    as_tibble()

pred$pred <- NA
pred$predu95 <- NA
pred$predl95 <- NA
for (i in 1:nrow(pred)) {
    pop <- pred$pop[i]
    dose <- pred$dose[i]
    logit <- e$a[, pop] + e$b[, pop] * log(dose)
    p <- 1 / (1 + exp(-logit))
    qs <- quantile(p, c(0.025, 0.5, 0.975))
    pred$predl95[i] <- qs[1]
    pred$pred[i] <- qs[2]
    pred$predu95[i] <- qs[3]
}

pred %>%
    mutate(pop = c("control", "resistant")[pop]) %>%
    ggplot(aes(x = dose, y = pred, fill = pop)) +
    geom_point(data = d, aes(dose, dead / total, color = pop)) +
    geom_line(aes(color=pop)) +
    geom_ribbon(aes(ymin = predl95, ymax = predu95), alpha = 0.5) +
    mydarktheme +
    coord_cartesian(expand=FALSE, ylim=c(0,1), xlim=c(1e-2, 1e4)) +
    xlab("omethoate (mg/L)") +
    ggtitle("Bayesian") +
    scale_x_log10()
ggsave("plots/plot4.png", height = 6, width = 10)


# As a side note it is generally better to use matrix notation so that the stan code need not change for minor variations

x <- model.matrix(~ pop / I(log(dose)) - 1, data = db)
stan_data_mat <- list(
    N = nrow(db),
    y = db$outcome,
    K = ncol(x),
    x = x
)

fit_mat <- stan(
    file = "model_matrix.stan",
    data = stan_data_mat,
    warmup = 500,
    iter = 3000,
    chains = 4,
    cores = 4,
    thin = 1,
    seed = 123
)

print(fit_mat, probs = c(0.025, 0.975))

########### SECOND BLOG POST START HERE ##############

-coef(m1)[1] / coef(m1)[2]

library(MASS)
# we can simply run dose.p(m1, cf = c(1,3) or more verbosely
cf <- c(1, 3) # The terms in the coefficient vector giving the intercept and coefficient of (log-)dose
p <- 0.5 # Probabilities at which to predict the dose needed.
eta <- family(m1)$linkfun(p)
b <- coef(m1)[cf]
x.p <- (eta - b[1L]) / b[2L]
names(x.p) <- paste("p = ", format(p), ":", sep = "")
pd <- -cbind(1, x.p) / b[2L]
SE <- sqrt(((pd %*% vcov(m1)[cf, cf]) * pd) %*% c(1, 1))

# LD50
exp(c(x.p - 1.96 * SE, x.p + 1.96 * SE))

# bayesian
post <- extract(fit_mat)
logx50 <- -post$beta[, 1] / post$beta[, 3]
exp(mean(logx50))
exp(quantile(logx50, c(0.025, 0.975)))

# LC50 using maximum liklihood
# control
con <- dose.p(m1, c(1, 3), p = 0.5)
exp(con)
t95 <- qnorm(0.025, lower.tail = FALSE)
exp(c(con - t95 * attr(con, "SE"), con + t95 * attr(con, "SE")))

# resistant
res <- dose.p(m1, c(2, 4), p = 0.5)
res
exp(c(res - t95 * attr(con, "SE"), res + t95 * attr(con, "SE")))

# resistance factor
exp(res[1]) / exp(con[1])

m1_full <- glm(cbind(dead, alive) ~ pop * log(dose),
    family = binomial("logit"),
    data = d
)
m1_no_pop <- glm(cbind(dead, alive) ~ log(dose),
    family = binomial("logit"),
    data = d
)
anova(m1_full, m1_no_pop, test = "LRT")

# control LC50
exp(mean(-post$beta[, 1] / post$beta[, 3]))

# resistant LC50
exp(mean(-post$beta[, 2] / post$beta[, 4]))

# resistance factor with 95% credible intervals
x50_diff <- exp(-post$beta[, 2] / post$beta[, 4] - 
    -post$beta[, 1] / post$beta[, 3])
mean(x50_diff)
quantile(x50_diff, c(0.025, 0.975))
