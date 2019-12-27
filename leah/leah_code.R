library(mirt)
library(MASS)
dat <- expand.table(LSAT7)

mod1 <- mirt(dat, 1, "2PL", constrain = list(c(1, 5, 9, 13, 17)), SE = TRUE)
mod2 <- mirt(dat, 1, "2PL", SE = TRUE)
mod3 <- mirt(dat, 1, "3PL", SE = TRUE, parprior = list(c((1:ncol(dat)) * 4 - 1, 'norm', -1.8, 1)))


tseq <- seq(-4, 4, by = .25)

probs1 <- sapply(1:5, function(i) probtrace(extract.item(mod1, i), tseq)[, 2])
probs2 <- sapply(1:5, function(i) probtrace(extract.item(mod2, i), tseq)[, 2])
probs3 <- sapply(1:5, function(i) probtrace(extract.item(mod3, i), tseq)[, 2])

## confidence envelopes

vcv1 <- extract.mirt(mod1, "vcov")
mv1 <- extract.mirt(mod1, "parvec")
vcv2 <- extract.mirt(mod2, "vcov")
mv2 <- extract.mirt(mod2, "parvec")
vcv3 <- extract.mirt(mod3, "vcov")
mv3 <- extract.mirt(mod3, "parvec")

pars1 <- mvrnorm(n = 1000, mu = mv1, Sigma = vcv1)
pars2 <- mvrnorm(n = 1000, mu = mv2, Sigma = vcv2)
pars3 <- mvrnorm(n = 1000, mu = mv3, Sigma = vcv3)

probs1_1 <- t(apply(pars1, 1, function(p) 1 / (1 + exp(-p[1] * tseq - p[2]))))
probs1_2 <- t(apply(pars1, 1, function(p) 1 / (1 + exp(-p[1] * tseq - p[3]))))
probs1_3 <- t(apply(pars1, 1, function(p) 1 / (1 + exp(-p[1] * tseq - p[4]))))
probs1_4 <- t(apply(pars1, 1, function(p) 1 / (1 + exp(-p[1] * tseq - p[5]))))
probs1_5 <- t(apply(pars1, 1, function(p) 1 / (1 + exp(-p[1] * tseq - p[6]))))

probs2_1 <- t(apply(pars2, 1, function(p) 1 / (1 + exp(-p[1] * tseq - p[2]))))
probs2_2 <- t(apply(pars2, 1, function(p) 1 / (1 + exp(-p[3] * tseq - p[4]))))
probs2_3 <- t(apply(pars2, 1, function(p) 1 / (1 + exp(-p[5] * tseq - p[6]))))
probs2_4 <- t(apply(pars2, 1, function(p) 1 / (1 + exp(-p[7] * tseq - p[8]))))
probs2_5 <- t(apply(pars2, 1, function(p) 1 / (1 + exp(-p[9] * tseq - p[10]))))

probs3_1 <- t(apply(pars3, 1, function(p){
    c <- 1 / (1 + exp(-p[3]))
    c + (1 - c) / (1 + exp(-p[1] * tseq - p[2]))}))
probs3_2 <- t(apply(pars3, 1, function(p){
    c <- 1 / (1 + exp(-p[6]))
    c + (1 - c) / (1 + exp(-p[4] * tseq - p[5]))}))
probs3_3 <- t(apply(pars3, 1, function(p){
    c <- 1 / (1 + exp(-p[9]))
    c + (1 - c) / (1 + exp(-p[7] * tseq - p[8]))}))
probs3_4 <- t(apply(pars3, 1, function(p){
    c <- 1 / (1 + exp(-p[12]))
    c + (1 - c) / (1 + exp(-p[10] * tseq - p[11]))}))
probs3_5 <- t(apply(pars3, 1, function(p){
    c <- 1 / (1 + exp(-p[15]))
    c + (1 - c) / (1 + exp(-p[13] * tseq - p[14]))}))

# item 1
curve(1 / (1 + exp(-mv1[1] * x - mv1[2])), xlim = c(-3, 3), ylim = c(0, 1))
points(tseq, apply(probs1_1, 2, quantile, probs = .025), lty = 2, type = 'l')
points(tseq, apply(probs1_1, 2, quantile, probs = .975), lty = 2, type = 'l')

curve(1 / (1 + exp(-mv2[1] * x - mv2[2])), add = TRUE, col = 4)
points(tseq, apply(probs2_1, 2, quantile, probs = .025), lty = 2, type = 'l', col = 4)
points(tseq, apply(probs2_1, 2, quantile, probs = .975), lty = 2, type = 'l', col = 4)

curve(1 / (1 + exp(-mv3[3])) + (1 - 1 / (1 + exp(-mv3[3])) ) / (1 + exp(-mv3[1] * x - mv3[2])), add = TRUE, col = 2)
points(tseq, apply(probs3_1, 2, quantile, probs = .025), lty = 2, type = 'l', col = 2)
points(tseq, apply(probs3_1, 2, quantile, probs = .975), lty = 2, type = 'l', col = 2)

# item 2
curve(1 / (1 + exp(-mv1[1] * x - mv1[3])), xlim = c(-3, 3), ylim = c(0, 1))
points(tseq, apply(probs1_2, 2, quantile, probs = .025), lty = 2, type = 'l')
points(tseq, apply(probs1_2, 2, quantile, probs = .975), lty = 2, type = 'l')

curve(1 / (1 + exp(-mv2[3] * x - mv2[4])), add = TRUE, col = 4)
points(tseq, apply(probs2_2, 2, quantile, probs = .025), lty = 2, type = 'l', col = 4)
points(tseq, apply(probs2_2, 2, quantile, probs = .975), lty = 2, type = 'l', col = 4)

curve(1 / (1 + exp(-mv3[6])) + (1 - 1 / (1 + exp(-mv3[6])) ) / (1 + exp(-mv3[4] * x - mv3[5])), add = TRUE, col = 2)
points(tseq, apply(probs3_2, 2, quantile, probs = .025), lty = 2, type = 'l', col = 2)
points(tseq, apply(probs3_2, 2, quantile, probs = .975), lty = 2, type = 'l', col = 2)

# item 3
curve(1 / (1 + exp(-mv1[1] * x - mv1[4])), xlim = c(-3, 3), ylim = c(0, 1))
points(tseq, apply(probs1_3, 2, quantile, probs = .025), lty = 2, type = 'l')
points(tseq, apply(probs1_3, 2, quantile, probs = .975), lty = 2, type = 'l')

curve(1 / (1 + exp(-mv2[5] * x - mv2[6])), add = TRUE, col = 4)
points(tseq, apply(probs2_3, 2, quantile, probs = .025), lty = 2, type = 'l', col = 4)
points(tseq, apply(probs2_3, 2, quantile, probs = .975), lty = 2, type = 'l', col = 4)

curve(1 / (1 + exp(-mv3[9])) + (1 - 1 / (1 + exp(-mv3[9])) ) / (1 + exp(-mv3[7] * x - mv3[8])), add = TRUE, col = 2)
points(tseq, apply(probs3_3, 2, quantile, probs = .025), lty = 2, type = 'l', col = 2)
points(tseq, apply(probs3_3, 2, quantile, probs = .975), lty = 2, type = 'l', col = 2)

# item 4
curve(1 / (1 + exp(-mv1[1] * x - mv1[5])), xlim = c(-3, 3), ylim = c(0, 1))
points(tseq, apply(probs1_4, 2, quantile, probs = .025), lty = 2, type = 'l')
points(tseq, apply(probs1_4, 2, quantile, probs = .975), lty = 2, type = 'l')

curve(1 / (1 + exp(-mv2[7] * x - mv2[8])), add = TRUE, col = 4)
points(tseq, apply(probs2_4, 2, quantile, probs = .025), lty = 2, type = 'l', col = 4)
points(tseq, apply(probs2_4, 2, quantile, probs = .975), lty = 2, type = 'l', col = 4)

curve(1 / (1 + exp(-mv3[12])) + (1 - 1 / (1 + exp(-mv3[12])) ) / (1 + exp(-mv3[10] * x - mv3[11])), add = TRUE, col = 2)
points(tseq, apply(probs3_4, 2, quantile, probs = .025), lty = 2, type = 'l', col = 2)
points(tseq, apply(probs3_4, 2, quantile, probs = .975), lty = 2, type = 'l', col = 2)

# item 5
curve(1 / (1 + exp(-mv1[1] * x - mv1[6])), xlim = c(-3, 3), ylim = c(0, 1))
points(tseq, apply(probs1_5, 2, quantile, probs = .025), lty = 2, type = 'l')
points(tseq, apply(probs1_5, 2, quantile, probs = .975), lty = 2, type = 'l')

curve(1 / (1 + exp(-mv2[9] * x - mv2[10])), add = TRUE, col = 4)
points(tseq, apply(probs2_5, 2, quantile, probs = .025), lty = 2, type = 'l', col = 4)
points(tseq, apply(probs2_5, 2, quantile, probs = .975), lty = 2, type = 'l', col = 4)

curve(1 / (1 + exp(-mv3[15])) + (1 - 1 / (1 + exp(-mv3[15])) ) / (1 + exp(-mv3[13] * x - mv3[14])), add = TRUE, col = 2)
points(tseq, apply(probs3_5, 2, quantile, probs = .025), lty = 2, type = 'l', col = 2)
points(tseq, apply(probs3_5, 2, quantile, probs = .975), lty = 2, type = 'l', col = 2)




####

p1a <- function(p, theta){
    c(1 / (1 + exp(-p[1] * theta - p[2])),
      1 / (1 + exp(-p[1] * theta - p[3])),
      1 / (1 + exp(-p[1] * theta - p[4])),
      1 / (1 + exp(-p[1] * theta - p[5])),
      1 / (1 + exp(-p[1] * theta - p[6])))
}
p2a <- function(p, theta){
    c(1 / (1 + exp(-p[1] * theta - p[2])),
      1 / (1 + exp(-p[3] * theta - p[4])),
      1 / (1 + exp(-p[5] * theta - p[6])),
      1 / (1 + exp(-p[7] * theta - p[8])),
      1 / (1 + exp(-p[9] * theta - p[10])))
}
p3a <- function(p, theta){
    c1 <- 1 / (1 + exp(-p[3]))
    c2 <- 1 / (1 + exp(-p[6]))
    c3 <- 1 / (1 + exp(-p[9]))
    c4 <- 1 / (1 + exp(-p[12]))
    c5 <- 1 / (1 + exp(-p[15]))
    c(c1 + (1 - c1) / (1 + exp(-p[1] * theta - p[2])),
      c2 + (1 - c2) / (1 + exp(-p[4] * theta - p[5])),
      c3 + (1 - c3) / (1 + exp(-p[7] * theta - p[8])),
      c4 + (1 - c4) / (1 + exp(-p[10] * theta - p[11])),
      c5 + (1 - c5) / (1 + exp(-p[13] * theta - p[14])))
}

f1a <- function(eta, pars, targetp){
    p <- p1a(pars, eta)
    sqrt(mean((p - targetp)^2))
}
f2a <- function(eta, pars, targetp){
    p <- p2a(pars, eta)
    sqrt(mean((p - targetp)^2))
}
f3a <- function(eta, pars, targetp){
    p <- p3a(pars, eta)
    sqrt(mean((p - targetp)^2))
}


res_1a <- lapply(tseq, function(t){
    targetp <- p1a(p = mv1, theta = t)
    apply(pars1, 1, function(pars){
        optimize(f = f1a, interval = c(-6, 6), pars = pars, targetp = targetp)
    })
})
res_2a <- lapply(tseq, function(t){
    targetp <- p2a(p = mv2, theta = t)
    apply(pars2, 1, function(pars){
        optimize(f = f2a, interval = c(-6, 6), pars = pars, targetp = targetp)
    })
})
res_3a <- lapply(tseq, function(t){
    targetp <- p3a(p = mv3, theta = t)
    apply(pars3, 1, function(pars){
        optimize(f = f3a, interval = c(-6, 6), pars = pars, targetp = targetp)
    })
})


res_1a_min <- t(sapply(res_1a, function(x) sapply(x, function(y) y$minimum)))
res_1a_obj <- t(sapply(res_1a, function(x) sapply(x, function(y) y$objective)))
res_2a_min <- t(sapply(res_2a, function(x) sapply(x, function(y) y$minimum)))
res_2a_obj <- t(sapply(res_2a, function(x) sapply(x, function(y) y$objective)))
res_3a_min <- t(sapply(res_3a, function(x) sapply(x, function(y) y$minimum)))
res_3a_obj <- t(sapply(res_3a, function(x) sapply(x, function(y) y$objective)))

save.image("illustration 1.RData")


#### ####

library(tidyverse)
library(mirt)

load("illustration 1.RData")


## Table 1

round(mv1, 2)[c(2, 3, 4, 5, 6, 1)]
round(sqrt(diag(vcv1)), 2)[c(2, 3, 4, 5, 6, 1)]
round(vcv1, 3)[c(2, 3, 4, 5, 6, 1), c(2, 3, 4, 5, 6, 1)]

round(mv2, 2)[c(2, 4, 6, 8, 10, 1, 3, 5, 7, 9)]
round(sqrt(diag(vcv2)), 2)[c(2, 4, 6, 8, 10, 1, 3, 5, 7, 9)]
round(vcv2, 3)[c(2, 4, 6, 8, 10, 1, 3, 5, 7, 9), c(2, 4, 6, 8, 10, 1, 3, 5, 7, 9)]

cbind(round(mv3, 2)[c(2, 5, 8, 11, 14, 1, 4, 7, 10, 13, 3, 6, 9, 12, 15)],
      round(sqrt(diag(vcv3)), 2)[c(2, 5, 8, 11, 14, 1, 4, 7, 10, 13, 3, 6, 9, 12, 15)])
round(vcv3, 3)[c(2, 5, 8, 11, 14, 1, 4, 7, 10, 13, 3, 6, 9, 12, 15), c(2, 5, 8, 11, 14, 1, 4, 7, 10, 13, 3, 6, 9, 12, 15)]

# check order of eta, etheta
mean(apply(res_1a_min, 2, function(x) cor(x, tseq, method = "spearman")) == 1) # 100%
mean(apply(res_2a_min, 2, function(x) cor(x, tseq, method = "spearman")) == 1) # 97.5%
mean(apply(res_3a_min, 2, function(x) cor(x, tseq, method = "spearman")) == 1) # 82.4%


LSATdat <- tibble("theta" = tseq,
                  "i1m" = 1 / (1 + exp(-mv2[1] * tseq - mv2[2])),
                  "i2m" = 1 / (1 + exp(-mv2[3] * tseq - mv2[4])),
                  "i3m" = 1 / (1 + exp(-mv2[5] * tseq - mv2[6])),
                  "i4m" = 1 / (1 + exp(-mv2[7] * tseq - mv2[8])),
                  "i5m" = 1 / (1 + exp(-mv2[9] * tseq - mv2[10])),
                  "i1l" = apply(probs2_1, 2, quantile, probs = .025),
                  "i2l" = apply(probs2_2, 2, quantile, probs = .025),
                  "i3l" = apply(probs2_3, 2, quantile, probs = .025),
                  "i4l" = apply(probs2_4, 2, quantile, probs = .025),
                  "i5l" = apply(probs2_5, 2, quantile, probs = .025),
                  "i1u" = apply(probs2_1, 2, quantile, probs = .975),
                  "i2u" = apply(probs2_2, 2, quantile, probs = .975),
                  "i3u" = apply(probs2_3, 2, quantile, probs = .975),
                  "i4u" = apply(probs2_4, 2, quantile, probs = .975),
                  "i5u" = apply(probs2_5, 2, quantile, probs = .975))

vcov <- extract.mirt(mod3, "vcov")
ses <- sqrt(diag(vcov))
vcov <- round(vcov, 2)
row.names(vcov) <- NULL
colnames(vcov) <- as.character(sapply(1:5, function(x) paste0(c("a", "d", "g"), x)))
for(i in 1:ncol(vcov)) vcov[, i] <- format(vcov[, i], digits = 2)
vcov[upper.tri(vcov, diag = TRUE)] <- " "

plot_partable <- tibble(parameter = as.character(sapply(1:5, function(x) paste0(c("a", "d", "g"), x))),
                        est = round(extract.mirt(mod3, "parvec"), 2),
                        se = round(ses, 2))
plot_partable <- cbind(plot_partable, vcov)


LSATdat <- LSATdat %>% gather(type, p, -theta)
LSATdat <- LSATdat %>% mutate(item = paste("Item", substr(type, start = 2, stop = 2)),
                              bound = substr(type, start = 3, stop = 3))

plot_LSAT_CB1 <-ggplot(LSATdat, aes(theta, p, col = item, linetype = bound)) + geom_path() + facet_wrap(~item) +
    scale_linetype_manual(values=c("dashed", "solid", "dashed")) + theme_light() +
    labs(x = expression(theta), y = "Probability") +
    theme(legend.position = "none")
plot_LSAT_CB1


LSATdat2 <- tibble("theta" = tseq,
                   "M1.50" =  apply(res_1a_obj, 1, quantile, probs = .5),
                   "M1.95" =  apply(res_1a_obj, 1, quantile, probs = .95),
                   "M2.50" =  apply(res_2a_obj, 1, quantile, probs = .5),
                   "M2.95" =  apply(res_2a_obj, 1, quantile, probs = .95),
                   "M3.50" =  apply(res_3a_obj, 1, quantile, probs = .5),
                   "M3.95" =  apply(res_3a_obj, 1, quantile, probs = .95))
LSATdat2 <- LSATdat2 %>% gather(type, RMSD, -theta)
LSATdat2 <- LSATdat2 %>% mutate(model = paste0(substr(type, start = 2, stop = 2), "PM"),
                                quantile = substr(type, start = 3, stop = 5),
                                light = (quantile == ".50" & RMSD > .02) | (quantile == ".95" & RMSD > .05))


plot_LSAT_RMSD1 <- ggplot(LSATdat2, aes(theta, RMSD, col = model, linetype = quantile, size = I(1))) +
    geom_path(na.rm=TRUE) +
    theme_minimal() + geom_hline(yintercept = c(.02, .05), col = I("gray"), linetype = "dotted") +
    scale_color_grey() +
    theme(legend.position = "bottom") +
    lims(y = c(0, .15)) + labs(x = expression(theta), y = expression(Delta)) +
    guides(color = FALSE, linetype = guide_legend(override.aes = list(color = "black", fill = "white")))
plot_LSAT_RMSD1 ## panel A of figure to include

LSATdat3 <- t(sapply(tseq, function(t){
    a <- mv1[1]
    b <- mv1[-1] / a
    info1 <- sum(catR::Ii(t, it = cbind(a, b, 0, 1))$Ii)
    a <- mv2[c(1, 3, 5, 7, 9)]
    b <- -mv2[c(2, 4, 6, 8, 10)] / a
    info2 <- sum(catR::Ii(t, cbind(a, b, 0, 1))$Ii)
    a <- mv3[c(1, 4, 7, 10, 13)]
    b <- -mv3[c(2, 5, 8, 11, 14)] / a
    c <- 1 / (1 + exp(-mv3[c(3, 6, 9, 12, 15)]))
    info3 <- sum(catR::Ii(t, cbind(a, b, c, 1))$Ii)
    c(info1, info2, info3)
}))
LSATdat4 <- tibble(theta = tseq, "1PM" = LSATdat3[, 1],
                   "2PM" = LSATdat3[, 2], "3PM" = LSATdat3[, 3])
LSATdat4 <- LSATdat4 %>% gather(model, info, -theta)

plot_LSAT_info <- ggplot(LSATdat4, aes(theta, info, col = model, size = I(1))) + geom_path() +
    scale_color_grey() +
    theme_minimal() +
    labs(x = expression(theta), y = "information") + theme(legend.position = "bottom")


jpeg(file = "LSATfig.jpeg", width = 7, height = 4, units = "in", res = 200)
gridExtra::grid.arrange(plot_LSAT_RMSD1, plot_LSAT_info, nrow = 1)
dev.off()
