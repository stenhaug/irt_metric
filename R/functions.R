f_to_min <- function(theta, model, observed_probs){
    rmse(
        observed_probs,
        theta_model_to_probs_vector(theta, model)
    )
}

get_best_theta <- function(observed_probs, model){
    optimize(
        f = f_to_min,
        model = model,
        observed_probs = observed_probs,
        interval = c(-6, 6)
    )
}

theta_model_to_probs_vector <- function(theta, model){
    n_items <- length(model@Data$K)
    1:n_items %>% map_dbl(~ probtrace(extract.item(model, .), theta)[, 2])
}

model_to_probs <- function(model){
    n_items <- length(model@Data$K)

    1:n_items %>%
        map(~ probtrace(extract.item(model, .), quads)[, 2]) %>%
        as_tibble(.name_repair = ~ paste0("item_", 1:n_items)) %>%
        mutate(quads = quads) %>%
        dplyr::select(quads, everything())
}

get_draws <- function(model, n_draws){
    n_items <- length(model@Data$K)

    mvrnorm(n = n_draws, mu = extract.mirt(model, "parvec"), Sigma = extract.mirt(model, "vcov")) %>%
        as_tibble(.name_repair = ~ paste0(c("disc_", "easy_"), rep(1:n_items, each = 2))) %>%
        mutate(draw = row_number()) %>%
        gather(var, value, -draw) %>%
        separate(var, into = c("parameter", "item")) %>%
        select(draw, item, parameter, value) %>%
        arrange(draw, item, parameter) %>%
        spread(parameter, value)
}

twopl <- function(disc, easy, theta) 1 / (1 + exp(-disc * theta - easy))

rmse <- function(x, y){sqrt(mean((x - y)^2))}
