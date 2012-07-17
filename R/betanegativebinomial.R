
dbetanbinom <- function(x, size, shape1, shape2, log = FALSE)
    .External("actuar_do_dpq", "dbetanbinom", x, size, shape1, shape2, log)
pbetanbinom <- function(q, size, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pbetanbinom", q, size, shape1, shape2, lower.tail, log.p)
qbetanbinom <- function(p, size, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qbetanbinom", p, size, shape1, shape2, lower.tail, log.p)
rbetanbinom <- function(n, size, shape1, shape2)
    .External("actuar_do_random", "rbetanbinom", n, size, shape1, shape2)