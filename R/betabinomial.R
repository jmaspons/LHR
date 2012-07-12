
dbetabinom <- function(x, size, shape1, shape2, log = FALSE)
    .External("dbetabinom", x, size, shape1, shape2, log)
pbetabinom <- function(q, size, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
    .External("pbetabinom", q, size, shape1, shape2, lower.tail, log.p)
qbetabinom <- function(p, size, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
    .External("qbetabinom", p, size, shape1, shape2, lower.tail, log.p)
rbetabinom <- function(n, size, shape1, shape2)
    .External("rbetabinom", n, size, shape1, shape2)