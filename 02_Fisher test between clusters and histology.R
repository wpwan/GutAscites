groups <- matrix(c(4, 6, 0,
                   3, 4, 1,
                   3, 1, 1),
                 nrow = 3,
                 dimnames = list(c("NOS", "SRC", "N/A"),
                                 c("A", "B", "C")))
groups
fisher.test(groups)  # p-value = 0.4883






















