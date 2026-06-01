# emmdt requires i0 >= 1

    Code
      emmdt(s0 = 10, i0 = 0, beta = 1.5, gamma = 1)
    Condition
      Error in `emmdt()`:
      ! i0 >= 1 is not TRUE

# emmdt requires s0 >= 1

    Code
      emmdt(s0 = 0, i0 = 1, beta = 1.5, gamma = 1)
    Condition
      Error in `emmdt()`:
      ! s0 >= 1 is not TRUE

