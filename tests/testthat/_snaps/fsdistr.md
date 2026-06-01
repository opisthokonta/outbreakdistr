# fsdistr validates ip_model

    Code
      fsdistr(s0 = 5, i0 = 1, beta = 1.5, ip_model = "invalid", ip_params = 1)
    Condition
      Error in `validate_ip_model_params()`:
      ! ip_model must be "constant", "exponential", or "gamma".

# fsdistr validates ip_params type

    Code
      fsdistr(s0 = 5, i0 = 1, beta = 1.5, ip_model = "constant", ip_params = "bad")
    Condition
      Error in `validate_ip_model_params()`:
      ! Infectious period parameters must be numeric.

