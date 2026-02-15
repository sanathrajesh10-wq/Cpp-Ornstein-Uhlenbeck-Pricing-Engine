# C-Ornstein-Uhlenbeck-Pricing-Engine
OU_process_2 is currently just the MLE itself, since until the numerical instability issue is resolved, you can not use its result for the actual numerical discretization or the analytical solution.
monte_carlo.cpp will be adjusted so that it uses the mean and standard deviation defined by the MLE's variable.
There exist some libraries, like map_ex, xin some of the files that serve no purpose and will be removed at a later date.
