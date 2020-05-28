sim_EM.Rdata
variable	description
alpha		fixed effect
beta		random effect
end_time	finish time of estimation
epsilon		error term
gamma		determine whether beta is zero or not
m		number of fixed effect
n		sample size
n_gamma_1	number of random are set to be 1
p		number of random effect
sigma2_beta	variance of beta given gamma = 1
sigma2_epsilon	error variance
start_time	start time of estimation
parameter	result from variational EM
--L		L(q) in updating the variation parameter
--m		variational mean
--s2		variational variance
--p		probability of being included in the model
--s2_eps	estimated error variance
--s2_beta	estimated beta variance
--pi		estimated probability of gamma = 1
--outer		expected complete data likelihood

sim_LASSO.Rdata
variable	description
LASSO_CV	Result from LASSO regression produced by cv.glmnet from package glmnet