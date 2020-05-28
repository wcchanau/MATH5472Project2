GWAS_result.Rdata
variable	description
n		sample size
p		number of random effect
Result		Variational EM result
--		Four lists corresponding to each phenotype
---		Under each of the four list:
---L		L(q) in updating the variation parameter
---m		variational mean
---s2		variational variance
---p		probability of being included in the model
---s2_eps	estimated error variance
---s2_beta	estimated beta variance
---pi		estimated probability of gamma = 1
---outer		expected complete data likelihood

GWAS_LASSO.Rdata
variable	description
LASSO		LASSO result
--		Four lists corresponding to each phenotype
---		Result from LASSO regression produced by cv.glmnet from package glmnet