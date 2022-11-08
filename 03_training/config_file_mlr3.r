# Config file
# ====
setwd('/mnt/netapp2/Store_uni/home/ulc/co/dfe/git/IBEROBDIA/')
# Data path
ExperimentName = 'second-experiment'
input.dir.path = '03_training/toRun/second-experiment/'
out.dir.path = '03_training/results/'
out.path = paste(out.dir.path, ExperimentName, '/', sep = '')
exec.dir.path = '03_training/exec/'
out_slurm.dir.path = '03_training/outs/'

if (dir.exists(out.path) == FALSE) {
  dir.create(out.path)
  message(paste('Creating ', ExperimentName, ' directory!'))
}

# Path of the algorithms
path.algs = '03_training/models/'
pattern.algs = '_mlr3.r'
out.filename.glmnet='glmnet.rds'
out.filename.rf='rf.rds'
out.filename.xgboost = "xgboost.rds"

# Cesga resources
part = 'thinnodes'
qos = 'default'
time = '02:00:00'
mem = '120G'
nodes = 1
ntasks = 24

# General Machine-Learns parameters
seed = 1342
workers = 24

# Glmnet hyperparameters
gl.alpha = c(0,1)
gl.nlambda = 500
gl.lambda.min.ratio = 10^(-2)
gl.s = c(0,1)

# Random Forest hyperparameters
rf.mtry = c(2,8)
rf.ntree = 400
rf.nodesize = c(1,3)

# XGboost hyperparameters
xg.booster = c("gbtree", "gblinear", "dart")
xg.alpha = c(0, 1)
xg.eta = c(0, 1)
xg.lambda = c(0.2, 0.8)
xg.gamma = c(0.2, 0.8)
xg.max_depth = c(3, 18)

# Cross validation parameters
## Supra parameters
cv.inner = 'Holdout'
cv.outer = 'RepCV'

## Sub parameters
## Estos parámetros no los estás pasando a los modelos 
predict = c('both') # train or both
iters = 10
reps = 5
folds = 3
strat= TRUE