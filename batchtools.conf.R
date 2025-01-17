#cluster.functions = makeClusterFunctionsSSH(list(Worker$new("localhost", ncpus = 20, max.load = 40)))
cluster.functions = makeClusterFunctionsMulticore(ncpus = 20)