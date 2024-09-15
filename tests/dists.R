#!/usr/bin/env Rscript
distName = commandArgs(TRUE)[1]
dat = scan("stdin", quiet = TRUE)
switch(distName,
    gamma = {
        alpha = dat[1]
        beta = dat[2]
        dat =  dat[3:length(dat)]
        res = ks.test(dat, "pgamma", shape=alpha, scale=beta)
    },
    exponential = {
        lambda = dat[1]
        dat =  dat[2:length(dat)]
        res = ks.test(dat, "pexp", rate = lambda)
    }
)

if(res$p.value < 0.01){
    quit(status = 1)
}
res = Box.test(dat, type = "Box-Pierce")
if(res$p.value < 0.01){
    quit(status = 2)
}
res = Box.test(dat, type = "Ljung-Box")
if(res$p.value < 0.01){
    quit(status = 3)
}

quit(status = 0)
