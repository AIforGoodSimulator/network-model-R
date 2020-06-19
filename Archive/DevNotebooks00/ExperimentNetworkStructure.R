## Figuring out net structure


## Loading libraries
library(EpiModel)

SimulateNetwork = function(n = 200, diss.time = 30){
  ## Initialising network
  ## Initialising network 
  n = n
  nw = network.initialize(n = n, directed = FALSE)
  
  prop.isobox = 0.43 # proportion of housing instances that are isoboxes
  prop.tent = 1 - prop.isobox # proportion of housing instances that are tents
  
  stopifnot(round(n*prop.isobox+n*prop.tent)==n)
  
  tent.capacity = 4 
  iso.capacity = 10
  #number of tents/isoboxes in camp, inferred from number of people, proportion of tents/isoboxes and tent/isobox capacity 
  num_of_tents = round(n*prop.tent/tent.capacity) 
  num_of_iso = round(n*prop.isobox/iso.capacity)
  
  # number of individuals in tents/isoboxes (that fit in pre-allocated number of tents/isoboxes)
  # remaining individuals will be placed in tents later on in the code. 
  #' These vector are place holders for the initial capacity given the current figures,
  #' They are unimportant when it comes to the acctualy network formation where the
  #' residence and the housing vectors (see below) are used
  num_in_tents = num_of_tents*tent.capacity 
  num_in_iso = num_of_iso*iso.capacity
  
  #' Housing attribute vector, corresponding to which tent or which isobox individual belong to
  housing = c(rep(paste0("tent", 1:num_of_tents),tent.capacity),
              rep(paste0("isobox", 1:num_of_iso), iso.capacity))
  #' [1] "tent1"  "tent2"  "tent3"  "tent4"  "tent5" 
  #'  "tent6"  "tent7"  "tent8"  "tent9"  "tent10"
  
  #' Residence attribute indicating whether individuals are in tents or in isoboxes
  residence = c(rep("tent", num_in_tents),
                rep("isobox", num_in_iso))
  
  #' If the number of allocated people is less than the number of people
  #' in the network allocate remaining individuals in a new tent and 
  #' update number of tents as well as number of individuals in tents
  #' This number in practice has been observed to be 4 or less, I haven't checked its behaviour though
  if (length(housing)<n){
    remaining = n-length(housing)
    print(paste(remaining, "nodes were left unassigned. Assigning them to new tent."))
    housing[length(housing)+(1:remaining)] = paste0("tent", num_of_tents+1)
    residence[length(residence)+(1:remaining)] = paste0("tent")
    num_of_tents = num_of_tents+1
    num_in_tents = num_in_tents+remaining
  }
  
  # if residence/housing vector are too long, cut them to size n
  residence = residence[1:n]
  housing = housing[1:n]
  
  # set attributes
  nw = set.vertex.attribute(nw, "housing", housing)
  nw = set.vertex.attribute(nw, "residence", residence)
  
  # set vectors to their network counterpart
  residence = get.vertex.attribute(nw, "residence")
  housing = get.vertex.attribute(nw, "housing")
  
  #################################################################
  ##               Setting up formation statistics               ##
  #################################################################
  
  formation <- ~edges+
    nodematch("housing", diff = FALSE) # amount of interaction with same class nodes
  
  #' number of ties outside, if set to 0, assumes mean 
  #' degree per housing is equal to housing capacity (unrealistic)
  external.friends = 5
  #' 1st term in RHS is equivalent to the "empirical" isobox/tent degree (accounting for initially unallocated
  #' individuals)
  #' 
  iso.occupancy = mean(table(housing)[grepl("isobox", names(table(housing)))])
  tent.occupancy = mean(table(housing)[grepl("tent", names(table(housing)))])
  
  mean_degree.iso =  iso.occupancy + external.friends
  mean_degree.tent =  tent.occupancy + external.friends
  
  # calculate mean degree empirically: degree_iso * num_in_iso + 
  mean_degree = (sum(residence == "isobox")*mean_degree.iso+
                   sum(residence == "tent")*mean_degree.tent)/
    length(residence) 
  # set expected number of edges
  edges = n*mean_degree/2 
  
  # set nodematch(housing) statistic as percentage of total number of edges
  mean_degree.match = (sum(residence == "isobox")*iso.occupancy+
                         sum(residence == "tent")*tent.occupancy)/
    length(residence) 
  
  max.in.housing.edges = (iso.capacity*(iso.capacity-1)/2) * num_of_iso + 
    (tent.capacity*(tent.capacity-1)/2)*num_of_tents
  
  housing.match =  max.in.housing.edges/2#n*mean_degree.match/2#edges*0.5  
  
  target.stats = c(edges,
                   housing.match)
  
  #################################################################
  ##              Setting up dissolution statistics              ##
  #################################################################
  
  d.rate = 0.0001
  coef.diss = dissolution_coefs(dissolution = ~offset(edges),
                                duration = diss.time, 
                                d.rate = d.rate) # this correspond to external deaths
  coef.diss
  
  #################################################################
  ##                  Network estimation (MCMC)                  ##
  #################################################################
  
  est1 <- netest(nw,
                 formation,
                 target.stats,
                 coef.diss,
                 edapprox = T,
                 verbose = F,
                 set.control.ergm = control.ergm(MCMLE.maxit = 100))
  
  summary(est1)
  
  #################################################################
  ##                         Diagnostics                         ##
  #################################################################
  
  cores = parallel::detectCores()-1
  
  dx = netdx(est1,
             nsims = 9,
             nsteps = 120, # simulating 6 months
             ncores = cores,
             nwstats.formula = ~edges+
               nodematch("housing", diff = FALSE))
  # dx
  
  perc.error.edge = dx$stats.table.formation[1,3]
  perc.error.nodematch = dx$stats.table.formation[2,3]
  
  return(list(perc.error.edge, perc.error.nodematch))
  
}

gridsearch = expand.grid(n = c(50,100,200,300,500,750,1000), 
                         dtime = c(7,14,21,28, 35, 42, 49, 60))

results = sapply(1:nrow(gridsearch), function(x){
  tempgrid = unlist(gridsearch[x,])
  print(tempgrid)
  SimulateNetwork(tempgrid[1], tempgrid[2])
})


results = data.frame(results)
results = t(results)

colnames(results) = c("EdgeError", "MatchError")

results = cbind(results, gridsearch)

library(tidyverse)

res2 = results %>% select(-MatchError) %>% pivot_wider(names_from = n, values_from = EdgeError)

res2 = data.matrix(res2)

res2 = res2[,-1]

library(pheatmap)
library(grid)

setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9,
                                                         name="vp",
                                                         just=c("right","top"))),
        action="prepend")
pheatmap(res2, cluster_rows = F, cluster_cols = F)
setHook("grid.newpage", NULL, "replace")
grid.text("Number of nodes, n", y=-0.07, gp=gpar(fontsize=16))
grid.text("Average edge duration", x=-0.07, rot=90, gp=gpar(fontsize=16))
grid.text("Pct. Difference target-to-estimate", y=1, x = 0.5, gp=gpar(fontsize=20))




