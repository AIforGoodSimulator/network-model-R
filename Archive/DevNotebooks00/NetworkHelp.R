## Loading libraries
library(EpiModel)

## Initialising network 
n = 50
nw = network.initialize(n = n, directed = FALSE)

prop.isobox = 0.43 # proportion of housing instances that are isoboxes
prop.tent = 1 - prop.isobox # proportion of housing instances that are tents

stopifnot(round(n*prop.isobox+n*prop.tent)==n)

tent.capacity = 4
iso.capacity = 10

num_of_iso = round(n*prop.isobox/iso.capacity) # number of isoboxes
num_of_tents = round(n*prop.tent/tent.capacity) # number of tents

num_in_iso = num_of_iso*iso.capacity
num_in_tents = n - num_in_iso # num_of_tents*tent.capacity 

iso_ids = 1:num_of_iso
tent_ids = 1:num_of_tents 

housing_iso = paste0("iso",
                     apportion_lr(
                       vector.length = num_in_iso,
                       values = iso_ids,
                       proportions = rep(1/(num_of_iso), num_of_iso)
                     )
)

housing_tents = paste0("tent",
                       apportion_lr(
                         vector.length = num_in_tents,
                         values = tent_ids,
                         proportions = rep(1/(num_of_tents), num_of_tents)
                       )
)

housing = c(housing_iso,housing_tents)

# residence vector to keep track of those in tents and isos
residence = c(rep("iso", num_in_iso), rep("tent", num_in_tents))

# set attributes
nw = set.vertex.attribute(nw, "housing", housing)

#################################################################
##               Setting up formation statistics               ##
#################################################################

formation = ~edges + 
  offset(nodematch("housing", diff = F))

# default degrees in housing units as per current occupancy
iso.default.degree = mean(table(housing)[iso_ids])
tent.default.degree = mean(table(housing)[tent_ids])

# number of external contacts per person in average
external.contacts = 0

mean_degree.iso =  iso.default.degree + external.contacts
mean_degree.tent =  tent.default.degree + external.contacts

# calculate mean degree
mean_degree = (num_in_iso*mean_degree.iso+
                 num_in_tents*mean_degree.tent)/n

expected.edges = n*mean_degree/2

target.stats = expected.edges



#################################################################
##              Setting up dissolution statistics              ##
#################################################################

d.rate = 0
coef.diss = dissolution_coefs(dissolution = ~offset(edges)+offset(nodematch("housing", diff = F)),
                              duration = c(1,1e6),
                              d.rate = d.rate)

#################################################################
##                  Network estimation (MCMC)                  ##
#################################################################

est <- netest(nw,
              formation,
              target.stats,
              coef.diss,
              coef.form = Inf,
              set.control.ergm = control.ergm(MCMLE.maxit = 500)
)
summary(est)

#################################################################
##                         Diagnostics                         ##
#################################################################

cores = parallel::detectCores()-1

dx = netdx(est,
           nsims = 9,
           nsteps = 120, # simulating 6 months
           ncores = cores,
           nwstats.formula = ~edges+
             nodematch("housing", diff = FALSE),
           dynamic = F)
dx

##################################################################
##                          Simulation                          ##
##################################################################

## random parameters, extracted from ?netsim
param <- param.net(inf.prob = 0.3, inf.prob.m2 = 0.15,
                   rec.rate = 0.02, rec.rate.m2 = 0.02,
                   a.rate = 0.002, a.rate.m2 = NA,
                   ds.rate = 0.001, ds.rate.m2 = 0.001,
                   di.rate = 0.001, di.rate.m2 = 0.001,
                   dr.rate = 0.001, dr.rate.m2 = 0.001)
init <- init.net(i.num = 10, i.num.m2 = 10,
                 r.num = 0, r.num.m2 = 0)
control <- control.net(type = "SIR", nsteps = 10, nsims = 2)

sim = netsim(est, param, init, control)

#################################################################
##                         Network GIF                         ##
#################################################################

library(animation)

nw_object = get_network(sim)
ani.record(reset = TRUE)  # clear history before recording
for (at in c(1:30)){
  
  net_at = network.collapse(nw_object, at = at)
  graph = intergraph::asIgraph(net_at)
  adj = igraph::as_adjacency_matrix(graph, sparse = F)
  
  colnames(adj) = igraph::V(graph)$housing
  rownames(adj) = igraph::V(graph)$housing
  
  adj = adj[order(rownames(adj)), order(colnames(adj))]
  
  pheatmap::pheatmap(adj,
                     color = c("grey50","black"),
                     border_color = "white",
                     angle_col = 45,
                     angle_row = 45,
                     fontsize = 6,
                     legend_breaks = c(0,1),
                     legend = F,
                     cluster_rows = F,
                     cluster_cols = F,
                     show_rownames = ifelse(n<100, T, F),
                     show_colnames = ifelse(n<100, T, F))
  ani.record()
}

oopts = ani.options(interval = 0.5)
ani.replay()

saveHTML(ani.replay(), img.name = "record_plot")
saveGIF(ani.replay(), movie.name = "networkdynamic.gif")

