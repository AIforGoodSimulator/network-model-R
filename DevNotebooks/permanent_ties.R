library(EpiModel)

dur = 1e6 # edge duration for nodes with the same value for the color attribute
n=60
nw = network.initialize(n = n,
                        directed = FALSE)

nw = set.vertex.attribute(nw, "color", c(rep("red", n/3),
                                         rep("blue", n/3),
                                        rep("green", n/3)))

formation = ~edges + 
  offset(nodematch("color", diff = FALSE))

max.edges = n*(n-1)/2

edges = 0*max.edges

target.stats = c(edges)

coef.diss = dissolution_coefs(dissolution = ~offset(edges)+
                                offset(nodematch("color", diff = FALSE)),
                              duration = c(2, dur))

est <- netest(nw,
              formation,
              target.stats,
              coef.diss,
              coef.form = Inf,
              set.control.ergm = control.ergm(MCMLE.maxit = 500)
)

cores = parallel::detectCores()-1

dx <- netdx(est,
            nsims = 1e3,
            nsteps = 90,
            ncores = cores,
            dynamic = FALSE,
            nwstats.formula = ~edges + nodematch("color", diff = FALSE),
            set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6))

print(dx)
## model
d.rate = 0

param <- param.net(inf.prob = 0.3, 
                   rec.rate = 0.02,
                   a.rate = 0,
                   ds.rate = d.rate,#0.001,
                   di.rate = d.rate,# 0.001,
                   dr.rate = d.rate# 0.001,
)

init <- init.net(i.num = 10,
                 r.num = 0
)

control <- control.net(type = "SIR",
                       nsteps = 30,
                       nsims = 3)

# Simulate the model with new network fit
sim <- netsim(est, param, init, control)

### GIF
library(animation)

nw_object = get_network(sim)
ani.record(reset = TRUE)  # clear history before recording
for (at in c(1:30)){
  
  net_at = network.collapse(nw_object, at = at)
  graph = intergraph::asIgraph(net_at)
  adj = igraph::as_adjacency_matrix(graph, sparse = F)
  
  colnames(adj) = igraph::V(graph)$color
  rownames(adj) = igraph::V(graph)$color
  
  adj = adj[order(rownames(adj)), order(colnames(adj))]
  pheatmap::pheatmap(adj,
                     color = c("grey50","black"),
                     border_color = "white",
                     angle_col = 45,
                     angle_row = 45,
                     fontsize = 6,
                     legend_breaks = c(0,1),
                     legend = FALSE,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     show_rownames = FALSE,
                     show_colnames = FALSE)
  ani.record()
}

oopts = ani.options(interval = 1)
# ani.replay()
#saveHTML(ani.replay(), img.name = "network")
saveGIF(ani.replay(), movie.name = paste0("network",dur,".gif"))


nw1 <- get_network(sim)
df <- as.data.frame(nw1)

color <- nw1 %v% "color"
summary(df$duration[which(color[df$head] == color[df$tail])])
summary(df$duration[which(color[df$head] != color[df$tail])])




