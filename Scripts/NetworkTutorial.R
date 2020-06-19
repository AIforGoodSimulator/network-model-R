library(EpiModel)


# first step is to initialize empty networks with desired number of
# individiuals/nodes (no edges for now)
# In this case the network does not need to be directed as relationships
# are bidirectional
# there is also a "race" attribute that can be used to define
# node categories, note this package was built predominantly to
# simulate STIs

number_of_individuals = 100
nw = network.initialize(n = number_of_individuals, directed = FALSE)

#  The next line creates a vertex attribute called  race , 
#which will have categories 0 and 1: there are 500 nodes of 
# Race 0 and 500 nodes of Race 1.
nw = set.vertex.attribute(nw, "race", 
                          rep(0:1, each = number_of_individuals/2))

## Next we specify the model parameters:
# 2 steps:
#     1. Define formation formula: the formation is the sum of 
#         several factors such as number of edges, as well as matching nodes
#         and councurrent number numbers 
formation <- ~edges+nodefactor("race")+nodematch("race")+concurrent

#     2. target statistics (from empirical data egocentric data):
#         These are things like mean degree (per catergory),
#         [we actually given expected number of edges; 
#         mean_degree*# of nodes]

mean_degree = 0.5
stat1 = (number_of_individuals/2)*mean_degree #number of expected edges
mean_degree_race1 = 0.75
mean_degree_race0 = 0.25
stat2 = (number_of_individuals/2)*mean_degree_race1
# same-race edge
same_race_edge = 0.9
stat3 = stat1*same_race_edge
# number of concurrent edges (concurrent degree)
stat4 = number_of_individuals*0.1
target.stats = c(stat1,stat2, stat3, stat4) 
# c(expecteed edges per category, expected edges Race1,
# Exected edges of same race, #nodes with concurrent degree)

## Neext define dissolution coefficient from given formula
# this coefficient represents how fast links are dissolved/broken
# Under the hood this is turned into a logit for estimation
coef.diss = dissolution_coefs(dissolution = ~offset(edges),
                              duration = 25)
coef.diss

#################################################################
##                          Model fit                          ##
#################################################################
# Here we use the netest function which is a wrapper around the 
# estimation function fro the ergm and tergm packages
# the function requires nw, formation, targe.stats and coef.diss

# the edapprox argument determines whether to use direct or indirect
# (approximation) method to find the network. Direct uses tergm 
# indirect uses ergm
est1 <- netest(nw, formation, target.stats,
               coef.diss, edapprox = TRUE)

## Moodel diagnostics
# nsims: number of model to simulate
# nsteps: 
dx <- netdx(est1, nsims = 15, nsteps = 500, 
            nwstats.formula = ~edges +
              nodefactor("race", base = 0)+
              nodematch("race")+
              concurrent)

dx

## plot diagnostics
par(mar = c(3,2,1,1))
# dashed line is expected numeber
# x-axis is timesteps
plot(dx)

# can plot other things
par(mfrow = c(1,2))
# duration of link
plot(dx, type = "duration")
# dissolution of link (1/duration)
plot(dx, type = "dissolution")

## EPIDEMIC SIMULATION

#' inf.prob: infection probability given contact
#' act.rate: number of acts (or contacts) per time unit
#' within a partnership
#' rec.rate: recovery rate

param = param.net(inf.prob = 0.1, 
                  act.rate = 5, 
                  rec.rate = 0.02) #average duration of disease
                                  # is 50 timesteps (1/rec.rate)

# for initial condition swe can use the argument i.num to set the initial number of 
# infected at the start or pass in a vector with disease status for eacch of the nodes
# EpiModel stores the individual-level disease status as a vector of lower-case 
# letters: “s” for susceptible, “i” for infected, and “r” for recovered. 
# In the assignment below we specify that Race 0 has a baseline prevalence oof 10% (assigned randomly)
status.vector <- c(rbinom(number_of_individuals/2, 1, 0.1), rep(0, number_of_individuals/2)) 
status.vector <- ifelse(status.vector == 1, "i", "s")
init <- init.net(status.vector = status.vector)
#next we specify the control type
# eppi.by will gives us results split by a certain parameter
control <- control.net(type = "SIS", nsteps = 500,
                       nsims = 10, epi.by = "race")

# our net is now parametrised and we can model the epidemic
sim1 = netsim(est1, #our net
              param, #our net params
              init, #initial state
              control) # control simulation settings

# we can now print sim1 and we'll get summaries about the simulation runs
# the model output are things like: s.numm - number of susceptiple, 
# i.num, number of infected for a SIS simulation

# we can get the summary statistics at a given timepoint(e.g. t = 500)
# like this:
summary(sim1, at =500)

# get simulation results in dataframe, default outputs mean outputs but can be changed
results = as.data.frame(sim1, out= "mean")

# one can get by-time network by:
nw_at_time1 = get_network(sim1, sim = 1)

# some details are stored for each transmission that occurs
#' the inf column shows the ID of transmitting node
tranmission_results = get_transmat(sim1, sim = 1)

#plotting the results
par(mfrow = c(1,1), mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim1)

# can also get per simulation plot
plot(sim1, mean.line = F, qnts = F, sim.lines = T)

# can also plot disease incidence and recoveries with flows
plot(sim1, y = c("si.flow", "is.flow"), qnts = 1, legend = T)

#by race or category
plot(sim1, y = c("i.num.race0", "i.num.race1"), legend = T)

# the default plot results is "epi" but we can also get a plot
# of the network

par(mfrow = c(1,2), mar = c(0,0,1,0))
plot(sim1, type = "network", at = 1,
     col.status = T, main = "Prevalence at t1")
plot(sim1, type = "network", at = 500,
     col.status = T, main = "Prevalence at t500")

############################################################################
############################################################################
###                                                                      ###
###                     DEPENDENT BIPARTITE SI MODEL                     ###
###                                                                      ###
############################################################################
############################################################################

# NETWORK ESTIMATION AND DIAGNOSTICS
#'  bipartite graph (or bigraph) is a graph whose vertices can be divided into two 
#'  disjoint and independent sets U and V such 
#'  that every edge connects a vertex in U to one in V. 

num.m1 <- 500 # could be females
num.m2 <- 500 # and males
nw <- network.initialize(num.m1 + num.m2,
                         bipartite = num.m1,
                         directed = FALSE)

# below is the degree distribution aka probabilies for every possible degree value
# we allow for, e.g:
# 40% prob no partner, 55% prob 1 partner, 4% prob 2 partners, 1% prob 3 partners
deg.dist.m1 = c(0.4, 0.55, 0.04, 0.01)
deg.dist.m2 = c(0.48, 0.41, 0.08, 0.03)

# we also specify an average degree tthat applies to everyone
# this is "simulated" with a poisson distribution btwn 0 and 2 and we sum up thee cumulative mass
# for degree large than 2 (strictly)
pois.dists = c(dpois(0:2, lambda = 0.66), ppois(2, lambda = 0.66, lower = F))
# this mean degree overdoes deg0 and deg2 a little

# let's cehck if the pois.dists matches our numbers
par(mar = c(3, 3, 2, 1), mfrow = c(1, 1)) 
cols <- transco(RColorBrewer::brewer.pal(4, "Set1"), 0.8)
barplot(cbind(deg.dist.m1, deg.dist.m2, pois.dists),beside = FALSE, ylim = c(0, 1), col = cols)
legend("topright", legend = paste0("deg", 3:0), pch = 15, col = rev(cols), bg = "white")

check_bip_degdist(num.m1, num.m2,
                  deg.dist.m1, deg.dist.m2)

## Model fit
# the formation basically tells the model which stats to track
# in this case it is trakcing number of edges (number of relationships in this case)
# the degree of nodes of mode 1 (b1degree) with degree 0 and 1
# the degree of nodes of mode 2 (b2degree) with degree 0 and 1
formation <- ~edges + b1degree(0:1) + b2degree(0:1) 
# target stats are:
#' overall number of partnerships at one point in time(aka #of edges)
#' numbers of nodes in the first mode with no partners
#' or only one partner
#' similar degree terms for the second mode nodes
#' CAN BE extracted from the check_bip_degdist
target.stats <- c(330, 200, 275, 240, 205)

coef.diss = dissolution_coefs(dissolution = ~offset(edges), 
                              duration = 25, 
                              d.rate = 0.005)
# d.rate expected death rate per unit of time (in this case:
# after 1000 unit of times we'd expect 5 deaths)
# d.rate adjusts the crude dissolution coefficient

#' THESE ARE THE ESSENTIALS TO BUILD A MODEL BASICALLY
est2 = netest(nw, formation, target.stats, coef.diss)

# Model diagnostics

dx = netdx(est2, nsims =5, nsteps = 500)
dx

# The  stats  argument of the  plot.netdx  
# function allows for sub- setting the statistics 
# for better visualization.
plot(dx, stats = "edges")

# Similar to epidemic plots, different elements of the plots
# may be toggled to highlight different components of the 
# diagnostics.
plot(dx, stats = c("b1deg1", "b2deg1"),
     sim.lines = TRUE, sim.lwd = 0.4,  qnts = FALSE,
     mean.lwd = 4, mean.smooth = FALSE)

## Epidemic simulation

# because we are making a dependent network, disease dynamics
# andn network dynamics are linked, we need to explicitly 
# change the params to account for this
param = param.net(inf.prob = 0.3,  # infection probability
                  inf.prob.m2 = 0.1,
                  a.rate = 0.005, # arrival rate
                  a.rate.m2 = NA, # aka birth rate in this case
                  ds.rate = 0.005, #departure rate for susceptible
                  ds.rate.m2 = 0.005, # rate
                  di.rate = 0.005, # departure rate for infected
                  di.rate.m2 = 0.005) # rate
?param.net

init = init.net(i.num = 50,
                i.num.m2 = 50)

control = control.net(type = "SI",
                      nsims = 5, 
                      nsteps = 500,
                      nwstats.formula = ~edges+meandeg,
                      delete.nodes = T)

# so now we get a model where the total population changes as
# births and deaths are allowed
sim2 = netsim(est2, param, init, control)

## DIAGNOSTICS

# CHECK THAT death/births did not affect (bias) expectted edges or mean
# degree too much
plot(sim2, type = "formation", plots.joined = FALSE)

# another plot; number of infected/susceptibles per mode(male/female)
plot(sim2, popfrac = F)
