library(EpiModel)

# First define population and summary statistics

n = 20000 # individuals in camp
# initialize empty net
nw = network.initialize(n = n, directed = FALSE)

# set age attribute according to real age breakdown
#' Age breakdown within the camp (20,000 people)
#   50+ 4.7%      class 3
#   19-49 54.8%   class 2
#   0-18 40.5%    class 1
num.age_cl1 = (40.5/100)*n
num.age_cl2 = (54.8/100)*n
num.age_cl3 = (4.7/100)*n

stopifnot(num.age_cl1+num.age_cl2+num.age_cl3 == n)

nw = set.vertex.attribute(nw, "age_cl", 
                          c(rep("young", num.age_cl1),
                            rep("adults", num.age_cl2),
                            rep("old", num.age_cl3)))

# formation to make network 
formation <- ~edges+nodefactor("age_cl", base = 1)+nodematch("age_cl", diff = FALSE)+concurrent
#' In this case we'll need
#' --> edges:[+1 statistic] expected # of edges in average (proxy for the mean degree)
#' --> nodefactor("age_cl", base = 1) [+2 statistics] this will ignore the stat for class with most members (adults), this
#' is recommended because each tie will count twice once for each node. nodefactor() allows us to set 
#' the expected number of edges per class/group. base = 1, sets class 1 to the base level (which I am at
#' this moment unsure of what it is, may be avg degree may not be)
#' --> nodematch(age_cl, diff = FALSE):[+1 statistic] expected number of edges of individuals from the same group, 
#' the diff flag allows you to set a single statistic that applies to everyone or you could also make 
#' it per class. 
#' --> concurrency: [+1 statistic] number of nodes expected to have 2 or more edges at any given time
#' 
mean_degree = 2 # how many connections a person has in average
stat1 = (n/3)*mean_degree #number of expected edges per age_class in general
#' [!!!!!] these need to average out to the mean degree above
mean_degree_age_cl1 = 2 # young
mean_degree_age_cl2 = 3.5 # adults
mean_degree_age_cl3 = 0.5 # old
stat2 = num.age_cl3*mean_degree_age_cl3 
stat3 = num.age_cl1*mean_degree_age_cl1
# same-race edge, let's say that only 50% of connections are between people of same age class
same_race_edge = 0.5
stat4 = stat1*same_race_edge # number of same class edges per person (see above)
# number of concurrent edges (concurrent degree), number of same connections at the same time
stat5 = n*0.5
target.stats = c(stat1,stat2, stat3, stat4, stat5) 

coef.diss = dissolution_coefs(dissolution = ~offset(edges),
                              duration = 100) # avg duration of each node -> 100 steps
coef.diss

# I found that the coefs are in a order that I was not expecting
est1 <- netest(nw, formation, target.stats,
               coef.diss, edapprox = T, verbose = F)

summary(est1)



