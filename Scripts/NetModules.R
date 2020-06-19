## Net functions

infect <- function(dat, at) {
  #' This function describes how people get infected, thus becoming
  #' exposed (in an SEIR model)
  
  # select those that are active (aka alive)
  active <- dat$attr$active
  # get their status (s, i, r...)
  status <- dat$attr$status
  # get the network
  nw <- dat$nw
  
  # get the ids of the infected
  idsInf <- which(active == 1 & status == "i")
  # get the number of active nodes
  nActive <- sum(active == 1)
  
  # get the number of infected (and set as eligible?)
  nElig <- length(idsInf)
  nInf <- 0
  
  if (nElig > 0 && nElig < nActive) {
    #discordant edgelist – a matrix of ID numbers of active dyads 
    # in the network in which one member of the dyad is susceptible
    # and the other is infected
    del <- discord_edgelist(dat, at)
    if (!(is.null(del))) { # if del exists
      # set transmission probability to infection prob
      del$transProb <- dat$param$inf.prob
      # get act rate (~contact rate)
      del$actRate <- dat$param$act.rate
      # calculate prob of infection from standard formula
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      # draw samples from binomial (Bernouilli) distribution with 
      # final prob
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      # select those that have transmitted according to the trial
      del <- del[which(transmit == 1), ]
      # Get new infected
      idsNewInf <- unique(del$sus)
      # update number of infected
      nInf <- length(idsNewInf)
      
      if (nInf > 0) {
        # set new infected as "e", for exposed
        dat$attr$status[idsNewInf] <- "e"
        # store infection time
        dat$attr$infTime[idsNewInf] <- at
      }
    }
  }
  
  if (at == 2) {
    # set flow of susceptiple to infected for at 1 and at 2
    # remember at 1 the model is initialised nothing actually happens
    dat$epi$se.flow <- c(0, nInf)
  }
  else {
    dat$epi$se.flow[at] <- nInf
  }
  dat$nw <- nw
  return(dat)
}



## New disease progression module
progress <- function(dat, at) {
  # This function dictates transition between states, in this case
  # tranisiton between exposed and infected and between infected and recovered
  # this function makes the recovery.net default function obsolete
  
  ####### get the eligible
  # active nodess
  active <- dat$attr$active
  # get status too (s, i, r, e)
  status <- dat$attr$status
  
  # get the ei rate and ir rate (Exposed to infected and infected to recovered)
  ei.rate <- dat$param$ei.rate
  ir.rate <- dat$param$ir.rate
  
  # also gonna make de.rate and dr.rate
  de.rate = dat$param$de.rate
  dr.rate = dat$param$dr.rate
  
  ## E to I progression
  nInf <- 0 # initialise number of infected variable
  # get IDs of those that are eligible for infection (aka active nodes that are exposed)
  # in this case exposure is a prerequisite to be infected
  idsEligInf <- which(active == 1 & status == "e")
  # get number of people that are eligible for infection
  nEligInf <- length(idsEligInf)
  
  
  if (nEligInf > 0) { # if anyone is eligible
    # pick those that will be infected accordging to a bernouilli trial
    vecInf <- which(rbinom(nEligInf, 1, ei.rate) == 1) # like the toss of a coin
    if (length(vecInf) > 0) { # if by the trial results somepeople were picked
      idsInf <- idsEligInf[vecInf] 
      # get the id for those and transition them from "e" to "i"
      nInf <- length(idsInf)
      status[idsInf] <- "i"
    }
  }
  
  ## I to R progression
  # same as before but f infected to recovered transition
  nRec <- 0
  idsEligRec <- which(active == 1 & status == "i")
  nEligRec <- length(idsEligRec)
  
  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, ir.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nRec <- length(idsRec)
      status[idsRec] <- "r"
    }
  }
  
  # dat$attr$status <- status
  
  # I should maybe modify the departure.net function insteead of this not sure
  
  ## E to D progression 
  # same as before but for exposed to dead transition
  nDE <- 0
  idsEligDE <- which(active == 1 & status == "e")
  nEligDE <- length(idsEligDE)
  
  if (nEligDE > 0) {
    vecDE <- which(rbinom(nEligDE, 1, de.rate) == 1)
    if (length(vecDE) > 0) {
      idsDE <- idsEligDE[vecDE]
      nDE <- length(idsDE)
      
      status[idsDE] = "d"
      
      dat$attr$active[idsDE] <- 0
      dat$attr$exitTime[idsDE] <- at
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                    v = idsDE, deactivate.edges = TRUE)
    }
  }
  
  ## R to D progression
  # same as before but for recovered to dead transition
  nDR <- 0
  idsEligDR <- which(active == 1 & status == "r")
  nEligDR <- length(idsEligDR)
  
  if (nEligDR > 0) {
    vecDR <- which(rbinom(nEligDR, 1, dr.rate) == 1)
    if (length(vecDR) > 0) {
      idsDR <- idsEligDR[vecDR]
      nDR <- length(idsDR)
      
      status[idsDR] = "d"
      
      dat$attr$active[idsDR] <- 0
      dat$attr$exitTime[idsDR] <- at
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                    v = idsDR, deactivate.edges = TRUE)
    }
  }
  
  dat$attr$status <- status
  
  if (at == 2) {
    dat$epi$ei.flow <- c(0, nInf)
    dat$epi$ir.flow <- c(0, nRec)
    dat$epi$de.flow <- c(0, nDE)
    dat$epi$dr.flow <- c(0, nDR)
    dat$epi$e.num <- c(0, sum(active == 1 & status == "e"))
    dat$epi$r.num <- c(0, sum(active == 1 & status == "r"))
  }
  else {
    dat$epi$ei.flow[at] <- nInf
    dat$epi$ir.flow[at] <- nRec
    dat$epi$de.flow[at] <- nDE
    dat$epi$dr.flow[at] <- nDR
    dat$epi$e.num[at] <- sum(active == 1 & status == "e")
    dat$epi$r.num[at] <- sum(active == 1 & status == "r")
  }
  
  return(dat)
}