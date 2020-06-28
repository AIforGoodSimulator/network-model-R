# Custom modules for Network Model

## States and Tokens

In this network model individuals can be in 7 seven different states:
* Susceptible - S
* Exposed - E
* Infected - I
* Quarantined - Q
* Require hospitalisation - H
* Recovered - R
* Fatality (Death) - F

## Exposure module

* __Purpose:__ Deal with transition from susceptible to exposed. Susceptible individuals can become exposed upon "succesful" contact with either exposed, infected or quarantined individuals.
* __Functions in use:__ `FromXToExposed()`, `custom_discord_edgelist()`  and `exposure()`
* __How it works:__ 
	* `custom_discord_edgelist()`: The original `discord_edgelist()` function gets a list of edges that connect susceptibles and infected individuals to model the transition from susceptible to infected upon contact between an infected and a suscpetible individual. This _custom_ function has been modified to return list of edges between susceptible and any other state that is fed into the function.
	* `FromXToExposed()`: models the transition between any susceptible and exposed upon contact with either quarantined, exposed or infected individuals. It makes use of `custom_discord_edgelist()`. The transmission probabality for each of the possible links (S<->E, S<->I, S<->Q) is passed to this function. Then the outcome is drawn for a Bernouilli distribution with the transmission probability as the p.
	* `exposure()`: calls `FromXToExposed()` to model the transitions before-mentioned.

* __Faults/Quirks:__ 
	* Interactions are modelled sequentially (i.e. first interactins between S and E individuals, then between S and Q and finally between S and I. This biases the transitions and makes people more susceptible than the parameters suggest. See (this github issue)[https://github.com/statnet/EpiModel/issues/417] for a discussion about it.

## Infection module

* __Purpose:__ model transition between exposed and infected state. 
* __Functions in use:__  `infect()`
* __How it works:__ For each exposed individual, whether they become infected is drawn from a Bernouilli distribution with _p_ equal to the `ei.rate` (exposed to infected rate).
* __Faults/Quirks:__ assumes everyone gets infected at the same rate, regardless of age/gender attributes, which may not be true

## Quarantining module

* __Purpose:__ model transition between infected and quarantined state.
* __Functions in use:__  `quarantinining()`
* __How it works:__ For each infected individual, whether they become quarantined is drawn from a Bernouilli distribution with _p_ equal to the `iq.rate` (infected to quarantined rate).
* __Faults/Quirks:__ 
	* Quarantined individuals are not actually quarantined in the actual network (i.e. they are still making new ties, though they decreased their contact rate). If this feature was to be implemented it would require tying the epidemic modelling to the network formation dynamics which is quite a bit more intricate. Calling the Q state quarantined may perhaps be misleading. 

## Hospitalisation module

* __Purpose:__ model transition between quarantined or infected states and requiring-hospitalisation state.
* __Functions in use:__  `requireHospitalization()`
* __How it works:__ For each quarantined individual, whether they require hospitalisation (i.e. become a critical patient) is drawn from a Bernouilli distribution with _p_ equal to the `ih.rate`. If the network nodes have an age attribute, the `ih.rate` is overwritten to values drawn from COVID19 research where the hospitalisation rate is proportional to age. For those that quarantine the rate of transition to requiring hospitalisation is not modulated by age. 
* __Faults/Quirks:__ 
	* It seems very arbitrary that the age-modularisation applies to the I->H transition but not to the Q->H one

## Recovery module

* __Purpose:__ model transition between either the infected, quarantined, requiring-hospitalisation states to the recovered state.
* __Functions in use:__  `recover()`
* __How it works:__ For each eligible individual, whether they recover from the infection is drawn from a Bernouilli distribution with _p_ equal to the rate of transition to the recovered state from the specific state of origin (i.e. `qr.rate`, `hr.rate` and `ir.rate`).
* __Faults/Quirks:__ 
	* Again we find a many-to-one transitions where the transitions are modelled sequentially instead of in parallel (with state-transitions matrices). This again effectively biases the transitions dynamics, making it more likely that individuals recover than what the parameters suggest.
	* Individuals of all ages and attributes are equally likely to recover, which is probably not true as older individuals or individuals with co-morbidities may struggle more to recover.

## Fatality module

* __Purpose:__ model transition between requiring hospitalisation and fatality state.
* __Functions:__ `fatality()`
* __How it works:__ For each individual requiring hospitalisation, whether they die of the infection is drawn from a Bernouilli distribution with _p_ equal to the `hf.rate`. In the same way as the requiring-hospitalisation rate, the `hf.rate` is age dependent with parameters based on research. In order for individuals to die they would become critical first, hence passing through the requiring hospitalisation rate. The `hf.rate` is further modulated based on the hospital capacity and the amount of days spent requiring hospitalisation.

## Other utility modules

I have added some further functionality to some of the core `EpiModel` functions mostly to be able to extracts information relating to time spent in each state. 
These include: `get_prev.net()`, `initilize.net()`, `saveout.net()` and `netsim()`.


