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
	* `exposure()`: calls `FromXToExposed()` to model the transitions 
