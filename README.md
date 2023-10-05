# m-ctp_heuristic

This repo contains the heuristic used to compute the initial upper bound for the m-ctp.

The heuristic works as follows:
0 - best <- \empty

1 - Identify a subset, denoted as S, of optional facilities that collectively cover all customers. Once all customers are covered, eliminate any redundant facilities.

2 - Apply the HGS-CVRP algorithm to the set S ∪ M to obtain a solution.

3 - Create a new subset, S', by introducing perturbations to S. The perturbation involves selecting a random facility from S and removing the k closest facilities to it, which includes the selected facility itself.

4 - Apply the HGS-CVRP algorithm on S' ∪ M.

5 - If the cost of S' is lower than that of S, replace S with S', and return to step 3. Else, best <- S, and return to step 1.

6 - Repeat these steps until the time limit is reached and return best.

