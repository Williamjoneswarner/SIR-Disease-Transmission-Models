1. Basic SIR Model Simulation (Varying Beta and Gamma)
Purpose: Simulates the dynamics of the SIR model using various transmission rates (beta) and recovery rates (gamma).
Key Parameters: Beta values range from 0.1 to 0.9, gamma values increase in 0.2 increments.
Output: Plots showing the evolution of susceptible, infected, and recovered populations over time for each combination of beta and gamma.
Visualization: Line plots, with different colors representing compartments, are created using ggplot2 and faceted by parameter combinations.
2. SIR Model with Effective Reproduction Number (R_t)
Purpose: Extends the SIR model by incorporating the effective reproduction number (R_t) over time.
Key Parameters: Same beta and gamma variations as the first model.
Output: Additional line plot for R_t is introduced, calculated as a function of beta, gamma, and the susceptible population at each time step.
Visualization: R_t is plotted alongside the susceptible, infected, and recovered curves using ggplot2, with faceting for different beta and gamma combinations.
3. EpiModel SIR Network Simulation (Multiple Beta and Gamma)
Purpose: Uses the EpiModel package to simulate SIR dynamics on a network structure (1000 nodes) for different beta and gamma values.
Key Parameters: Beta values from 0.1 to 0.9, gamma values in increments of 0.2.
Output: Epidemic data is extracted, and the number of infected individuals over time is plotted.
Visualization: Faceted plots show results for each beta-gamma pair, using a network simulation framework.
4. Single Beta and Gamma Simulation (EpiModel)
Purpose: Runs the SIR model on a network with fixed beta (0.3) and gamma (0.1).
Key Parameters: Single values for beta (0.3) and gamma (0.1).
Output: Time-series plot of infected individuals is generated.
Visualization: A simple line plot of the number of infected individuals over time for the fixed parameter values.
5. Multiple Simulations with Mean and Standard Deviation (10 Iterations)
Purpose: Repeats the network simulation (with beta = 0.3 and gamma = 0.1) ten times to analyze variability.
Key Parameters: Beta = 0.3, Gamma = 0.1; 10 iterations.
Output: Mean and standard deviation of infected individuals across all simulations are computed.
Visualization: Line plot with error bars representing the mean infected individuals over time, plus standard deviation.
