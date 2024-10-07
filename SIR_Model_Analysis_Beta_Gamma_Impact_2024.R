#SIR Model Dynamics: Impact of Varying Transmission and Recovery Rates on Disease Spread - single pathogen/infection event
# Load required packages
install.packages("deSolve")
install.packages("ggplot2")
library(deSolve)
library(ggplot2)

# SIR model function
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

# Parameters for beta and gamma
beta_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)  # Specific beta values
gamma_values <- seq(0.1, 0.9, by = 0.2)  # Gamma values from 0.1 to 0.9 in 0.2 increments

# Initial state
initial_state <- c(S = 990, I = 10, R = 0)
N <- sum(initial_state)  # Total population

# Time settings
time_seq <- seq(0, 100, by = 1)  # Time sequence from 0 to 100 days

# Prepare an empty data frame to store results
results <- data.frame()  # Initialize an empty data frame

# Loop over each beta value
for (beta in beta_values) {
  # Loop over each gamma value
  for (gamma in gamma_values) {
    
    # Define parameters for the current combination
    parameters <- c(beta = beta, gamma = gamma)
    
    # Run the simulation
    output <- ode(y = initial_state, times = time_seq, func = sir_model, parms = parameters)
    
    # Convert output to data frame and add beta and gamma columns
    output_df <- as.data.frame(output)
    colnames(output_df) <- c("Time", "Susceptible", "Infected", "Recovered")
    output_df$Beta <- beta  # Add beta value
    output_df$Gamma <- gamma  # Add gamma value
    output_df$Parameter_Combination <- paste("Beta =", beta, ", Gamma =", gamma)  # Combine parameters
    
    # Combine results with previous data
    results <- rbind(results, output_df)  # Append the current result to the results data frame
  }
}

# Plotting results using facet wrap
ggplot(results, aes(x = Time)) +
  geom_line(aes(y = Susceptible, color = "Susceptible"), size = 1) +
  geom_line(aes(y = Infected, color = "Infected"), size = 1) +
  geom_line(aes(y = Recovered, color = "Recovered"), size = 1) +
  labs(title = "SIR Model Dynamics for Different Beta and Gamma Combinations",
       y = "Number of Individuals",
       x = "Time (days)") +
  scale_color_manual(name = "Compartment", values = c("Susceptible" = "blue", 
                                                      "Infected" = "red", 
                                                      "Recovered" = "green")) +
  facet_wrap(~ Parameter_Combination, scales = "free_y") +  # Create facets for each parameter combination
  theme_minimal() +
  theme(legend.position = "bottom")  # Move legend to the bottom for clarity


#### with Rt

# Load necessary libraries
library(deSolve)
library(ggplot2)

# SIR model function
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I / N  # Change in Susceptible
    dI <- beta * S * I / N - gamma * I  # Change in Infected
    dR <- gamma * I  # Change in Recovered
    
    return(list(c(dS, dI, dR)))  # Return only the state derivatives
  })
}

# Parameters for beta and gamma
beta_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)  # Specific beta values
gamma_values <- seq(0.1, 0.9, by = 0.2)  # Gamma values from 0.1 to 0.9 in 0.2 increments

# Initial state
initial_state <- c(S = 990, I = 10, R = 0)
N <- sum(initial_state)  # Total population

# Time settings
time_seq <- seq(0, 100, by = 1)  # Time sequence from 0 to 100 days

# Prepare an empty data frame to store results
results <- data.frame()  # Initialize an empty data frame

# Loop over each beta value
for (beta in beta_values) {
  # Loop over each gamma value
  for (gamma in gamma_values) {
    
    # Define parameters for the current combination
    parameters <- c(beta = beta, gamma = gamma)
    
    # Run the simulation
    output <- ode(y = initial_state, times = time_seq, func = sir_model, parms = parameters)
    
    # Convert output to data frame and add beta and gamma columns
    output_df <- as.data.frame(output)
    colnames(output_df) <- c("Time", "Susceptible", "Infected", "Recovered")
    output_df$Beta <- beta  # Add beta value
    output_df$Gamma <- gamma  # Add gamma value
    
    # Create a combined title for facets
    output_df$Title <- paste("Beta =", beta, ", Gamma =", gamma)
    
    # Combine results with previous data
    results <- rbind(results, output_df)  # Append the current result to the results data frame
  }
}

# Calculate R as a new column based on current number of susceptible individuals
results$R_t <- results$Beta / results$Gamma * (results$Susceptible / N)  # Effective reproduction number

# Plotting results using facet wrap
ggplot(results, aes(x = Time)) +
  geom_line(aes(y = Susceptible, color = "Susceptible"), size = 1) +
  geom_line(aes(y = Infected, color = "Infected"), size = 1) +
  geom_line(aes(y = Recovered, color = "Recovered"), size = 1) +
  geom_line(aes(y = R_t * N, color = "R_t"), size = 1, linetype = "dashed") +  # R_t scaled back to total population
  labs(title = "SIR Model Dynamics with Effective R Calculation",
       y = "Number of Individuals / R_t",
       x = "Time (days)") +
  scale_color_manual(name = "Compartment", values = c("Susceptible" = "blue", 
                                                      "Infected" = "red", 
                                                      "Recovered" = "green",
                                                      "R_t" = "purple")) +
  facet_wrap(~ Title, scales = "free_y") +  # Use the combined title for facets
  coord_cartesian(ylim = c(0, 1000)) +  # Cap y-axis at 1000
  theme_minimal() +
  theme(legend.position = "bottom")  # Move legend to the bottom for clarity

#EpiModel
install.packages("EpiModel")  # Run this line only if EpiModel is not installed

# Load necessary libraries
library(EpiModel)
library(ggplot2)
library(dplyr)

# Function to run the SIR model for given beta and gamma
run_sir_model <- function(beta, gamma) {
  # Define the network with 1000 nodes
  nw <- network.initialize(1000, directed = FALSE)
  
  # Estimate the network model (simple network with constant degree)
  formation <- ~edges
  target.stats <- 500  # Average number of edges in the network
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation, target.stats, coef.diss)
  
  # Set the model parameters
  param <- param.net(inf.prob = beta, act.rate = 1, rec.rate = gamma)
  
  # Define the initial conditions
  init <- init.net(i.num = 10, r.num = 0)  # 10 infected, 0 recovered initially
  
  # Define the control settings
  control <- control.net(type = "SIR", nsteps = 100, nsims = 1)
  
  # Run the simulation
  out <- netsim(est, param, init, control)
  
  # Extract the epidemic data from the epi object and return
  df <- as.data.frame(out$epi)
  df$time <- 1:nrow(df)  # Create a time variable
  df$Beta <- beta  # Add beta value
  df$Gamma <- gamma  # Add gamma value
  
  return(df)
}

# Define beta and gamma values
beta_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)
gamma_values <- seq(0.1, 0.9, by = 0.2)

# Prepare an empty list to store results
results <- list()

# Loop over each combination of beta and gamma
for (beta in beta_values) {
  for (gamma in gamma_values) {
    sim_result <- run_sir_model(beta, gamma)
    
    # Append to results list
    results[[paste("Beta =", beta, "Gamma =", gamma)]] <- sim_result
  }
}

# Combine all results into one data frame
all_results <- bind_rows(results, .id = "Simulation")

# Check the structure of the resulting data
str(all_results)  # Ensure 'i.num' is present in the dataset

# Renaming the columns of all_results to reflect compartments
colnames(all_results)[2:4] <- c("pop_size", "susceptible", "infected", "recovered")

# Check the structure again to confirm
head(all_results)

# Now, plot the results using ggplot
ggplot(all_results, aes(x = time, y = infected, color = as.factor(Beta))) +
  geom_line(linewidth = 1) +  # Use linewidth instead of size
  facet_wrap(~ Beta + Gamma, scales = "free_y", 
             labeller = label_both) +  # Automatically label each plot with Beta and Gamma
  labs(title = "SIR Model Simulation with Varying Beta and Gamma",
       y = "Infected Individuals (i.num)",
       x = "Time Steps",
       color = "Beta Values") +
  theme_minimal() +
  theme(legend.position = "bottom")

#single beta gamma value in epimodel

# Load necessary libraries
library(EpiModel)
library(ggplot2)
library(dplyr)

# Function to run the SIR model for given beta and gamma
run_sir_model <- function(beta, gamma) {
  # Define the network with 1000 nodes
  nw <- network.initialize(1000, directed = FALSE)
  
  # Estimate the network model (simple network with constant degree)
  formation <- ~edges
  target.stats <- 500  # Average number of edges in the network
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation, target.stats, coef.diss)
  
  # Set the model parameters
  param <- param.net(inf.prob = beta, act.rate = 1, rec.rate = gamma)
  
  # Define the initial conditions
  init <- init.net(i.num = 10, r.num = 0)  # 10 infected, 0 recovered initially
  
  # Define the control settings
  control <- control.net(type = "SIR", nsteps = 100, nsims = 1)
  
  # Run the simulation
  out <- netsim(est, param, init, control)
  
  # Extract the epidemic data from the epi object and return
  df <- as.data.frame(out$epi)
  df$time <- 1:nrow(df)  # Create a time variable
  df$Beta <- beta  # Add beta value
  df$Gamma <- gamma  # Add gamma value
  
  return(df)
}

# Define a single beta and gamma value
beta <- 0.3
gamma <- 0.1

# Run the simulation for this specific beta and gamma
sim_result <- run_sir_model(beta, gamma)

# Rename columns of sim_result to reflect compartments
colnames(sim_result)[1:4] <- c("pop_size", "susceptible", "infected", "recovered")

# Now, plot the results using ggplot
ggplot(sim_result, aes(x = time, y = infected)) +
  geom_line(linewidth = 1, color = "red") +  # Use linewidth and color for the infected individuals
  labs(title = paste("SIR Model Simulation (Beta =", beta, ", Gamma =", gamma, ")"),
       y = "Infected Individuals (i.num)",
       x = "Time Steps") +
  theme_minimal() +
  theme(legend.position = "none")

#10 interations of the above

# Load necessary libraries
library(EpiModel)
library(ggplot2)
library(dplyr)

# Function to run the SIR model for given beta and gamma
run_sir_model <- function(beta, gamma) {
  # Define the network with 1000 nodes
  nw <- network.initialize(1000, directed = FALSE)
  
  # Estimate the network model (simple network with constant degree)
  formation <- ~edges
  target.stats <- 500  # Average number of edges in the network
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation, target.stats, coef.diss)
  
  # Set the model parameters
  param <- param.net(inf.prob = beta, act.rate = 1, rec.rate = gamma)
  
  # Define the initial conditions
  init <- init.net(i.num = 10, r.num = 0)  # 10 infected, 0 recovered initially
  
  # Define the control settings
  control <- control.net(type = "SIR", nsteps = 250, nsims = 1)
  
  # Run the simulation
  out <- netsim(est, param, init, control)
  
  # Extract the epidemic data from the epi object and return
  df <- as.data.frame(out$epi)
  df$time <- 1:nrow(df)  # Create a time variable
  df$Beta <- beta  # Add beta value
  df$Gamma <- gamma  # Add gamma value
  
  return(df)
}

# Define a single beta and gamma value
beta <- 0.3
gamma <- 0.1

# Number of simulations
num_simulations <- 10  # Change to 10
sim_results <- vector("list", num_simulations)

# Run the simulation multiple times
for (i in 1:num_simulations) {
  sim_results[[i]] <- run_sir_model(beta, gamma)
}

# Combine all simulations into a single data frame
combined_results <- bind_rows(sim_results, .id = "simulation")

# Rename columns to reflect compartments
colnames(combined_results)[2:5] <- c("pop_size", "susceptible", "infected", "recovered")

# Calculate mean and standard deviation of infected individuals over time
mean_infected <- combined_results %>%
  group_by(time) %>%
  summarise(mean_infected = mean(infected, na.rm = TRUE),  # Use 'infected' instead of 'i.num'
            sd_infected = sd(infected, na.rm = TRUE))

# Plot the results with error bars and set x-axis limits
ggplot(mean_infected, aes(x = time, y = mean_infected)) +
  geom_line(linewidth = 1, color = "red") +  # Mean line
  geom_errorbar(aes(ymin = mean_infected - sd_infected, 
                    ymax = mean_infected + sd_infected),
                width = 0.1, color = "black") +  # Error bars
  labs(title = paste("SIR Model Simulation (Beta =", beta, ", Gamma =", gamma, ")"),
       y = "Mean Infected Individuals",
       x = "Time Steps") +
  #ylim(0, 50) +  # Set x-axis limits
  theme_minimal() +
  theme(legend.position = "none")

#10 interations of the above

# Load necessary libraries
library(EpiModel)
library(ggplot2)
library(dplyr)

# Function to run the SIR model for given beta and gamma
run_sir_model <- function(beta, gamma) {
  # Define the network with 1000 nodes
  nw <- network.initialize(1000, directed = FALSE)
  
  # Estimate the network model (simple network with constant degree)
  formation <- ~edges
  target.stats <- 500  # Average number of edges in the network
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation, target.stats, coef.diss)
  
  # Set the model parameters
  param <- param.net(inf.prob = beta, act.rate = 1, rec.rate = gamma)
  
  # Define the initial conditions
  init <- init.net(i.num = 10, r.num = 0)  # 10 infected, 0 recovered initially
  
  # Define the control settings
  control <- control.net(type = "SIR", nsteps = 250, nsims = 10)
  
  # Run the simulation
  out <- netsim(est, param, init, control)
  
  # Extract the epidemic data from the epi object and return
  df <- as.data.frame(out$epi)
  df$time <- 1:nrow(df)  # Create a time variable
  df$Beta <- beta  # Add beta value
  df$Gamma <- gamma  # Add gamma value
  
  return(df)
}

# Define a single beta and gamma value
beta <- 0.3
gamma <- 0.1


sim_results <- run_sir_model(beta, gamma)

# Combine all simulations into a single data frame
combined_results <- bind_rows(sim_results, .id = "simulation")

# Rename columns to reflect compartments
colnames(combined_results)[2:5] <- c("pop_size", "susceptible", "infected", "recovered")

# Calculate mean and standard deviation of infected individuals over time
mean_infected <- combined_results %>%
  group_by(time) %>%
  summarise(mean_infected = mean(infected, na.rm = TRUE),  # Use 'infected' instead of 'i.num'
            sd_infected = sd(infected, na.rm = TRUE))

# Plot the results with error bars and set x-axis limits
ggplot(mean_infected, aes(x = time, y = mean_infected)) +
  geom_line(linewidth = 1, color = "red") +  # Mean line
  geom_errorbar(aes(ymin = mean_infected - sd_infected, 
                    ymax = mean_infected + sd_infected),
                width = 0.1, color = "black") +  # Error bars
  labs(title = paste("SIR Model Simulation (Beta =", beta, ", Gamma =", gamma, ")"),
       y = "Mean Infected Individuals",
       x = "Time Steps") +
  #ylim(0, 50) +  # Set x-axis limits
  theme_minimal() +
  theme(legend.position = "none")

#### Agent based bodel with increased suseptability of those over 40

# Load necessary libraries
library(EpiModel)
library(ggplot2)

# Function to create an agent-based SIR model with age-based susceptibility
run_agent_based_model <- function(pop_size = 1000, beta_under_40 = 0.1, beta_over_40 = 0.3, gamma = 0.1, nsteps = 100) {
  
  # Step 1: Initialize the population (1000 individuals)
  nw <- network.initialize(n = pop_size, directed = FALSE)
  
  # Step 2: Assign random attributes (gender and age) to individuals
  # Gender: 0 for Female, 1 for Male (0.5 probability)
  gender <- rbinom(pop_size, 1, 0.5) 
  
  # Age: Normal distribution centered around 40 with SD = 15 (capped at 18:75)
  age <- round(rnorm(pop_size, mean = 40, sd = 15))
  age[age < 18] <- 18  # Set minimum age to 18
  age[age > 75] <- 75  # Set maximum age to 75
  
  # Step 3: Create the network model (simple undirected random network)
  formation <- ~edges
  target.stats <- 500  # Average number of edges in the network
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation, target.stats, coef.diss)
  
  # Step 4: Define the infection probability based on age (higher for those over 40)
  # This modifies beta dynamically based on each individual's age
  inf_prob <- ifelse(age > 40, beta_over_40, beta_under_40)
  
  # Step 5: Set the model parameters (using age-based infection probability)
  param <- param.net(inf.prob = inf_prob, act.rate = 1, rec.rate = gamma)
  
  # Step 6: Define the initial conditions (10 infected individuals, rest susceptible)
  init <- init.net(i.num = 10, r.num = 0)
  
  # Step 7: Define the control settings for the SIR model
  control <- control.net(type = "SIR", nsteps = nsteps, nsims = 1)
  
  # Step 8: Add gender and age attributes to the network
  nw <- set.vertex.attribute(nw, "gender", gender)
  nw <- set.vertex.attribute(nw, "age", age)
  
  # Step 9: Run the agent-based simulation
  out <- netsim(est, param, init, control)
  
  # Step 10: Extract the epidemic data and return
  df <- as.data.frame(out$epi)
  df$time <- 1:nrow(df)  # Create a time variable
  
  return(df)
}

# Define parameters
pop_size <- 1000
beta_under_40 <- 0.1  # Infection probability for individuals under 40
beta_over_40 <- 0.3   # Infection probability for individuals over 40
gamma <- 0.1          # Recovery rate
nsteps <- 200         # Number of timesteps

# Run the agent-based simulation
sim_result <- run_agent_based_model(pop_size, beta_under_40, beta_over_40, gamma, nsteps)

# Rename columns of sim_result to reflect compartments
colnames(sim_result)[1:4] <- c("pop_size", "susceptible", "infected", "recovered")

# Plot the results using ggplot
ggplot(sim_result, aes(x = time, y = infected)) +
  geom_line(linewidth = 1, color = "red") +  # Use linewidth and color for infected individuals
  labs(title = paste("Agent-Based SIR Model (Beta Under 40 =", beta_under_40, ", Beta Over 40 =", beta_over_40, ", Gamma =", gamma, ")"),
       y = "Infected Individuals (i.num)",
       x = "Time Steps") +
  theme_minimal()

####10x repeat of the above. 

# Load necessary libraries
library(EpiModel)
library(ggplot2)

# Function to create an agent-based SIR model with age-based susceptibility
run_agent_based_model <- function(pop_size = 1000, beta_under_40 = 0.1, beta_over_40 = 0.3, gamma = 0.1, nsteps = 100) {
  
  # Step 1: Initialize the population (1000 individuals)
  nw <- network.initialize(n = pop_size, directed = FALSE)
  
  # Step 2: Assign random attributes (gender and age) to individuals
  # Gender: 0 for Female, 1 for Male (0.5 probability)
  gender <- rbinom(pop_size, 1, 0.5) 
  
  # Age: Normal distribution centered around 40 with SD = 15 (capped at 18:75)
  age <- round(rnorm(pop_size, mean = 40, sd = 15))
  age[age < 18] <- 18  # Set minimum age to 18
  age[age > 75] <- 75  # Set maximum age to 75
  
  # Step 3: Create the network model (simple undirected random network)
  formation <- ~edges
  target.stats <- 500  # Average number of edges in the network
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation, target.stats, coef.diss)
  
  # Step 4: Define the infection probability based on age (higher for those over 40)
  # This modifies beta dynamically based on each individual's age
  inf_prob <- ifelse(age > 40, beta_over_40, beta_under_40)
  
  # Step 5: Set the model parameters (using age-based infection probability)
  param <- param.net(inf.prob = inf_prob, act.rate = 1, rec.rate = gamma)
  
  # Step 6: Define the initial conditions (10 infected individuals, rest susceptible)
  init <- init.net(i.num = 10, r.num = 0)
  
  # Step 7: Define the control settings for the SIR model
  control <- control.net(type = "SIR", nsteps = nsteps, nsims = 1)
  
  # Step 8: Add gender and age attributes to the network
  nw <- set.vertex.attribute(nw, "gender", gender)
  nw <- set.vertex.attribute(nw, "age", age)
  
  # Step 9: Run the agent-based simulation
  out <- netsim(est, param, init, control)
  
  # Step 10: Extract the epidemic data and return
  df <- as.data.frame(out$epi)
  df$time <- 1:nrow(df)  # Create a time variable
  
  return(df)
}

# Define parameters
pop_size <- 1000
beta_under_40 <- 0.1  # Infection probability for individuals under 40
beta_over_40 <- 0.3   # Infection probability for individuals over 40
gamma <- 0.1          # Recovery rate
nsteps <- 200         # Number of timesteps

# Number of simulations
num_simulations <- 10  # Change to 10
sim_results <- vector("list", num_simulations)

# Run the simulation multiple times
for (i in 1:num_simulations) {
  sim_results[[i]] <- run_agent_based_model(pop_size, beta_under_40, beta_over_40, gamma, nsteps)
}

# Combine all simulations into a single data frame
combined_results <- bind_rows(sim_results, .id = "simulation")

# Rename columns to reflect compartments
colnames(combined_results)[2:5] <- c("pop_size", "susceptible", "infected", "recovered")

# Calculate mean and standard deviation of infected individuals over time
mean_infected <- combined_results %>%
  group_by(time) %>%
  summarise(mean_infected = mean(infected, na.rm = TRUE),  # Use 'infected' instead of 'i.num'
            sd_infected = sd(infected, na.rm = TRUE))

# Plot the results with error bars and set x-axis limits
ggplot(mean_infected, aes(x = time, y = mean_infected)) +
  geom_line(linewidth = 1, color = "red") +  # Mean line
  geom_errorbar(aes(ymin = mean_infected - sd_infected, 
                    ymax = mean_infected + sd_infected),
                width = 0.1, color = "black") +  # Error bars
  labs(title = paste("SIR Model Simulation (Beta =", beta, ", Gamma =", gamma, ")"),
       y = "Mean Infected Individuals",
       x = "Time Steps") +
  #ylim(0, 50) +  # Set x-axis limits
  theme_minimal() +
  theme(legend.position = "none")

# Plot the simulation results using ggplot2 and facet_wrap
ggplot(combined_results, aes(x = time, y = infected)) +
  geom_line(linewidth = 1, color = "red") +  # Plot the infected individuals over time
  facet_wrap(~simulation, scales = "free_y") +  # Create a facet for each simulation
  labs(title = "SIR Model Simulation Results by Simulation",
       y = "Infected Individuals",
       x = "Time Steps") +
  theme_minimal() +
  theme(legend.position = "none")
