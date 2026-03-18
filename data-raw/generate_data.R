# Script to generate example dataset for fbardl package
# This script creates a simulated cointegrated time series with structural breaks

set.seed(12345)

n <- 150

# Generate I(1) independent variables
e1 <- rnorm(n, 0, 0.5)
e2 <- rnorm(n, 0, 0.5)

x1 <- cumsum(e1)
x2 <- cumsum(e2)

# Generate structural break using Fourier terms
ttrend <- 1:n
k <- 1.5  # Fourier frequency
fourier_component <- 0.5 * sin(2 * pi * k * ttrend / n) + 0.3 * cos(2 * pi * k * ttrend / n)

# Cointegrating relationship: y = 2 + 0.8*x1 - 0.5*x2 + fourier + error
# with error correction dynamics

# Long-run equilibrium
equilibrium <- 2 + 0.8 * x1 - 0.5 * x2 + fourier_component

# Generate y with error correction mechanism
alpha <- -0.3  # Speed of adjustment
y <- numeric(n)
y[1] <- equilibrium[1] + rnorm(1, 0, 0.3)

for (t in 2:n) {
  # Error correction term
  ecm <- y[t-1] - equilibrium[t-1]

  # Short-run dynamics
  dy <- alpha * ecm + 0.2 * (x1[t] - x1[t-1]) - 0.1 * (x2[t] - x2[t-1]) + rnorm(1, 0, 0.3)

  y[t] <- y[t-1] + dy
}

# Create data frame
fbardl_data <- data.frame(
  y = y,
  x1 = x1,
  x2 = x2
)

# Save as RDA file
usethis::use_data(fbardl_data, overwrite = TRUE)
