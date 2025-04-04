# Install required packages if not already installed
# install.packages("kernlab")
# install.packages("MASS")

# Set seed for reproducibility
set.seed(123)

# Generate synthetic data in 3 dimensions
n <- 1000  # Number of observations
d <- 5    # Number of dimension of the input space

# Generate input features from multivariate normal distribution
X <- MASS::mvrnorm(n, mu = rep(0, d), Sigma = 5* diag(d))
colnames(X) <- paste0("x", 1:d)

# Generate target values using a non-linear function with noise
true_f <- function(x) {
  return(sin(x[1]) + 0.5*cos(x[2]) + 0.3*x[3]^2) * cos(x[4+5])
}
y <- apply(X, 1, true_f) + rnorm(n, sd = 0.2)

# Create data frame
data <- data.frame(X, y = y)

# Fit Gaussian process model with explicit package reference
gp_model <- kernlab::gausspr(
  x = as.matrix(data[, 1:d]),
  y = data$y,
  kernel = "rbf",
  kpar = list(sigma = 0.5),
  var = 0.05
)

# Generate test data on a grid (for visualization purposes)
# We'll create a grid in the first two dimensions and fix the third dimension
grid_size <- 20
x1_seq <- seq(min(X[, 1]), max(X[, 1]), length.out = grid_size)
x2_seq <- seq(min(X[, 2]), max(X[, 2]), length.out = grid_size)
x3_fixed <- median(X[, 3])  # Fix the third dimension at its median value
x4_fixed <- median(X[, 4])  # Fix the third dimension at its median value
x5_fixed <- median(X[, 5])  # Fix the third dimension at its median value

# Create grid of test points
test_grid <- expand.grid(x1 = x1_seq, x2 = x2_seq, x3 = x3_fixed, x4 = x4_fixed, x5 = x5_fixed)

# Make predictions on test grid
predictions <- kernlab::predict(gp_model, as.matrix(test_grid))

# Visualize the results (2D slice of 3D function)
# Convert predictions to matrix for contour plot
z_matrix <- matrix(predictions, nrow = grid_size, ncol = grid_size)

# Plot the true function values on the same grid for comparison
true_values <- apply(as.matrix(test_grid), 1, true_f)
true_matrix <- matrix(true_values, nrow = grid_size, ncol = grid_size)

# Create plots
par(mfrow = c(1, 2))
contour(x1_seq, x2_seq, z_matrix, main = "GP Prediction",
        xlab = "x1", ylab = "x2", nlevels = 15)
points(X[, 1], X[, 2], pch = 19, cex = 0.5, col = "darkgray")

contour(x1_seq, x2_seq, true_matrix, main = "True Function",
        xlab = "x1", ylab = "x2", nlevels = 15)
points(X[, 1], X[, 2], pch = 19, cex = 0.5, col = "darkgray")
