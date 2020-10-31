# Seperable Filters -------------------------------------------------------
library(raster)
library(spatialEco)
source("scripts/functions_for_seperable_filters.R")
  
# Matrix Multiplication ---------------------------------------------------

# Sobel Filter for edge detetction
sobel_v <- matrix(c(1,2,1), nrow = 3, ncol = 1)
sobel_h <- matrix(c(-1,0,1), nrow = 1, ncol = 3)

sobel_m <- sobel_v %*% sobel_h

# is the matrix seperable? if so the result is 1
qr(sobel_m)$rank


# test seperable filter on raster -----------------------------------------

# setup raster
m_r <- matrix(runif(10000),
              ncol = 100, nrow = 100)

r <- raster(m_r)
plot(r)

# weight matrix
weight_m <- gaussian.kernel(sigma = 10, n = 11)
outer_vecs <- get_outer_product_vectors(weight_m)

r_conv_mat <- focal(r, weight_m)
r_conv_sep_1 <- focal(r, outer_vecs$vertical)
r_conv_sep_final <- focal(r_conv_sep_1, outer_vecs$horizontal)

par(mfrow=c(1,2))
plot(r_conv_sep_final)
plot(r_conv_mat)

identical(round(r_conv_mat, 5), round(r_conv_sep_final, 5))



# Benchmarking ------------------------------------------------------------
library(microbenchmark)
library(tidyverse)


# setup raster
r <- raster(matrix(runif(1000000),
                   ncol = 1000, nrow = 1000))

# weight matrices
wm_11 <- gaussian.kernel(sigma = 11, n = 11)
wm_21 <- gaussian.kernel(sigma = 21, n = 21)
wm_31 <- gaussian.kernel(sigma = 31, n = 31)
wm_51 <- gaussian.kernel(sigma = 51, n = 51)
wm_81 <- gaussian.kernel(sigma = 81, n = 81)
wm_101 <- gaussian.kernel(sigma = 101, n = 101)
wm_501 <- gaussian.kernel(sigma = 501, n = 501)


ov_11 <- get_outer_product_vectors(wm_11)
ov_21 <- get_outer_product_vectors(wm_21)
ov_31 <- get_outer_product_vectors(wm_31)
ov_51 <- get_outer_product_vectors(wm_51)
ov_81 <- get_outer_product_vectors(wm_81)
ov_101 <- get_outer_product_vectors(wm_101)
ov_501 <- get_outer_product_vectors(wm_501)



focal_bench <- microbenchmark(
  focal(r, wm_11),
  focal(r, wm_21),
  focal(r, wm_31),
  focal(r, wm_51),
  focal(r, wm_81),
  focal(r, wm_101),
  focal_seperable(r, ov_11),
  focal_seperable(r, ov_21),
  focal_seperable(r, ov_31),
  focal_seperable(r, ov_51),
  focal_seperable(r, ov_81),
  focal_seperable(r, ov_101),
  times = 30L
)


focal_bench_plot <- 
  focal_bench %>%
  group_by(expr) %>%
  summarise(time = median(time)) %>%
  mutate(expr = str_replace(expr,
                            pattern = "focal_seperable",
                            replacement = ""),
         expr = str_replace(expr,
                                pattern = "focal",
                                replacement = ""),
         expr = str_replace(expr,
                                pattern = c("\\(r,"),
                                replacement = ""),
         expr = str_replace(expr,
                            pattern = c("\\)"),
                            replacement = ""),
         expr = str_trim(expr)
         ) %>%
  separate(col = expr,
           into = c("method", "matrix_dimension"),
           sep = "_") %>%
  mutate(method = case_when(
    method == "wm" ~ "normal",
    method == "ov" ~ "seperable filter"
  ),
  matrix_dimension = as.integer(matrix_dimension),
  time = time / 1000000000)

ggplot(focal_bench_plot, aes(x = matrix_dimension, y = time, group = method, color = method)) +
  geom_line() +
  theme_bw() +
  labs(x = "Dimension of Symmetrical Kernel",
       y = "Time in Seconds",
       title = "Performance of Weights Matrix vs. Seperable Filter",
       subtitle = "Gaussian Blur on a 1000x1000 Raster")
  


# Profiling ---------------------------------------------------------------

library(profvis)

profvis(
  focal(r, wm_101)
)

profvis(
  focal_seperable(r, ov_101)
)



# Kernels -----------------------------------------------------------------


