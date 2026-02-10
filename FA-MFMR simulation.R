# ===================================================================
# FA-MFMR 对比仿真代码 
# ===================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(fda)
  library(fdapace)
  library(glmnet)
  library(AER)
  library(ncvreg)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(TVMR)
  library(parallel)   
})

# ===================================================================
# -----------参数配置-----------
# ===================================================================

config <- list(
  # 基本参数
  n = 5000,              # 样本量
  G = 5,                # 组数 (扩展到20个)
  m = 6,                 # 基函数维数
  K_true = 3,            # 真实因子数
  J = 100,               # 工具变量个数
  
  # 时间网格
  t_range = c(0, 50),
  n_grid = 51,           # 参考网格点数
  
  # 随机观测时间参数
  min_obs = 5,          # 每个个体最少观测次数
  max_obs = 15,         # 每个个体最多观测次数
  obs_prob_center = 0.7, # 中心区域观测概率
  obs_prob_edge = 0.3,   # 边缘区域观测概率
  
  # 噪声参数
  sigma_eps = 0.1,       # 观测噪声
  sigma_y = 0.05,        # 响应变量噪声
  
  # 评估参数
  key_times = c(10, 20, 30, 40),  # 关键时间点
  
  # 种子
  seed = 2025
)

# 计算派生参数
config$p <- config$G * config$m
config$t_grid <- seq(config$t_range[1], config$t_range[2], length.out = config$n_grid)
config$dt <- mean(diff(config$t_grid))

cat("=== FA-MFMR vs MFMR 对比仿真参数 ===\n")
cat(sprintf("样本量: %d, 组数: %d, 基函数维数: %d\n", config$n, config$G, config$m))
cat(sprintf("真实因子数: %d, 工具变量数: %d\n", config$K_true, config$J))
cat(sprintf("时间范围: [%.1f, %.1f], 参考网格点数: %d\n", 
            config$t_range[1], config$t_range[2], config$n_grid))
cat("Beta函数类型: 前4个非零(二次、线性、三角、高斯) + 16个零函数\n")
cat("对比方法: FA-MFMR (因子结构+IV) vs MFMR\n")
cat("==========================\n\n")

# ===================================================================
# -----------辅助信息显示函数-----------
# ===================================================================

#' 生成Fourier基函数
# generate_fourier_basis <- function(t_grid, m) {
#   # 归一化时间到[0,1]
#   t_normalized <- (t_grid - min(t_grid)) / (max(t_grid) - min(t_grid))
# 
#   # 创建Fourier基
#   fourier_basis <- create.fourier.basis(c(0, 1), nbasis = m, dropind = 1)
#   basis_eval <- eval.basis(t_normalized, fourier_basis)
# 
#   basis_eval <- scale(basis_eval, center = FALSE,
#                       scale = sqrt(colSums(basis_eval^2)))
# 
#   return(basis_eval)
# }

#' 生成B-spline基函数
generate_bspline_basis <- function(t_grid, m) {
  t_min <- min(t_grid)
  t_max <- max(t_grid)
  
  # 创建B-spline基（3次样条）
  bspline_basis <- create.bspline.basis(c(t_min, t_max), nbasis = m, norder = 4)
  basis_eval <- eval.basis(t_grid, bspline_basis)
  
  basis_eval <- scale(basis_eval, center = FALSE,
                      scale = sqrt(colSums(basis_eval^2)))
  return(basis_eval)
}

#' 生成固定的真实回归系数函数 (20个函数: 前4个非零 + 16个零函数)
generate_true_beta_functions <- function(t_grid) {
  n_grid <- length(t_grid)
  G <- 20  # 扩展到20个函数
  
  # 归一化时间到[0,1]
  t_min <- min(t_grid)
  t_max <- max(t_grid)
  t_normalized <- (t_grid - t_min) / (t_max - t_min)
  
  # 初始化beta函数矩阵
  beta_true_mat <- matrix(0, nrow = n_grid, ncol = G)
  
  # Beta1: 线性函数 beta2(t) = t - 0.5
  beta_true_mat[, 1] <- 2* (t_normalized - 0.5)
  
  # Beta2: 二次函数 beta1(t) = 2*(t-0.5)² - c1_theoretical
  c1_theoretical <- 1/3
  beta_true_mat[,2] <- 4 * (t_normalized - 0.5)^2 - c1_theoretical
  
  # Beta3: 三角函数 beta3(t) = 0.3*sin(2π*t)
  beta_true_mat[, 3] <- sin(2 * pi * t_normalized)
  
  # Beta4: 高斯函数 beta4(t) = exp(-(t-0.5)²/(0.2²)) - c4_theoretical
  mu <- 0.5
  sigma <- 0.2
  t_integration <- seq(0, 1, length.out = 10001)
  c4_theoretical <- mean(exp(-(t_integration - mu)^2 / (sigma^2)))
  beta_true_mat[, 4] <- exp(-(t_normalized - mu)^2 / (sigma^2)) - c4_theoretical 
  
  # Beta5-Beta20: 零函数 (已经在初始化时设为0)
  
  return(beta_true_mat)
}

#' 生成工具变量和相关矩阵 
generate_instruments <- function(n, J, K_true, p) {
  
  # 生成工具变量（遗传变异）
  allele_freqs <- runif(J, 0.2, 0.4)
  Z <- matrix(0, nrow = n, ncol = J)
  for (j in 1:J) {
    Z[, j] <- rbinom(n, size = 2, prob = allele_freqs[j])
  }
  Z <- scale(Z)
  
  # 生成因子矩阵
  Pi_F <- matrix(runif(J * K_true, -2, 2), nrow = J)
  F_mat <- Z %*% Pi_F + matrix(rnorm(n * K_true, sd = 0.1), n, K_true)
  
  F_mat <- scale(F_mat)
  
  # 生成特质矩阵
  Pi_U <- matrix(runif(J * p, -1.5, 1.5), nrow = J)
  U_mat <- Z %*% Pi_U + matrix(rnorm(n * p, sd = 0.2), n, p)
  U_mat <- scale(U_mat)
  
  # 生成系数矩阵
  B_mat <- matrix(rnorm(p * K_true, sd = 0.5), nrow = p, ncol = K_true)
  A_mat <- F_mat %*% t(B_mat) + U_mat
  
  return(list(
    Z = Z,
    A_mat = A_mat,
    F_mat = F_mat,
    U_mat = U_mat,
    Pi_F = Pi_F,
    Pi_U = Pi_U
  ))
}

#' 为每个个体生成随机观测时间
generate_random_obs_times <- function(n, t_range, min_obs, max_obs, center_prob, edge_prob) {
  obs_times_list <- vector("list", n)
  
  for (i in 1:n) {
    n_obs_i <- sample(min_obs:max_obs, 1)
    candidate_times <- seq(t_range[1], t_range[2], length.out = 200)
    
    t_mid <- mean(t_range)
    t_span <- diff(t_range)
    center_region <- abs(candidate_times - t_mid) <= 0.3 * t_span
    
    obs_probs <- ifelse(center_region, center_prob, edge_prob)
    obs_probs <- obs_probs / sum(obs_probs)
    
    obs_indices <- sample(length(candidate_times), n_obs_i, 
                          prob = obs_probs, replace = FALSE)
    obs_times_i <- sort(candidate_times[obs_indices])
    
    if (runif(1) < 0.3 && !any(abs(obs_times_i - t_range[1]) < 2)) {
      obs_times_i[1] <- t_range[1] + runif(1, 0, 3)
    }
    if (runif(1) < 0.3 && !any(abs(obs_times_i - t_range[2]) < 2)) {
      obs_times_i[length(obs_times_i)] <- t_range[2] - runif(1, 0, 3)
    }
    
    obs_times_list[[i]] <- sort(obs_times_i)
  }
  
  return(obs_times_list)
}

#' 生成完整函数轨迹 
generate_complete_trajectories <- function(config, A_mat, basis_eval) {
  X_complete_array <- array(0, dim = c(config$n, config$G, length(config$t_grid)))
  bm_t50_values <- numeric(config$n)
  
  
  # === 新增：每个个体都有自己的混淆路径 ===
  common_bm_indiv <- matrix(0, nrow = config$n, ncol = length(config$t_grid))
  for (i in 1:config$n) {
    common_bm_i <- cumsum(rnorm(length(config$t_grid), mean = 0, sd = 1))
    common_bm_i <- (common_bm_i - mean(common_bm_i)) / sd(common_bm_i)
    common_bm_indiv[i, ] <- common_bm_i
  }
  
  # 归一化时间
  t_normalized <- (config$t_grid - min(config$t_grid)) / 
    (max(config$t_grid) - min(config$t_grid))
  
  # 找到t=50对应的索引
  idx_t50 <- which.min(abs(config$t_grid - 50))
  n_time <- length(t_normalized)
  
  for (i in 1:config$n) {
    
    bm_t50_sum <- 0  # 累加该个体所有组的BM在t=50的值
    
    for (g in 1:config$G) {
      idx_start <- (g - 1) * config$m + 1
      idx_end <- g * config$m
      a_g <- A_mat[i, idx_start:idx_end]
      
      # 生成基础轨迹
      base_traj <- as.numeric(basis_eval %*% a_g)
      
      # 生成噪声BM路径（代替原来的独立噪声）
      bm_noise <- cumsum(rnorm(n_time, mean = 0, sd = 1))
      bm_noise <- (bm_noise - mean(bm_noise)) / sd(bm_noise)
      noise_traj <- sqrt(config$sigma_eps) * bm_noise
      
      # # 为每个(i,g)对生成独立的 Brownian motion 路径
      # bm_path <- cumsum(rnorm(n_time, mean = 0, sd = 1))
      # bm_path <- (bm_path - mean(bm_path)) / sd(bm_path)
      
      # === 原本独立的 Brownian path ===
      indiv_bm <- cumsum(rnorm(n_time))
      indiv_bm <- (indiv_bm - mean(indiv_bm)) / sd(indiv_bm)
      
      # === 新增：混入全局成分，产生个体间相关结构 ===
      bm_path <-  0.6 * common_bm_indiv[i, ] + 0.4 * indiv_bm
      bm_path <- (bm_path - mean(bm_path)) / sd(bm_path)
      
      # 累加该组在t=50的BM值
      bm_t50_sum <- bm_t50_sum + bm_path[idx_t50]
      
      confounder_effect <-  bm_path 
      
      # confounder_effect <- Confounder[i] * sin(pi * t_normalized)
      
      X_complete_array[i, g, ] <- base_traj +  noise_traj + confounder_effect 
    }
    
    # 计算该个体所有组在t=50的BM平均值
    # bm_t50_values[i] <- bm_t50_sum / config$G
    bm_t50_values[i] <- bm_t50_sum
  }
  
  return(list(
    X_complete_array = X_complete_array,
    bm_t50_values = bm_t50_values,  # 新增返回值
    common_bm_indiv = common_bm_indiv 
  ))
}

#' 从完整轨迹中抽取观测数据
extract_observed_data <- function(X_complete_array, obs_times_list, config) {
  Ly_list <- vector("list", config$G)
  Lt_list <- vector("list", config$G)
  
  for (g in 1:config$G) {
    Ly_list[[g]] <- vector("list", config$n)
    Lt_list[[g]] <- vector("list", config$n)
  }
  
  for (i in 1:config$n) {
    obs_times_i <- obs_times_list[[i]]
    
    for (g in 1:config$G) {
      complete_traj <- X_complete_array[i, g, ]
      observed_values <- approx(x = config$t_grid, 
                                y = complete_traj, 
                                xout = obs_times_i, 
                                rule = 2)$y
      
      Ly_list[[g]][[i]] <- observed_values
      Lt_list[[g]][[i]] <- obs_times_i
    }
  }
  
  return(list(Ly_list = Ly_list, Lt_list = Lt_list))
}

#' 生成响应变量 (使用真实beta函数)
generate_response_from_complete_trajectories <- function(config, X_complete_array, beta_true_mat, bm_t50_values, common_bm_indiv, A_mat) {
  Y <- numeric(config$n)
  
  # 找到t=50对应的索引
  idx_t50 <- which.min(abs(config$t_grid - 50))
  
  #还原尺度
  t_span <- max(config$t_grid) - min(config$t_grid)
  
  for (i in 1:config$n) {
    val <- 0
    for (g in 1:config$G) {
      X_ig_complete <- X_complete_array[i, g, ]
      beta_g_true <- beta_true_mat[, g]
      
      integrand <- X_ig_complete * beta_g_true
      val <- val + sum(integrand) * config$dt 
    }
    
    # 使用该个体在t=50处的BM值作为混淆效应
    confounder_direct_effect <- bm_t50_values[i]
    
    # 为Y生成独立的BM路径作为噪声
    n_time <- length(config$t_grid)
    bm_Y_noise <- cumsum(rnorm(n_time, mean = 0, sd = 1))
    bm_Y_noise <- (bm_Y_noise - mean(bm_Y_noise)) / sd(bm_Y_noise)
    noise_effect <- config$sigma_y * bm_Y_noise[idx_t50]
    
    Y[i] <- val + confounder_direct_effect + noise_effect 
  }
  
  return(Y)
}

# ===================================================================
# -----------FPCA计算函数-----------
# ===================================================================

#' 计算所有组的FPCA
compute_fpca_for_all_groups <- function(func_data, config) {
  tryCatch({
    fpca_list <- vector("list", config$G)
    A_est_list <- list()
    m_vec <- integer(config$G)
    
    cat("  执行FPCA计算...")
    for (g in 1:config$G) {
      fpca_g <- FPCA(Ly = func_data$Ly_list[[g]], 
                     Lt = func_data$Lt_list[[g]], 
                     list(dataType = 'Sparse', 
                          error = TRUE,
                          methodSelectK = "FVE",
                          FVEthreshold = 0.95,
                          maxK = 15,
                          userBwMu = 1.5,
                          userBwCov = 2.0,
                          nRegGrid = 101,
                          verbose = FALSE))
      
      fpca_list[[g]] <- fpca_g
      
      FVE <- fpca_g$cumFVE
      K_g <- which(FVE >= 0.95)[1]
      if (is.na(K_g)) K_g <- min(10, length(FVE))
      m_vec[g] <- K_g
      A_est_list[[g]] <- fpca_g$xiEst[, 1:K_g, drop = FALSE]
    }
    cat("完成\n")
    
    # 构造估计的系数矩阵
    m_total <- sum(m_vec)
    A_est <- matrix(0, nrow = config$n, ncol = m_total)
    col_start <- 1
    for (g in 1:config$G) {
      K_g <- m_vec[g]
      col_end <- col_start + K_g - 1
      A_est[, col_start:col_end] <- A_est_list[[g]]
      col_start <- col_end + 1
    }
    
    return(list(
      success = TRUE,
      fpca_list = fpca_list,
      A_est_list = A_est_list,
      m_vec = m_vec,
      A_est = A_est
    ))
    
  }, error = function(e) {
    cat("失败:", e$message, "\n")
    return(list(success = FALSE, error = e$message))
  })
}

# ===================================================================
# -----------Post-Lasso去偏估计核心函数-----------
# ===================================================================

#' Projection + SCAD估计核心函数
projection_scad_estimation <- function(Y, F_iv, U_iv, cv_folds = 10) {
  
  # 第一步：投影方法 - 把因子效应先去掉
  cat("    执行投影方法...")
  
  # Y对F_iv回归，得到残差
  if (ncol(F_iv) > 0) {
    fit_Y_F <- lm(Y ~ F_iv - 1)  # 不包含截距，因为数据已中心化
    Y_residual <- residuals(fit_Y_F)
    gamma_hat <- coef(fit_Y_F)
  } else {
    Y_residual <- Y
    gamma_hat <- numeric(0)
  }
  
  # U_iv对F_iv回归，得到残差矩阵
  if (ncol(F_iv) > 0) {
    U_iv_residual <- matrix(0, nrow = nrow(U_iv), ncol = ncol(U_iv))
    for (j in 1:ncol(U_iv)) {
      fit_U_F <- lm(U_iv[, j] ~ F_iv - 1)
      U_iv_residual[, j] <- residuals(fit_U_F)
    }
  } else {
    U_iv_residual <- U_iv
  }
  
  cat("完成\n")
  
  # 第二步：对残差执行SCAD回归
  cat("    执行SCAD回归...")
  
  # 使用ncvreg包进行SCAD回归，交叉验证选择最优lambda
  cv_fit <- cv.ncvreg(U_iv_residual, Y_residual, penalty = "SCAD", 
                      nfolds = cv_folds, seed = 123)
  lambda_opt <- cv_fit$lambda.min
  
  # 直接从cv_fit中提取对应lambda.min的系数，避免重复拟合
  lambda_index <- which.min(abs(cv_fit$lambda - lambda_opt))
  beta_scad <- as.vector(cv_fit$fit$beta[-1, lambda_index])  # 去掉截距
  
  cat("完成\n")
  
  # 返回结果
  return(list(
    gamma_hat = gamma_hat,          # 因子系数
    beta_scad = beta_scad,          # 特质系数（SCAD估计）
    lambda_opt = lambda_opt,        # 最优正则化参数
    n_selected = sum(abs(beta_scad) > 1e-8)  # 选中的变量数量
  ))
}

#' 使用ratio method选择因子个数
select_factors_ratio_method <- function(X, K_max) {
  # 计算协方差矩阵的特征值
  cov_X <- cov(X)
  eigenvals <- eigen(cov_X, only.values = TRUE)$values
  eigenvals <- sort(eigenvals, decreasing = TRUE)
  
  # 确保特征值为正
  eigenvals <- pmax(eigenvals, 1e-10)
  
  # 计算ratio准则
  ratios <- numeric(K_max)
  for (k in 1:K_max) {
    if (k < length(eigenvals)) {
      ratios[k] <- eigenvals[k] / eigenvals[k + 1]
    } else {
      ratios[k] <- eigenvals[k] / 1e-10
    }
  }
  
  # 选择使ratio最大的k
  K_selected <- which.max(ratios)
  
  return(K_selected)
}


#' #' 使用信息准则选择因子个数
#' select_factors_ic_method <- function(X, K_max) {
#'   n <- nrow(X)
#'   p <- ncol(X)
#' 
#'   cov_X <- cov(X)
#'   eigenvals <- eigen(cov_X, only.values = TRUE)$values
#'   eigenvals <- sort(eigenvals, decreasing = TRUE)
#'   eigenvals <- pmax(eigenvals, 1e-10)
#' 
#'   IC_values <- numeric(K_max)
#' 
#'   for (k in 1:K_max) {
#'     # 计算残差方差
#'     if (k < length(eigenvals)) {
#'       V_k <- mean(eigenvals[(k+1):length(eigenvals)])
#'     } else {
#'       V_k <- 1e-10
#'     }
#' 
#'     # IC准则 (参考Bai & Ng 2002)
#'     penalty <- k * ((n + p)/(n * p)) * log(min(n, p))
#'     IC_values[k] <- log(V_k) + penalty
#'   }
#' 
#'   K_selected <- which.min(IC_values)
#' 
#'   cat(sprintf("  IC values (前10个): %s\n",
#'               paste(sprintf("%.4f", IC_values[1:min(10, K_max)]), collapse=", ")))
#' 
#'   return(K_selected)
#' }


#' 因子模型估计 
perform_factor_model_estimation <- function(fpca_results, config) {
  # 因子模型估计
  A_centered <- scale(fpca_results$A_est, center = TRUE, scale = FALSE)
  K_max_search <- min(10, floor(ncol(fpca_results$A_est)/2))
  
  # 选择因子个数
  K_selected <- select_factors_ratio_method(A_centered, K_max_search)
  
  # SVD分解和因子估计
  svd_A <- svd(A_centered)
  F_hat <- sqrt(config$n) * svd_A$u[, 1:K_selected, drop = FALSE]
  
  # 计算载荷矩阵
  D_k <- diag(svd_A$d[1:K_selected], nrow = K_selected)
  V_k <- svd_A$v[, 1:K_selected, drop = FALSE]
  B_hat <- sqrt(config$n) * V_k %*% D_k / config$n
  
  # 计算特质成分
  low_rank_part <- F_hat %*% t(B_hat)
  U_hat <- A_centered - low_rank_part
  
  return(list(
    K_selected = K_selected,
    F_hat = F_hat,
    U_hat = U_hat
  ))
}

# ===================================================================
# -----------主要算法: FA-MFMR vs MFMR vs 2SFIR vs MPCMR-----------
# ===================================================================

#' MPCMR方法实现 (改进为多暴露联合MPCMR)
run_mpcmr <- function(fpca_results, Y, instruments, config, 
                      calculate_ci = FALSE, nLM = 20, Parallel = TRUE,
                      gmm_file_path = "gmm_lm_onesample.R") {
  tryCatch({
    cat("  运行多暴露MPCMR (直接拼接系数后GMM估计)...\n")
    
    # ===== 步骤0: 加载原作者的GMM核心函数 =====
    cat("    步骤0: 加载GMM核心函数...", gmm_file_path, "\n")
    if (!file.exists(gmm_file_path)) {
      stop(sprintf("找不到文件: %s\n请指定正确的gmm_lm_onesample.R文件路径", gmm_file_path))
    }
    source(gmm_file_path)
    cat("    加载完成\n")
    
    # ===== 步骤1: 直接拼接所有暴露的FPCA系数矩阵 =====
    cat("    步骤1: 拼接所有暴露的FPCA系数矩阵...")
    
    Xi_combined <- matrix(0, nrow = config$n, ncol = sum(fpca_results$m_vec))
    col_idx <- 1
    
    for (g in 1:config$G) {
      K_g <- fpca_results$m_vec[g]
      fpca_g <- fpca_results$fpca_list[[g]]
      Xi_combined[, col_idx:(col_idx + K_g - 1)] <- fpca_g$xiEst[, 1:K_g]
      col_idx <- col_idx + K_g
    }
    
    cat("完成\n")
    
    # ===== 步骤2: 调用gmm_lm_onesample进行联合GMM估计 =====
    cat("    步骤2: 执行GMM估计 (联合处理所有暴露)...")
    
    Z <- instruments$Z
    
    # 核心：直接对拼接的系数矩阵进行GMM估计
    gmm_res <- gmm_lm_onesample(
      X = Xi_combined,
      Y = Y,
      Z = Z
    )
    
    cat("完成\n")
    
    # ===== 步骤3: 提取系数并分配回各暴露 =====
    cat("    步骤3: 重构各暴露的系数函数...\n")
    
    coeff_combined <- gmm_res$gmm_est
    n_components <- sum(fpca_results$m_vec)
    
    # 检查系数长度
    if (length(coeff_combined) != n_components) {
      cat(sprintf("\n警告: 系数长度(%d)与主成分数(%d)不匹配\n", 
                  length(coeff_combined), n_components))
      if (length(coeff_combined) < n_components) {
        coeff_combined <- c(coeff_combined, rep(0, n_components - length(coeff_combined)))
      } else {
        coeff_combined <- coeff_combined[1:n_components]
      }
    }
    
    # 分解系数回各暴露
    beta_est_functions <- list()
    col_idx <- 1
    
    for (g in 1:config$G) {
      K_g <- fpca_results$m_vec[g]
      fpca_g <- fpca_results$fpca_list[[g]]
      
      # 提取该暴露对应的系数
      coeff_g <- coeff_combined[col_idx:(col_idx + K_g - 1)]
      
      # 通过原FPCA基函数重构系数函数
      time_points <- fpca_g$workGrid
      phi_g <- fpca_g$phi[, 1:K_g, drop = FALSE]
      beta_values <- as.vector(phi_g %*% coeff_g)
      
      # 计算Wald类型的置信区间 (基于GMM方差估计)
      pointwise_var <- diag(phi_g %*% gmm_res$variance_matrix[(col_idx):(col_idx+K_g-1), 
                                                              (col_idx):(col_idx+K_g-1)] %*% t(phi_g))
      beta_low_wald <- beta_values - 1.96 * sqrt(pointwise_var)
      beta_up_wald <- beta_values + 1.96 * sqrt(pointwise_var)
      
      beta_est_functions[[g]] <- list(
        time = time_points,
        beta = beta_values,
        coeff = coeff_g,
        K_g = K_g,
        # Wald类型置信区间
        ci_low_wald = beta_low_wald,
        ci_up_wald = beta_up_wald,
        ci_var = pointwise_var,
        # LM置信区间（如果计算）
        ci_low_lm = NA,
        ci_up_lm = NA
      )
      
      col_idx <- col_idx + K_g
    }
    
    cat("完成\n")
    
    # ===== 步骤4: 计算LM置信区间（可选，通常耗时） =====
    if (calculate_ci) {
      cat("    步骤4: 计算LM置信区间...\n")
      
      # 准备候选点
      nn <- nLM
      K_total <- n_components
      
      beta_candidates <- c()
      for (k in 1:K_total) {
        seqs <- seq(gmm_res$gmm_est[k] - 4*gmm_res$gmm_se[k], 
                    gmm_res$gmm_est[k] + 4*gmm_res$gmm_se[k], 
                    length = nn)
        beta_candidates <- cbind(beta_candidates, 
                                 rep(seqs, each = nn^(K_total - k)))
      }
      
      # LM检验函数
      vector_to_LM <- function(vector) {
        LMres <- getLM(X = Xi_combined, Y = Y, Z = Z, beta0 = vector)
        return(LMres$lm_pval > 0.05)  # TRUE表示在CI内部
      }
      
      # 并行或顺序计算
      if (!Parallel) {
        LMres_vector <- apply(beta_candidates, 1, vector_to_LM)
      } else {
        # 并行计算
        Cores_used <- detectCores() - 1
        cl <- makeCluster(Cores_used)
        clusterExport(cl = cl, 
                      varlist = c('Xi_combined', 'Y', 'Z', 'getLM', 'vector_to_LM'),
                      envir = environment())
        LMres_vector <- parApply(cl, beta_candidates, 1, vector_to_LM)
        stopCluster(cl)
      }
      
      # 提取CI内部的点
      LMres_used <- beta_candidates[LMres_vector, ]
      if (sum(LMres_vector) == 1) {
        LMres_used <- t(as.matrix(LMres_used))
      }
      
      # 为各暴露重构LM置信区间
      col_idx <- 1
      for (g in 1:config$G) {
        K_g <- fpca_results$m_vec[g]
        fpca_g <- fpca_results$fpca_list[[g]]
        
        # 提取该暴露对应的系数范围
        coeff_range <- LMres_used[, col_idx:(col_idx + K_g - 1), drop = FALSE]
        
        # 重构时间轴上的置信区间
        phi_g <- fpca_g$phi[, 1:K_g, drop = FALSE]
        time_points <- fpca_g$workGrid
        
        pointwise_LM_range <- phi_g %*% t(coeff_range)
        ci_low_lm <- apply(pointwise_LM_range, 1, min)
        ci_up_lm <- apply(pointwise_LM_range, 1, max)
        
        beta_est_functions[[g]]$ci_low_lm <- ci_low_lm
        beta_est_functions[[g]]$ci_up_lm <- ci_up_lm
        
        col_idx <- col_idx + K_g
      }
      
      cat("完成\n")
    }
    
    # ===== 步骤5: 计算IV有效性检验 =====
    cat("    步骤5: IV有效性检验...\n")
    
    IV_validity <- c(
      Q_statistic = gmm_res$Q_stat,
      df = ncol(Z) - n_components,
      p_value = gmm_res$Q_pval
    )
    
    cat("完成\n")
    
    cat("  多暴露MPCMR估计完成！\n")
    
    # ===== 返回结果 =====
    return(list(
      success = TRUE,
      beta_est_functions = beta_est_functions,
      coeff_combined = coeff_combined,
      Xi_combined = Xi_combined,
      gmm_res = gmm_res,
      IV_validity = IV_validity,
      method_type = "MPCMR_multi",
      strategy = "direct_concatenation_gmm",
      n_components_total = n_components,
      n_instruments = ncol(Z),
      # 返回核心的置信区间计算函数，以便后续使用
      get_coefficient_ci = function(g, ci_type = "wald") {
        if (g < 1 || g > config$G) {
          stop(sprintf("暴露组号应在1-%d之间", config$G))
        }
        result <- list(
          time = beta_est_functions[[g]]$time,
          beta = beta_est_functions[[g]]$beta,
          ci_low = if(ci_type == "wald") beta_est_functions[[g]]$ci_low_wald 
          else beta_est_functions[[g]]$ci_low_lm,
          ci_up = if(ci_type == "wald") beta_est_functions[[g]]$ci_up_wald 
          else beta_est_functions[[g]]$ci_up_lm
        )
        return(result)
      },
      note = "所有暴露系数直接拼接后通过gmm_lm_onesample进行联合GMM估计。置信区间包括Wald型和LM型（若calculate_ci=TRUE）"
    ))
    
  }, error = function(e) {
    cat("失败:", e$message, "\n")
    cat("详细错误信息:\n")
    print(e)
    return(list(success = FALSE, error = e$message))
  })
}


#' 多元版 PACE + 2SFRI (两阶段都是功能型回归)
run_pace_2sfri <- function(fpca_results, Y, instruments, config) {
  tryCatch({
    cat("  运行PACE+2SFRI方法...\n")
    
    # ===== 步骤1: PACE重建各暴露曲线 =====
    cat("    步骤1: PACE重建各暴露曲线...")
    # 已经在fpca_results中完成了PACE，直接使用结果
    # x_ik(t) = μ_k(t) + Σ ξ_ikm φ_km(t)
    cat("完成\n")
    
    # ===== 步骤2: 定义多元去均值曲线 =====
    cat("    步骤2: 计算去均值曲线...")
    # D_ik(t) = x_ik(t) - μ_k(t)
    D_list <- vector("list", config$G)
    
    for (g in 1:config$G) {
      fpca_g <- fpca_results$fpca_list[[g]]
      mu_g <- fpca_g$mu  # 均值函数
      work_grid <- fpca_g$workGrid
      
      # 对每个个体，计算去均值后的曲线
      D_ig <- matrix(0, nrow = config$n, ncol = length(work_grid))
      for (i in 1:config$n) {
        # 重建曲线: x_ik(t) = μ_k(t) + Σ ξ_ikm φ_km(t)
        K_g <- fpca_results$m_vec[g]
        xi_ig <- fpca_g$xiEst[i, 1:K_g]
        phi_g <- fpca_g$phi[, 1:K_g, drop = FALSE]
        x_ig_reconstructed <- mu_g + as.vector(phi_g %*% xi_ig)
        
        # 去均值
        D_ig[i, ] <- x_ig_reconstructed - mu_g
      }
      
      D_list[[g]] <- list(
        values = D_ig,
        grid = work_grid
      )
    }
    cat("完成\n")
    
    # ===== 步骤3: 第一阶段（多元功能线性回归）=====
    cat("    步骤3: 第一阶段功能型回归 (IV -> D_i(t))...\n")
    # D_i(t) = β_0(t) + B_1(t)G_i + B_2(t)Z_i + r_i(t)
    # 这里 B_1(t) 和 B_2(t) 都是时变系数
    
    Z <- instruments$Z
    q <- ncol(Z)
    
    # 对每个暴露k，进行功能型线性回归
    residuals_list <- vector("list", config$G)
    B2_estimates <- vector("list", config$G)
    
    for (g in 1:config$G) {
      work_grid <- D_list[[g]]$grid
      n_time <- length(work_grid)
      D_ig <- D_list[[g]]$values
      
      # 存储残差和时变系数估计
      r_ig <- matrix(0, nrow = config$n, ncol = n_time)
      B2_gt <- matrix(0, nrow = n_time, ncol = q)  # Z_i的系数
      
      # 对每个时间点t进行OLS回归，得到时变系数
      for (t_idx in 1:n_time) {
        D_it <- D_ig[, t_idx]  # 所有个体在时间t的值
        
        # 回归: D_i(t) ~ Z_i
        fit <- lm(D_it ~ Z)
        r_ig[, t_idx] <- residuals(fit)
        B2_gt[t_idx, ] <- coef(fit)[-1]  # 去掉截距，得到B_2(t)
      }
      
      residuals_list[[g]] <- list(
        values = r_ig,
        grid = work_grid
      )
      B2_estimates[[g]] <- B2_gt
      
      cat(sprintf("      暴露%d: 第一阶段完成\n", g))
    }
    cat("    步骤3完成\n")
    
    # ===== 步骤4: 第二阶段（功能线性回归，连续结局）=====
    cat("    步骤4: 第二阶段功能型回归 (D_i(t), r_i(t) -> y_i)...\n")
    # y_i = γ_0 + Σ_k ∫ γ_1k(t) D_ik(t) dt + Σ_k ∫ γ_3k(t) r_ik(t) dt + ε_i
    # (不包含Z项)
    
    # 使用FPC展开来表示时变系数 γ_1k(t) 和 γ_3k(t)
    # 即: γ_1k(t) = Σ_m α_1km φ_km(t)
    
    # 构造设计矩阵
    # 第一部分: D_ik(t) 的FPC系数
    D_fpc_scores <- matrix(0, nrow = config$n, ncol = sum(fpca_results$m_vec))
    col_idx <- 1
    for (g in 1:config$G) {
      K_g <- fpca_results$m_vec[g]
      fpca_g <- fpca_results$fpca_list[[g]]
      D_fpc_scores[, col_idx:(col_idx + K_g - 1)] <- fpca_g$xiEst[, 1:K_g]
      col_idx <- col_idx + K_g
    }
    
    # 第二部分: r_ik(t) 做FPCA以获得系数
    r_fpc_scores <- matrix(0, nrow = config$n, ncol = 0)
    r_fpca_list <- vector("list", config$G)
    r_m_vec <- integer(config$G)
    
    for (g in 1:config$G) {
      work_grid <- residuals_list[[g]]$grid
      r_ig <- residuals_list[[g]]$values
      
      # 将残差转换为Ly/Lt格式
      Ly_r <- vector("list", config$n)
      Lt_r <- vector("list", config$n)
      for (i in 1:config$n) {
        Ly_r[[i]] <- r_ig[i, ]
        Lt_r[[i]] <- work_grid
      }
      
      # 对残差做FPCA
      fpca_r <- FPCA(Ly = Ly_r, Lt = Lt_r,
                     list(dataType = 'Dense',
                          error = FALSE,
                          methodSelectK = "FVE",
                          FVEthreshold = 0.95,
                          maxK = 10,
                          verbose = FALSE))
      
      K_r <- ncol(fpca_r$xiEst)
      r_m_vec[g] <- K_r
      r_fpca_list[[g]] <- fpca_r
      
      # 添加到总的FPC系数矩阵
      r_fpc_scores <- cbind(r_fpc_scores, fpca_r$xiEst[, 1:K_r])
      
      cat(sprintf("      残差%d: FPCA保留%d个主成分\n", g, K_r))
    }
    
    # 构造完整的设计矩阵: [D的FPC系数, r的FPC系数] (不包含Z)
    X_design <- cbind(D_fpc_scores, r_fpc_scores)
    
    # OLS回归
    fit_final <- lm(Y ~ X_design)
    coef_final <- coef(fit_final)[-1]  # 去掉截距
    
    # 提取系数
    n_D_coef <- ncol(D_fpc_scores)
    n_r_coef <- ncol(r_fpc_scores)
    
    alpha_D <- coef_final[1:n_D_coef]  # D的FPC系数
    alpha_r <- coef_final[(n_D_coef + 1):(n_D_coef + n_r_coef)]  # r的FPC系数
    
    cat("    步骤4完成\n")
    
    # ===== 重构时变系数函数 γ_1k(t) =====
    cat("    重构时变系数函数...\n")
    beta_est_functions <- vector("list", config$G)
    col_idx <- 1
    
    for (g in 1:config$G) {
      K_g <- fpca_results$m_vec[g]
      fpca_g <- fpca_results$fpca_list[[g]]
      work_grid <- fpca_g$workGrid
      
      # 提取该组对应的α系数
      alpha_g <- alpha_D[col_idx:(col_idx + K_g - 1)]
      
      # 重构: γ_1k(t) = Σ_m α_1km φ_km(t)
      phi_g <- fpca_g$phi[, 1:K_g, drop = FALSE]
      gamma_1k_t <- as.vector(phi_g %*% alpha_g)
      
      beta_est_functions[[g]] <- list(
        time = work_grid,
        beta = gamma_1k_t
      )
      
      col_idx <- col_idx + K_g
    }
    cat("完成\n")
    
    return(list(
      success = TRUE,
      beta_est_functions = beta_est_functions,
      alpha_D = alpha_D,           # D的FPC系数
      alpha_r = alpha_r,           # 残差的FPC系数
      D_list = D_list,
      residuals_list = residuals_list,
      r_fpca_list = r_fpca_list,
      r_m_vec = r_m_vec,
      B2_estimates = B2_estimates,  # 第一阶段的时变系数
      method_type = "PACE+2SFRI"
    ))
    
  }, error = function(e) {
    cat("失败:", e$message, "\n")
    return(list(success = FALSE, error = e$message))
  })
}


#' MFMR方法 - 直接对FPCA系数矩阵进行估计
run_MFMR <- function(fpca_results, Y, instruments, config) {
  tryCatch({
    cat("  运行MFMR...")
    
    # 直接使用FPCA得到的A矩阵作为内生变量
    A_est <- fpca_results$A_est
    
    # # 最小二乘估计
    # XtX <- t(A_est) %*% A_est
    # XtY <- t(A_est) %*% Y
    # coef_final <- solve(XtX, XtY)
    Z <- instruments$Z
    
    # 用AER包的ivreg 做2SLS
    iv_fit <- ivreg(Y ~ A_est | Z)
    coef_final <- coef(iv_fit)[-1]
    
    # 系数函数重构
    est_beta_mat <- vector("list", config$G)
    col_idx <- 1
    
    for (g in 1:config$G) {
      K_g <- fpca_results$m_vec[g]
      if (col_idx <= length(coef_final)) {
        end_idx <- min(col_idx + K_g - 1, length(coef_final))
        est_beta_mat[[g]] <- coef_final[col_idx:end_idx]
        if (length(est_beta_mat[[g]]) < K_g) {
          est_beta_mat[[g]] <- c(est_beta_mat[[g]], rep(0, K_g - length(est_beta_mat[[g]])))
        }
      } else {
        est_beta_mat[[g]] <- rep(0, K_g)
      }
      col_idx <- col_idx + K_g
    }
    
    cat("完成\n")
    
    return(list(
      success = TRUE,
      fpca_list = fpca_results$fpca_list,
      est_beta_mat = est_beta_mat,
      m_vec = fpca_results$m_vec,
      method_type = "MFMR"
    ))
    
  }, error = function(e) {
    cat("失败:", e$message, "\n")
    return(list(success = FALSE, error = e$message))
  })
}


#' FA-MFMR方法 - Post-Lasso核心
run_FA_MFMR <- function(fpca_results, Y, instruments, config, factor_results) {
  tryCatch({
    cat("  运行FA-MFMR (因子结构 + Post-Lasso + 工具变量)...")
    
    # 第一阶段：2SLS
    F_iv <- matrix(0, nrow = config$n, ncol = factor_results$K_selected)
    F_stat <- numeric(factor_results$K_selected)
    
    for (k in 1:factor_results$K_selected) {
      fit <- lm(factor_results$F_hat[, k] ~ instruments$Z)
      F_iv[, k] <- fitted(fit)
      fit_summary <- summary(fit)
      F_stat[k] <- fit_summary$fstatistic[1]
    }
    
    U_iv <- matrix(0, nrow = config$n, ncol = ncol(factor_results$U_hat))
    for (j in 1:ncol(factor_results$U_hat)) {
      fit <- lm(factor_results$U_hat[, j] ~ instruments$Z)
      U_iv[, j] <- fitted(fit)
    }
    
    # 第二阶段回归：Projection + SCAD估计
    cat("\n")
    projection_results <- projection_scad_estimation(Y, F_iv, U_iv)
    coef_selected <- projection_results$beta_scad
    
    # 系数函数重构
    est_beta_mat <- vector("list", config$G)
    col_idx <- 1
    
    for (g in 1:config$G) {
      K_g <- fpca_results$m_vec[g]
      if (col_idx <= length(coef_selected)) {
        end_idx <- min(col_idx + K_g - 1, length(coef_selected))
        est_beta_mat[[g]] <- coef_selected[col_idx:end_idx]
        if (length(est_beta_mat[[g]]) < K_g) {
          est_beta_mat[[g]] <- c(est_beta_mat[[g]], rep(0, K_g - length(est_beta_mat[[g]])))
        }
      } else {
        est_beta_mat[[g]] <- rep(0, K_g)
      }
      col_idx <- col_idx + K_g
    }
    
    cat("完成\n")
    
    return(list(
      success = TRUE,
      fpca_list = fpca_results$fpca_list,
      est_beta_mat = est_beta_mat,
      m_vec = fpca_results$m_vec,
      F_stat = F_stat,
      K_selected = factor_results$K_selected,
      method_type = "FA-MFMR"
    ))
    
  }, error = function(e) {
    cat("失败:", e$message, "\n")
    return(list(success = FALSE, error = e$message))
  })
}

# ===================================================================
# -----------结果评估函数-----------
# ===================================================================

#' 计算MSE (统一处理所有方法)
calculate_mse <- function(method_results, beta_true_mat, config) {
  key_times <- config$key_times
  mse_vec <- numeric(config$G)
  
  if (!method_results$success) {
    return(list(mse_vec = rep(Inf, config$G), total_mse = Inf, 
                nonzero_mse = Inf, zero_mse = Inf))
  }
  
  # 判断是MPCMR方法还是其他方法
  is_mpcmr <- !is.null(method_results$beta_est_functions)
  
  for (g in 1:config$G) {
    if (is_mpcmr) {
      # MPCMR方法：直接使用beta_est_functions
      if (!is.null(method_results$beta_est_functions[[g]])) {
        time_points <- method_results$beta_est_functions[[g]]$time
        beta_est <- method_results$beta_est_functions[[g]]$beta
        
        # 插值真实beta到MPCMR的时间网格
        beta_true_interp <- approx(config$t_grid, beta_true_mat[, g], 
                                   xout = time_points, rule = 2)$y
        
        eval_idx <- sapply(key_times, function(tt) which.min(abs(time_points - tt)))
        eval_idx <- eval_idx[eval_idx <= length(time_points) & eval_idx > 0]
        
        if (length(eval_idx) > 0) {
          mse_vec[g] <- mean((beta_true_interp[eval_idx] - beta_est[eval_idx])^2)
        } else {
          mse_vec[g] <- mean((beta_true_interp - beta_est)^2)
        }
      } else {
        mse_vec[g] <- Inf
      }
    } else {
      # FA-MFMR和MFMR方法：使用FPCA重构
      phi_kl_g <- method_results$fpca_list[[g]]$phi[, 1:method_results$m_vec[g], drop = FALSE]
      beta_est_fpca_grid <- as.vector(phi_kl_g %*% method_results$est_beta_mat[[g]])
      
      fpca_grid <- method_results$fpca_list[[g]]$workGrid
      beta_true_interp <- approx(config$t_grid, beta_true_mat[, g], 
                                 xout = fpca_grid, rule = 2)$y
      
      eval_idx <- sapply(key_times, function(tt) which.min(abs(fpca_grid - tt)))
      eval_idx <- eval_idx[eval_idx <= length(fpca_grid) & eval_idx > 0]
      
      if (length(eval_idx) > 0) {
        mse_vec[g] <- mean((beta_true_interp[eval_idx] - beta_est_fpca_grid[eval_idx])^2)
      } else {
        mse_vec[g] <- mean((beta_true_interp - beta_est_fpca_grid)^2)
      }
    }
  }
  
  # 分别计算非零函数和零函数的MSE
  nonzero_mse <- sum(mse_vec[1:4])  # 前4个非零函数
  zero_mse <- sum(mse_vec[5:config$G])  # 后16个零函数
  
  return(list(
    mse_vec = mse_vec, 
    total_mse = sum(mse_vec),
    nonzero_mse = nonzero_mse,
    zero_mse = zero_mse
  ))
}

# ===================================================================
# -----------数据生成函数-----------
# ===================================================================

#' 完整的数据生成流程
generate_simulation_data <- function(config) {
  
  cat("开始数据生成流程...\n")
  
  # 1. 生成随机观测时间
  cat("1/6 生成随机观测时间...")
  obs_times_list <- generate_random_obs_times(
    config$n, config$t_range, config$min_obs, config$max_obs, 
    config$obs_prob_center, config$obs_prob_edge
  )
  cat("完成\n")
  
  # 2. 生成Fourier基函数
  cat("2/6 生成Fourier基函数...")
  basis_eval <- generate_bspline_basis(config$t_grid, config$m)
  cat("完成\n")
  
  # 3. 生成工具变量
  cat("3/6 生成工具变量...")
  instruments <- generate_instruments(config$n, config$J, config$K_true, config$p)
  cat("完成\n")
  
  # 4. 生成完整函数轨迹
  cat("4/6 生成完整函数轨迹...")
  trajectories <- generate_complete_trajectories(config, instruments$A_mat, basis_eval)
  X_complete_array <- trajectories$X_complete_array
  bm_t50_values <- trajectories$bm_t50_values 
  common_bm_indiv <- trajectories$common_bm_indiv
  cat("完成\n")
  
  # 5. 抽取观测数据
  cat("5/6 抽取观测数据...")
  func_data <- extract_observed_data(X_complete_array, obs_times_list, config)
  cat("完成\n")
  
  # 6. 生成真实beta函数和响应变量
  cat("6/6 生成真实beta函数和响应变量...")
  beta_true_mat <- generate_true_beta_functions(config$t_grid)
  Y <- generate_response_from_complete_trajectories(config, X_complete_array, beta_true_mat, bm_t50_values, common_bm_indiv, instruments$A_mat) 
  cat("完成\n")
  
  cat("✓ 数据生成完成！\n\n")
  
  return(list(
    func_data = func_data,
    Y = Y,
    instruments = instruments,
    beta_true_mat = beta_true_mat,
    obs_times_list = obs_times_list,
    X_complete_array = X_complete_array,
    bm_t50_values = bm_t50_values,
    common_bm_indiv = common_bm_indiv
  ))
}

# ===================================================================
# -----------可视化函数-----------
# ===================================================================

#' 可视化结果 - 所有函数在一张大图展示 (2行排列)
plot_results <- function(method_results, data, config, method_name = "Method") {
  if (!method_results$success) {
    cat(sprintf("%s失败，无法绘图\n", method_name))
    return(NULL)
  }
  
  # 判断是MPCMR方法还是其他方法
  is_mpcmr <- !is.null(method_results$beta_est_functions)
  
  # ===== 系数函数对比图：所有函数在一张图上 =====
  # 设置绘图参数：第一行3张，第二行2张
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 2), mgp = c(2.5, 1, 0))
  group_names <- c("线性", "二次函数", "三角", "高斯", "零函数")
  mse_vec <- numeric(config$G)
  
  # 只绘制前5个函数（4个非零 + 1个零函数作为代表）
  plot_indices <- c(1, 2, 3, 4, 5)
  
  for (plot_idx in 1:5) {
    g <- plot_indices[plot_idx]
    
    if (is_mpcmr) {
      # MPCMR方法
      if (!is.null(method_results$beta_est_functions[[g]])) {
        time_points <- method_results$beta_est_functions[[g]]$time
        beta_est <- method_results$beta_est_functions[[g]]$beta
        
        beta_true_interp <- approx(config$t_grid, data$beta_true_mat[, g], 
                                   xout = time_points, rule = 2)$y
        
        eval_idx <- sapply(config$key_times, function(tt) which.min(abs(time_points - tt)))
        eval_idx <- eval_idx[eval_idx <= length(time_points) & eval_idx > 0]
        
        if (length(eval_idx) > 0) {
          mse_vec[g] <- mean((beta_true_interp[eval_idx] - beta_est[eval_idx])^2)
        } else {
          mse_vec[g] <- mean((beta_true_interp - beta_est)^2)
        }
        
        # 绘图 - 纵轴范围固定为[-1, 1]
        plot(time_points, beta_true_interp, type = 'l', col = 'navy', lwd = 3,
             ylim = c(-1, 1),
             ylab = expression(beta[g](t)), 
             xlab = "时间 t",
             main = sprintf("%s (MSE=%.4f)", group_names[plot_idx], mse_vec[g]),
             cex.main = 1.0, cex.lab = 1.0)
        
        lines(time_points, beta_est, col = 'red', lwd = 2.5, lty = 2)
        
        if (length(eval_idx) > 0) {
          points(time_points[eval_idx], beta_true_interp[eval_idx], 
                 col = 'navy', pch = 16, cex = 1.0)
          points(time_points[eval_idx], beta_est[eval_idx], 
                 col = 'red', pch = 17, cex = 1.0)
        }
      }
    } else {
      # FA-MFMR和MFMR方法
      phi_kl_g <- method_results$fpca_list[[g]]$phi[, 1:method_results$m_vec[g], drop = FALSE]
      beta_est_fpca_grid <- as.vector(phi_kl_g %*% method_results$est_beta_mat[[g]])
      
      fpca_grid <- method_results$fpca_list[[g]]$workGrid
      beta_true_interp <- approx(config$t_grid, data$beta_true_mat[, g], 
                                 xout = fpca_grid, rule = 2)$y
      
      eval_idx <- sapply(config$key_times, function(tt) which.min(abs(fpca_grid - tt)))
      eval_idx <- eval_idx[eval_idx <= length(fpca_grid) & eval_idx > 0]
      
      if (length(eval_idx) > 0) {
        mse_vec[g] <- mean((beta_true_interp[eval_idx] - beta_est_fpca_grid[eval_idx])^2)
      } else {
        mse_vec[g] <- mean((beta_true_interp - beta_est_fpca_grid)^2)
      }
      
      # 绘图 - 纵轴范围固定为[-1, 1]
      plot(fpca_grid, beta_true_interp, type = 'l', col = 'navy', lwd = 3,
           ylim = c(-1, 1),
           ylab = expression(beta[g](t)), 
           xlab = "时间 t",
           main = sprintf("%s (MSE=%.4f)", group_names[plot_idx], mse_vec[g]),
           cex.main = 1.0, cex.lab = 1.0)
      
      lines(fpca_grid, beta_est_fpca_grid, col = 'red', lwd = 2.5, lty = 2)
      
      if (length(eval_idx) > 0) {
        points(fpca_grid[eval_idx], beta_true_interp[eval_idx], 
               col = 'navy', pch = 16, cex = 1.0)
        points(fpca_grid[eval_idx], beta_est_fpca_grid[eval_idx], 
               col = 'red', pch = 17, cex = 1.0)
      }
    }
    
    # 添加网格
    grid(col = "lightgray", lty = 3)
    
    # 只在第一个子图添加图例
    if (plot_idx == 1) {
      legend("topright", 
             legend = c("真实值", "估计值", "评估点"),
             col = c("navy", "red", "gray"), 
             lty = c(1, 2, 0), 
             pch = c(NA, NA, 16),
             lwd = c(3, 2.5, NA),
             bty = "n", cex = 0.8)
    }
  }
  
  # 第二行第三个位置添加图例和统计信息
  par(mar = c(1, 1, 1, 1))
  plot.new()
  
  # 计算整体统计信息
  total_mse <- sum(mse_vec)
  nonzero_mse <- sum(mse_vec[1:4])
  zero_mse <- sum(mse_vec[5:config$G])
  
  # 添加标题
  text(0.5, 0.95, paste(method_name, "结果总结", sep = "\n"), 
       cex = 1.2, font = 2, adj = c(0.5, 0.5))
  
  # 添加详细图例
  legend(0.5, 0.75, 
         legend = c("真实函数 (蓝色实线)", 
                    "估计函数 (红色虚线)", 
                    "关键评估点 (菱形标记)"),
         col = c("navy", "red", "black"), 
         lty = c(1, 2, 0), 
         pch = c(NA, NA, 17),
         lwd = c(3, 2.5, NA),
         bty = "round", cex = 1.0, 
         x.intersp = 1.5, y.intersp = 1.5,
         xjust = 0.5, yjust = 0.5)
  
  # 添加统计信息
  info_text <- sprintf(
    "总体MSE: %.6f\n\n非零函数MSE :\n%.6f\n\n零函数MSE :\n%.6f",
    total_mse, nonzero_mse, zero_mse
  )
  
  text(0.5, 0.25, info_text, 
       cex = 0.95, adj = c(0.5, 0.5))
  
  return(list(
    mse_vec = mse_vec,
    total_mse = total_mse,
    nonzero_mse = nonzero_mse,
    zero_mse = zero_mse
  ))
}

#' 可视化真实Beta函数 (重点显示前4个非零函数)
plot_true_beta_functions <- function(config) {
  # 生成真实beta函数
  beta_true_mat <- generate_true_beta_functions(config$t_grid)
  
  # 设置绘图参数 - 显示前4个非零函数
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  group_names <- c("Beta1 (线性)", "Beta2 (二次函数)", "Beta3 (三角)", "Beta4 (高斯)")
  colors <- c("red", "blue", "green", "orange")
  
  for (g in 1:4) {  # 只显示前4个非零函数
    plot(config$t_grid, beta_true_mat[, g], 
         type = 'l', col = colors[g], lwd = 3,
         ylab = expression(beta[g](t)), 
         xlab = "时间 t",
         main = group_names[g],
         cex.main = 1.1, cex.lab = 1.1)
    
    # 标记关键时间点
    key_idx <- sapply(config$key_times, function(tt) which.min(abs(config$t_grid - tt)))
    points(config$t_grid[key_idx], beta_true_mat[key_idx, g], 
           col = colors[g], pch = 16, cex = 1.2)
    
    grid(col = "lightgray", lty = 3)
  }
  
  cat(sprintf("注意: Beta5-Beta20 全部为零函数 (共%d个零函数)\n", config$G - 4))
  
  return(beta_true_mat)
}


# ===================================================================
# -----------CI计算函数-----------
# ===================================================================

#' 计算MFMR的Wald型理论CI
compute_ci_mfmr <- function(Y, A_est, Z, fpca_results, config, alpha = 0.05) {
  
  cat("    计算MFMR的Wald型理论置信区间...\n")
  
  # 2SLS回归
  iv_fit <- ivreg(Y ~ A_est | Z)
  coef_2sls <- coef(iv_fit)[-1]
  
  # 计算方差
  n <- length(Y)
  k <- ncol(A_est)
  
  residuals_iv <- residuals(iv_fit)
  sigma2_hat <- sum(residuals_iv^2) / (n - k)
  
  # 第一阶段拟合值
  first_stage_fit <- lm(A_est ~ Z)
  A_hat <- fitted(first_stage_fit)
  
  # 方差-协方差矩阵
  var_cov <- sigma2_hat * solve(crossprod(A_hat))
  se_2sls <- sqrt(diag(var_cov))
  
  # 计算CI
  z_critical <- qnorm(1 - alpha/2)
  ci_low <- coef_2sls - z_critical * se_2sls
  ci_up <- coef_2sls + z_critical * se_2sls
  
  # 按Beta函数分组构造CI
  beta_ci_functions <- vector("list", config$G)
  col_idx <- 1
  
  for (g in 1:config$G) {
    K_g <- fpca_results$m_vec[g]
    fpca_g <- fpca_results$fpca_list[[g]]
    work_grid <- fpca_g$workGrid
    
    if (col_idx <= length(coef_2sls)) {
      end_idx <- min(col_idx + K_g - 1, length(coef_2sls))
      coef_g <- coef_2sls[col_idx:end_idx]
      se_g <- se_2sls[col_idx:end_idx]
      
      if (length(coef_g) < K_g) {
        coef_g <- c(coef_g, rep(0, K_g - length(coef_g)))
        se_g <- c(se_g, rep(Inf, K_g - length(se_g)))
      }
    } else {
      coef_g <- rep(0, K_g)
      se_g <- rep(Inf, K_g)
    }
    
    # 重构时间轴上的CI - 误差传播
    phi_g <- fpca_g$phi[, 1:K_g, drop = FALSE]
    beta_g <- as.vector(phi_g %*% coef_g)
    
    se_pointwise <- numeric(nrow(phi_g))
    for (t_idx in 1:nrow(phi_g)) {
      se_pointwise[t_idx] <- sqrt(sum((phi_g[t_idx, ] * se_g)^2))
    }
    
    beta_ci_functions[[g]] <- list(
      time = work_grid,
      beta = beta_g,
      se = se_pointwise,
      ci_low = beta_g - z_critical * se_pointwise,
      ci_up = beta_g + z_critical * se_pointwise
    )
    
    col_idx <- col_idx + K_g
  }
  
  cat("      完成\n")
  
  return(list(
    beta_ci_functions = beta_ci_functions,
    alpha = alpha,
    method = "MFMR_Wald"
  ))
}

#' 转换MPCMR的CI到系数函数
convert_mpcmr_ci <- function(mpcmr_results, fpca_results, config, alpha = 0.05) {
  
  cat("    转换MPCMR的GMM理论CI到系数函数空间...\n")
  
  beta_ci_functions <- vector("list", config$G)
  
  for (g in 1:config$G) {
    if (!is.null(mpcmr_results$beta_est_functions[[g]])) {
      beta_ci_functions[[g]] <- list(
        time = mpcmr_results$beta_est_functions[[g]]$time,
        beta = mpcmr_results$beta_est_functions[[g]]$beta,
        se = sqrt(mpcmr_results$beta_est_functions[[g]]$ci_var),
        ci_low = mpcmr_results$beta_est_functions[[g]]$ci_low_wald,
        ci_up = mpcmr_results$beta_est_functions[[g]]$ci_up_wald
      )
    } else {
      beta_ci_functions[[g]] <- list(
        time = NA,
        beta = NA,
        se = NA,
        ci_low = NA,
        ci_up = NA
      )
    }
  }
  
  cat("      完成\n")
  
  return(list(
    beta_ci_functions = beta_ci_functions,
    alpha = alpha,
    method = "MPCMR_GMM"
  ))
}

#' Bootstrap FPCA - FA-MFMR和PACE+2SFRI共用
bootstrap_fpca_results <- function(func_data_original, obs_times_list_original, 
                                   config, n_bootstrap = 50,
                                   Y_original = NULL,
                                   instruments_original = NULL) {
  
  cat(sprintf("    执行%d次功能型Bootstrap (重计算FPCA + 重抽样Y和instruments)...\n", n_bootstrap))
  
  fpca_bootstrap_list <- vector("list", n_bootstrap)
  
  # 【修复】直接初始化为列表，不用 if-else
  Y_bootstrap_list <- vector("list", n_bootstrap)
  instruments_bootstrap_list <- vector("list", n_bootstrap)
  
  for (b in 1:n_bootstrap) {
    if (b %% 50 == 0) cat(sprintf("      Bootstrap迭代 %d/%d\n", b, n_bootstrap))
    
    # ========== 【修复】正确的有放回抽样 ==========
    idx_boot <- sample(1:config$n, size = config$n, replace = TRUE)  # 注意: 应该是 replace = TRUE
    
    # 构造Bootstrap样本的Ly/Lt数据
    Ly_boot <- vector("list", config$G)
    Lt_boot <- vector("list", config$G)
    
    for (g in 1:config$G) {
      Ly_boot[[g]] <- vector("list", config$n)
      Lt_boot[[g]] <- vector("list", config$n)
    }
    
    for (i in 1:config$n) {
      i_boot <- idx_boot[i]
      
      for (g in 1:config$G) {
        Ly_boot[[g]][[i]] <- func_data_original$Ly_list[[g]][[i_boot]]
        Lt_boot[[g]][[i]] <- func_data_original$Lt_list[[g]][[i_boot]]
      }
    }
    
    # ========== 【新增】Bootstrap 样本对应的 Y ==========
    if (!is.null(Y_original)) {
      Y_boot <- Y_original[idx_boot]
      Y_bootstrap_list[[b]] <- Y_boot
    }
    
    # ========== 【新增】Bootstrap 样本对应的 instruments ==========
    if (!is.null(instruments_original)) {
      Z_boot <- instruments_original$Z[idx_boot, , drop = FALSE]
      instruments_boot <- list(
        Z = Z_boot,
        A_mat = instruments_original$A_mat[idx_boot, , drop = FALSE],
        F_mat = instruments_original$F_mat[idx_boot, , drop = FALSE],
        U_mat = instruments_original$U_mat[idx_boot, , drop = FALSE],
        Pi_F = instruments_original$Pi_F,
        Pi_U = instruments_original$Pi_U
      )
      instruments_bootstrap_list[[b]] <- instruments_boot
    }
    
    # 对Bootstrap样本计算FPCA
    tryCatch({
      fpca_g_boot <- vector("list", config$G)
      m_vec_boot <- integer(config$G)
      
      for (g in 1:config$G) {
        fpca_g_boot[[g]] <- FPCA(
          Ly = Ly_boot[[g]], 
          Lt = Lt_boot[[g]], 
          list(dataType = 'Sparse', 
               error = TRUE,
               methodSelectK = "FVE",
               FVEthreshold = 0.95,
               maxK = 15,
               userBwMu = 1.5,
               userBwCov = 2.0,
               nRegGrid = 101,
               verbose = FALSE)
        )
        
        FVE <- fpca_g_boot[[g]]$cumFVE
        K_g <- which(FVE >= 0.95)[1]
        if (is.na(K_g)) K_g <- min(10, length(FVE))
        m_vec_boot[g] <- K_g
      }
      
      # 构造Bootstrap样本的FPCA系数矩阵
      m_total_boot <- sum(m_vec_boot)
      A_est_boot <- matrix(0, nrow = config$n, ncol = m_total_boot)
      col_start <- 1
      
      for (g in 1:config$G) {
        K_g <- m_vec_boot[g]
        col_end <- col_start + K_g - 1
        A_est_boot[, col_start:col_end] <- fpca_g_boot[[g]]$xiEst[, 1:K_g, drop = FALSE]
        col_start <- col_end + 1
      }
      
      fpca_bootstrap_list[[b]] <- list(
        fpca_list = fpca_g_boot,
        m_vec = m_vec_boot,
        A_est = A_est_boot,
        success = TRUE
      )
      
    }, error = function(e) {
      fpca_bootstrap_list[[b]] <<- list(success = FALSE, error = e$message)
    })
  }
  
  cat("      完成\n")
  
  return(list(
    fpca_bootstrap_list = fpca_bootstrap_list,
    Y_bootstrap_list = Y_bootstrap_list,
    instruments_bootstrap_list = instruments_bootstrap_list
  ))
}


#' FA-MFMR的Bootstrap CI
compute_ci_FA_MFMR <- function(Y_original, F_iv_original, U_iv_original,
                               fpca_results, config, instruments,
                               fpca_bootstrap_list, factor_results,
                               alpha = 0.05) {
  
  cat("    计算FA-MFMR的Bootstrap置信区间...\n")
  
  n <- config$n
  
  # ========== 【修复】检查 fpca_bootstrap_list 的结构 ==========
  # 如果是新的返回格式（包含Y和instruments）
  if (!is.null(fpca_bootstrap_list$fpca_bootstrap_list)) {
    fpca_list_actual <- fpca_bootstrap_list$fpca_bootstrap_list
    Y_bootstrap_list <- fpca_bootstrap_list$Y_bootstrap_list
    instruments_bootstrap_list <- fpca_bootstrap_list$instruments_bootstrap_list
  } else {
    # 向后兼容：如果是旧的格式
    fpca_list_actual <- fpca_bootstrap_list
    Y_bootstrap_list <- NULL
    instruments_bootstrap_list <- NULL
  }
  
  n_bootstrap <- length(fpca_list_actual)
  
  coef_boot_list <- vector("list", n_bootstrap)
  
  for (b in 1:n_bootstrap) {
    if (b %% 50 == 0) cat(sprintf("      FA-MFMR Bootstrap %d/%d\n", b, n_bootstrap))
    
    if (!fpca_list_actual[[b]]$success) next
    
    fpca_b <- fpca_list_actual[[b]]
    A_est_b <- fpca_b$A_est
    
    # ========== 【修复】使用 Bootstrap 的 Y 和 instruments ==========
    Y_use <- if (!is.null(Y_bootstrap_list)) Y_bootstrap_list[[b]] else Y_original
    instruments_use <- if (!is.null(instruments_bootstrap_list)) instruments_bootstrap_list[[b]] else instruments
    
    tryCatch({
      # Bootstrap样本的因子分解
      A_centered_b <- scale(A_est_b, center = TRUE, scale = FALSE)
      K_max_search <- min(10, floor(ncol(A_est_b)/2))
      K_selected_b <- select_factors_ratio_method(A_centered_b, K_max_search)
      
      # SVD分解
      svd_A_b <- svd(A_centered_b)
      F_hat_b <- sqrt(n) * svd_A_b$u[, 1:K_selected_b, drop = FALSE]
      
      D_k_b <- diag(svd_A_b$d[1:K_selected_b], nrow = K_selected_b)
      V_k_b <- svd_A_b$v[, 1:K_selected_b, drop = FALSE]
      B_hat_b <- sqrt(n) * V_k_b %*% D_k_b / n
      
      low_rank_part_b <- F_hat_b %*% t(B_hat_b)
      U_hat_b <- A_centered_b - low_rank_part_b
      
      # 2SLS 第一阶段 - 【修复】使用 Bootstrap instruments
      F_iv_b <- matrix(0, nrow = n, ncol = K_selected_b)
      for (k in 1:K_selected_b) {
        fit_k <- lm(F_hat_b[, k] ~ instruments_use$Z)
        F_iv_b[, k] <- fitted(fit_k)
      }
      
      U_iv_b <- matrix(0, nrow = n, ncol = ncol(U_hat_b))
      for (j in 1:ncol(U_hat_b)) {
        fit_j <- lm(U_hat_b[, j] ~ instruments_use$Z)
        U_iv_b[, j] <- fitted(fit_j)
      }
      
      # 投影方法 + SCAD - 【修复】使用 Bootstrap Y
      fit_Y_F_b <- lm(Y_use ~ F_iv_b - 1)
      Y_residual_b <- residuals(fit_Y_F_b)
      
      U_iv_residual_b <- matrix(0, nrow = n, ncol = ncol(U_iv_b))
      for (j in 1:ncol(U_iv_b)) {
        fit_U_F_b <- lm(U_iv_b[, j] ~ F_iv_b - 1)
        U_iv_residual_b[, j] <- residuals(fit_U_F_b)
      }
      
      # SCAD回归
      cv_fit_b <- cv.ncvreg(U_iv_residual_b, Y_residual_b, penalty = "SCAD", 
                            nfolds = 10, seed = 123)
      lambda_opt_b <- cv_fit_b$lambda.min
      lambda_index_b <- which.min(abs(cv_fit_b$lambda - lambda_opt_b))
      beta_scad_b <- as.vector(cv_fit_b$fit$beta[-1, lambda_index_b])
      
      coef_boot_list[[b]] <- list(
        coef = beta_scad_b,
        fpca = fpca_b,
        success = TRUE
      )
      
    }, error = function(e) {
      coef_boot_list[[b]] <<- list(success = FALSE, error = e$message)
    })
  }
  
  # 从Bootstrap系数重构时间轴上的CI
  z_critical <- qnorm(1 - alpha/2)
  beta_ci_functions <- vector("list", config$G)
  
  for (g in 1:config$G) {
    time_points <- fpca_results$fpca_list[[g]]$workGrid
    n_time <- length(time_points)
    K_g_original <- fpca_results$m_vec[g]
    
    beta_g_boot_mat <- matrix(NA, nrow = n_bootstrap, ncol = n_time)
    
    for (b in 1:n_bootstrap) {
      if (!coef_boot_list[[b]]$success) next
      
      coef_b <- coef_boot_list[[b]]$coef
      fpca_b <- coef_boot_list[[b]]$fpca
      
      # 定位第g组的系数
      if (g == 1) {
        col_start <- 1
      } else {
        col_start <- sum(fpca_b$m_vec[1:(g-1)]) + 1
      }
      col_end <- sum(fpca_b$m_vec[1:g])
      K_g_b <- fpca_b$m_vec[g]
      
      if (col_start <= length(coef_b)) {
        end_idx <- min(col_end, length(coef_b))
        coef_g_b <- coef_b[col_start:end_idx]
        
        if (length(coef_g_b) < K_g_b) {
          coef_g_b <- c(coef_g_b, rep(0, K_g_b - length(coef_g_b)))
        }
        
        # 使用原始FPCA的基函数在该时间点重构
        phi_g_boot <- fpca_b$fpca_list[[g]]$phi[, 1:K_g_b, drop = FALSE]
        grid_boot <- fpca_b$fpca_list[[g]]$workGrid
        
        # 对每个原始时间点插值Bootstrap基函数
        beta_g_boot_t <- numeric(n_time)
        for (t_idx in 1:n_time) {
          phi_interp_values <- numeric(K_g_b)
          for (k in 1:K_g_b) {
            phi_interp_values[k] <- approx(grid_boot, phi_g_boot[, k], 
                                           xout = time_points[t_idx], rule = 2)$y
          }
          beta_g_boot_t[t_idx] <- sum(phi_interp_values * coef_g_b[1:K_g_b])
        }
        
        beta_g_boot_mat[b, ] <- beta_g_boot_t
      }
    }
    
    # 计算百分位法CI
    ci_low <- apply(beta_g_boot_mat, 2, function(x) quantile(x, alpha/2, na.rm = TRUE))
    ci_up <- apply(beta_g_boot_mat, 2, function(x) quantile(x, 1 - alpha/2, na.rm = TRUE))
    beta_point <- apply(beta_g_boot_mat, 2, mean, na.rm = TRUE)
    se_boot <- apply(beta_g_boot_mat, 2, sd, na.rm = TRUE)
    
    beta_ci_functions[[g]] <- list(
      time = time_points,
      beta = beta_point,
      se = se_boot,
      ci_low = ci_low,
      ci_up = ci_up
    )
  }
  
  cat("      完成\n")
  
  return(list(
    beta_ci_functions = beta_ci_functions,
    alpha = alpha,
    n_bootstrap = n_bootstrap,
    method = "FA_MFMR_Bootstrap"
  ))
}

#' PACE+2SFRI的Bootstrap CI
compute_ci_pace_2sfri <- function(Y_original, func_data_original, 
                                  obs_times_list_original, config, 
                                  instruments, fpca_bootstrap_list, 
                                  fpca_results, alpha = 0.05) {
  
  cat("    计算PACE+2SFRI的Bootstrap置信区间...\n")
  
  n <- config$n
  
  # ========== 【修复】检查 fpca_bootstrap_list 的结构 ==========
  if (!is.null(fpca_bootstrap_list$fpca_bootstrap_list)) {
    fpca_list_actual <- fpca_bootstrap_list$fpca_bootstrap_list
    Y_bootstrap_list <- fpca_bootstrap_list$Y_bootstrap_list
    instruments_bootstrap_list <- fpca_bootstrap_list$instruments_bootstrap_list
  } else {
    fpca_list_actual <- fpca_bootstrap_list
    Y_bootstrap_list <- NULL
    instruments_bootstrap_list <- NULL
  }
  
  n_bootstrap <- length(fpca_list_actual)
  
  coef_boot_list <- vector("list", n_bootstrap)
  
  for (b in 1:n_bootstrap) {
    if (b %% 50 == 0) cat(sprintf("      PACE+2SFRI Bootstrap %d/%d\n", b, n_bootstrap))
    
    if (!fpca_list_actual[[b]]$success) next
    
    fpca_b <- fpca_list_actual[[b]]
    
    # ========== 【修复】使用 Bootstrap 的 Y 和 instruments ==========
    Y_use <- if (!is.null(Y_bootstrap_list)) Y_bootstrap_list[[b]] else Y_original
    instruments_use <- if (!is.null(instruments_bootstrap_list)) instruments_bootstrap_list[[b]] else instruments
    
    tryCatch({
      # 第一阶段：计算D(t)和r(t)
      D_list_b <- vector("list", config$G)
      residuals_list_b <- vector("list", config$G)
      
      for (g in 1:config$G) {
        fpca_g_b <- fpca_b$fpca_list[[g]]
        mu_g <- fpca_g_b$mu
        work_grid <- fpca_g_b$workGrid
        K_g <- fpca_b$m_vec[g]
        
        # 重建D(t)
        D_ig_b <- matrix(0, nrow = n, ncol = length(work_grid))
        for (i in 1:n) {
          xi_ig <- fpca_g_b$xiEst[i, 1:K_g]
          phi_g <- fpca_g_b$phi[, 1:K_g, drop = FALSE]
          x_ig <- mu_g + as.vector(phi_g %*% xi_ig)
          D_ig_b[i, ] <- x_ig - mu_g
        }
        
        D_list_b[[g]] <- list(values = D_ig_b, grid = work_grid)
        
        # 第一阶段回归：D(t) ~ Z - 【修复】使用 Bootstrap instruments
        r_ig_b <- matrix(0, nrow = n, ncol = length(work_grid))
        for (t_idx in 1:length(work_grid)) {
          D_it <- D_ig_b[, t_idx]
          fit_t <- lm(D_it ~ instruments_use$Z)
          r_ig_b[, t_idx] <- residuals(fit_t)
        }
        
        residuals_list_b[[g]] <- list(values = r_ig_b, grid = work_grid)
      }
      
      # 第二阶段：D(t)的FPC系数矩阵
      D_fpc_scores_b <- matrix(0, nrow = n, ncol = sum(fpca_b$m_vec))
      col_idx <- 1
      for (g in 1:config$G) {
        K_g <- fpca_b$m_vec[g]
        D_fpc_scores_b[, col_idx:(col_idx + K_g - 1)] <- fpca_b$fpca_list[[g]]$xiEst[, 1:K_g]
        col_idx <- col_idx + K_g
      }
      
      # r(t)的FPCA
      r_fpc_scores_b <- matrix(0, nrow = n, ncol = 0)
      for (g in 1:config$G) {
        work_grid <- residuals_list_b[[g]]$grid
        r_ig <- residuals_list_b[[g]]$values
        
        Ly_r <- vector("list", n)
        Lt_r <- vector("list", n)
        for (i in 1:n) {
          Ly_r[[i]] <- r_ig[i, ]
          Lt_r[[i]] <- work_grid
        }
        
        fpca_r_b <- FPCA(Ly = Ly_r, Lt = Lt_r,
                         list(dataType = 'Dense',
                              error = FALSE,
                              methodSelectK = "FVE",
                              FVEthreshold = 0.95,
                              maxK = 10,
                              verbose = FALSE))
        
        K_r <- ncol(fpca_r_b$xiEst)
        r_fpc_scores_b <- cbind(r_fpc_scores_b, fpca_r_b$xiEst[, 1:K_r])
      }
      
      # 第二阶段OLS - 【修复】使用 Bootstrap Y
      X_design_b <- cbind(D_fpc_scores_b, r_fpc_scores_b)
      fit_final_b <- lm(Y_use ~ X_design_b)
      coef_final_b <- coef(fit_final_b)[-1]
      
      coef_boot_list[[b]] <- list(
        coef = coef_final_b,
        fpca = fpca_b,
        D_fpc_dim = ncol(D_fpc_scores_b),
        success = TRUE
      )
      
    }, error = function(e) {
      coef_boot_list[[b]] <<- list(success = FALSE, error = e$message)
    })
  }
  
  # 从Bootstrap系数重构时间轴上的CI
  z_critical <- qnorm(1 - alpha/2)
  beta_ci_functions <- vector("list", config$G)
  
  for (g in 1:config$G) {
    time_points <- fpca_results$fpca_list[[g]]$workGrid
    n_time <- length(time_points)
    K_g_original <- fpca_results$m_vec[g]
    
    beta_g_boot_mat <- matrix(NA, nrow = n_bootstrap, ncol = n_time)
    
    for (b in 1:n_bootstrap) {
      if (!coef_boot_list[[b]]$success) next
      
      coef_b <- coef_boot_list[[b]]$coef
      fpca_b <- coef_boot_list[[b]]$fpca
      D_fpc_dim_b <- coef_boot_list[[b]]$D_fpc_dim
      
      # 定位第g组的D的系数
      if (g == 1) {
        col_start <- 1
      } else {
        col_start <- sum(fpca_b$m_vec[1:(g-1)]) + 1
      }
      col_end <- sum(fpca_b$m_vec[1:g])
      K_g_b <- fpca_b$m_vec[g]
      
      if (col_start <= D_fpc_dim_b) {
        end_idx <- min(col_end, D_fpc_dim_b)
        coef_g_b <- coef_b[col_start:end_idx]
        
        if (length(coef_g_b) < K_g_b) {
          coef_g_b <- c(coef_g_b, rep(0, K_g_b - length(coef_g_b)))
        }
        
        # 使用原始FPCA基函数重构
        phi_g_boot <- fpca_b$fpca_list[[g]]$phi[, 1:K_g_b, drop = FALSE]
        grid_boot <- fpca_b$fpca_list[[g]]$workGrid
        
        beta_g_boot_t <- numeric(n_time)
        for (t_idx in 1:n_time) {
          phi_interp_values <- numeric(K_g_b)
          for (k in 1:K_g_b) {
            phi_interp_values[k] <- approx(grid_boot, phi_g_boot[, k], 
                                           xout = time_points[t_idx], rule = 2)$y
          }
          beta_g_boot_t[t_idx] <- sum(phi_interp_values * coef_g_b[1:K_g_b])
        }
        
        beta_g_boot_mat[b, ] <- beta_g_boot_t
      }
    }
    
    # 计算百分位法CI
    ci_low <- apply(beta_g_boot_mat, 2, function(x) quantile(x, alpha/2, na.rm = TRUE))
    ci_up <- apply(beta_g_boot_mat, 2, function(x) quantile(x, 1 - alpha/2, na.rm = TRUE))
    beta_point <- apply(beta_g_boot_mat, 2, mean, na.rm = TRUE)
    se_boot <- apply(beta_g_boot_mat, 2, sd, na.rm = TRUE)
    
    beta_ci_functions[[g]] <- list(
      time = time_points,
      beta = beta_point,
      se = se_boot,
      ci_low = ci_low,
      ci_up = ci_up
    )
  }
  
  cat("      完成\n")
  
  return(list(
    beta_ci_functions = beta_ci_functions,
    alpha = alpha,
    n_bootstrap = n_bootstrap,
    method = "PACE_2SFRI_Bootstrap"
  ))
}

#' 计算CI宽度
calculate_ci_width <- function(method_ci_results, config, key_times = NULL) {
  
  if (is.null(key_times)) {
    key_times <- config$key_times
  }
  
  width_by_group <- numeric(config$G)
  
  for (g in 1:config$G) {
    if (!is.null(method_ci_results$beta_ci_functions[[g]])) {
      
      ci_info <- method_ci_results$beta_ci_functions[[g]]
      time_points <- ci_info$time
      ci_low <- ci_info$ci_low
      ci_up <- ci_info$ci_up
      
      eval_idx <- sapply(key_times, function(tt) which.min(abs(time_points - tt)))
      eval_idx <- eval_idx[eval_idx <= length(time_points) & eval_idx > 0]
      eval_idx <- unique(eval_idx)
      
      if (length(eval_idx) > 0) {
        widths <- ci_up[eval_idx] - ci_low[eval_idx]
        width_by_group[g] <- mean(widths, na.rm = TRUE)
      } else {
        widths <- ci_up - ci_low
        width_by_group[g] <- mean(widths, na.rm = TRUE)
      }
    } else {
      width_by_group[g] <- NA
    }
  }
  
  width_nonzero <- mean(width_by_group[1:4], na.rm = TRUE)
  width_zero <- mean(width_by_group[5:config$G], na.rm = TRUE)
  width_overall <- mean(width_by_group, na.rm = TRUE)
  
  return(list(
    width_by_group = width_by_group,
    width_nonzero = width_nonzero,
    width_zero = width_zero,
    width_overall = width_overall
  ))
}

#' 绘制带置信区间的系数函数
plot_beta_with_ci <- function(method_ci_results, true_beta_mat, config, 
                              method_name = "Method", plot_indices = 1:5) {
  
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 2), mgp = c(2.5, 1, 0))
  group_names <- c("线性", "二次函数", "三角", "高斯", "零函数")
  
  for (plot_idx in plot_indices) {
    if (plot_idx > config$G) break
    
    g <- plot_idx
    
    if (!is.null(method_ci_results$beta_ci_functions[[g]])) {
      
      ci_info <- method_ci_results$beta_ci_functions[[g]]
      time_points <- ci_info$time
      beta_est <- ci_info$beta
      ci_low <- ci_info$ci_low
      ci_up <- ci_info$ci_up
      
      true_beta_interp <- approx(config$t_grid, true_beta_mat[, g], 
                                 xout = time_points, rule = 2)$y
      
      # ===== 【修改】计算合适的纵坐标范围 =====
      # 收集所有需要展示的值
      all_values <- c(true_beta_interp, beta_est, ci_low, ci_up)
      all_values <- all_values[!is.na(all_values) & is.finite(all_values)]
      
      if (length(all_values) > 0) {
        y_min <- min(all_values)
        y_max <- max(all_values)
        y_range <- y_max - y_min
        
        # 添加10%的边距，确保曲线不贴边
        margin <- max(0.1 * y_range, 0.1)  # 至少0.1的边距
        ylim_use <- c(y_min - margin, y_max + margin)
      } else {
        # 如果没有有效值，使用默认范围
        ylim_use <- c(-1, 1)
      }
      
      plot(time_points, true_beta_interp, type = 'l', col = 'navy', lwd = 3,
           ylim = ylim_use,  # 使用自适应范围
           ylab = expression(beta[g](t)), 
           xlab = "时间 t",
           main = sprintf("%s (%s)", group_names[g], method_name),
           cex.main = 1.0, cex.lab = 1.0)
      
      # 添加CI区间（阴影）
      polygon(c(time_points, rev(time_points)), 
              c(ci_low, rev(ci_up)),
              col = rgb(1, 0, 0, 0.2), border = NA)
      
      # 添加估计值和真实值
      lines(time_points, beta_est, col = 'red', lwd = 2.5, lty = 2)
      lines(time_points, true_beta_interp, col = 'navy', lwd = 3)
      
      grid(col = "lightgray", lty = 3)
      
      if (plot_idx == 1) {
        legend("topright", 
               legend = c("真实值", "估计值", "95% CI"),
               col = c("navy", "red", "red"), 
               lty = c(1, 2, 1),
               fill = c(NA, NA, rgb(1, 0, 0, 0.2)),
               lwd = c(3, 2.5, NA),
               bty = "n", cex = 0.8)
      }
    }
  }
  
  par(mar = c(1, 1, 1, 1))
  plot.new()
  
  width_info <- calculate_ci_width(method_ci_results, config)
  
  info_text <- sprintf(
    "CI平均宽度:\n\n总体: %.4f\n\n非零函数: %.4f\n\n零函数: %.4f",
    width_info$width_overall,
    width_info$width_nonzero,
    width_info$width_zero
  )
  
  text(0.5, 0.5, info_text, cex = 1.2, adj = c(0.5, 0.5))
}

# ===================================================================
# -----------主执行函数：单次仿真（简化版）-----------
# ===================================================================

#' 单次仿真（简化版，不绘图）
#' 单次仿真（优化版 - 共享bootstrap样本）
run_single_simulation <- function(sim_id, config, verbose = FALSE, 
                                  compute_ci = TRUE, n_bootstrap = 5) {
  
  if (verbose) {
    cat(sprintf("\n运行第 %d 次仿真...\n", sim_id))
  }
  
  # 更新种子确保每次不同
  set.seed(config$seed + sim_id, kind = "L'Ecuyer-CMRG")
  
  tryCatch({
    # 1. 生成数据
    if (verbose) cat("  生成数据...\n")
    data <- generate_simulation_data(config)
    
    # 2. 计算FPCA
    if (verbose) cat("  计算FPCA...\n")
    fpca_results <- compute_fpca_for_all_groups(data$func_data, config)
    
    if (!fpca_results$success) {
      return(list(success = FALSE, error = "FPCA失败"))
    }
    
    # 3. 因子模型估计
    if (verbose) cat("  因子模型估计...\n")
    factor_results <- perform_factor_model_estimation(fpca_results, config)
    
    # 4. 运行四种方法
    if (verbose) cat("  运行FA-MFMR...\n")
    FA_MFMR_results <- run_FA_MFMR(fpca_results, data$Y, data$instruments, 
                                   config, factor_results)
    
    if (verbose) cat("  运行MFMR...\n")
    MFMR_results <- run_MFMR(fpca_results, data$Y, data$instruments, config)
    
    if (verbose) cat("  运行MPCMR...\n")
    MPCMR_results <- run_mpcmr(fpca_results, data$Y, data$instruments, config)
    
    if (verbose) cat("  运行PACE+2SFRI...\n")
    PACE_2SFRI_results <- run_pace_2sfri(fpca_results, data$Y, data$instruments, config)
    
    # 5. 计算MSE
    if (verbose) cat("  计算MSE...\n")
    lm_mse <- calculate_mse(FA_MFMR_results, data$beta_true_mat, config)
    MFMR_mse <- calculate_mse(MFMR_results, data$beta_true_mat, config)
    MPCMR_mse <- calculate_mse(MPCMR_results, data$beta_true_mat, config)
    PACE_2SFRI_mse <- calculate_mse(PACE_2SFRI_results, data$beta_true_mat, config)
    
    # ========== 【优化】6. 生成一次共享的bootstrap样本 ==========
    FA_MFMR_coverage <- NULL
    MFMR_coverage <- NULL
    MPCMR_coverage <- NULL
    PACE_2SFRI_coverage <- NULL
    
    fpca_boot_shared <- NULL  # 共享的bootstrap结果
    
    if (compute_ci) {
      if (verbose) cat("  计算置信区间...\n")
      
      # 【关键改进】只生成一次bootstrap样本和对应的Y、instruments
      if (verbose) cat("    生成共享的Bootstrap样本...\n")
      fpca_boot_shared <- bootstrap_fpca_results(
        data$func_data, 
        data$obs_times_list,
        config, 
        n_bootstrap = n_bootstrap,
        Y_original = data$Y,
        instruments_original = data$instruments
      )
      
      # ========== FA-MFMR: Bootstrap CI ==========
      tryCatch({
        if (verbose) cat("    FA-MFMR Bootstrap CI...\n")
        FA_MFMR_ci <- compute_ci_FA_MFMR(
          data$Y, 
          matrix(0, config$n, 1), 
          matrix(0, config$n, 1),
          fpca_results, 
          config, 
          data$instruments,
          fpca_boot_shared,  # 【改进】使用共享的bootstrap样本
          factor_results
        )
        FA_MFMR_coverage <- calculate_coverage(FA_MFMR_ci, data$beta_true_mat, config)
      }, error = function(e) {
        if (verbose) cat("    FA-MFMR CI 计算失败:", e$message, "\n")
      })
      
      # ========== MFMR: Wald型CI ==========
      tryCatch({
        if (verbose) cat("    MFMR Wald CI...\n")
        MFMR_ci <- compute_ci_mfmr(
          data$Y, 
          fpca_results$A_est, 
          data$instruments$Z,
          fpca_results, 
          config
        )
        MFMR_coverage <- calculate_coverage(MFMR_ci, data$beta_true_mat, config)
      }, error = function(e) {
        if (verbose) cat("    MFMR CI 计算失败:", e$message, "\n")
      })
      
      # ========== MPCMR: 使用结果中的Wald型CI ==========
      tryCatch({
        if (verbose) cat("    MPCMR GMM CI...\n")
        MPCMR_ci <- convert_mpcmr_ci(MPCMR_results, fpca_results, config)
        MPCMR_coverage <- calculate_coverage(MPCMR_ci, data$beta_true_mat, config)
      }, error = function(e) {
        if (verbose) cat("    MPCMR CI 计算失败:", e$message, "\n")
      })
      
      # ========== PACE+2SFRI: Bootstrap CI (使用共享样本) ==========
      tryCatch({
        if (verbose) cat("    PACE+2SFRI Bootstrap CI...\n")
        PACE_2SFRI_ci <- compute_ci_pace_2sfri(
          data$Y, 
          data$func_data, 
          data$obs_times_list, 
          config,
          data$instruments, 
          fpca_boot_shared,  # 【改进】使用同一批bootstrap样本
          fpca_results
        )
        PACE_2SFRI_coverage <- calculate_coverage(PACE_2SFRI_ci, data$beta_true_mat, config)
      }, error = function(e) {
        if (verbose) cat("    PACE+2SFRI CI 计算失败:", e$message, "\n")
      })
    }
    
    # 7. 返回结果
    return(list(
      success = TRUE,
      sim_id = sim_id,
      
      # MSE结果
      FA_MFMR = list(
        total_mse = lm_mse$total_mse,
        nonzero_mse = lm_mse$nonzero_mse,
        zero_mse = lm_mse$zero_mse,
        mse_by_group = lm_mse$mse_vec,
        coverage = FA_MFMR_coverage
      ),
      
      MFMR = list(
        total_mse = MFMR_mse$total_mse,
        nonzero_mse = MFMR_mse$nonzero_mse,
        zero_mse = MFMR_mse$zero_mse,
        mse_by_group = MFMR_mse$mse_vec,
        coverage = MFMR_coverage
      ),
      
      MPCMR = list(
        total_mse = MPCMR_mse$total_mse,
        nonzero_mse = MPCMR_mse$nonzero_mse,
        zero_mse = MPCMR_mse$zero_mse,
        mse_by_group = MPCMR_mse$mse_vec,
        coverage = MPCMR_coverage
      ),
      
      PACE_2SFRI = list(
        total_mse = PACE_2SFRI_mse$total_mse,
        nonzero_mse = PACE_2SFRI_mse$nonzero_mse,
        zero_mse = PACE_2SFRI_mse$zero_mse,
        mse_by_group = PACE_2SFRI_mse$mse_vec,
        coverage = PACE_2SFRI_coverage
      ),
      
      # 因子选择
      K_selected = factor_results$K_selected,
      K_true = config$K_true,
      
      # F统计量（FA-MFMR）
      F_stat_mean = if(FA_MFMR_results$success) mean(FA_MFMR_results$F_stat) else NA
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      sim_id = sim_id,
      error = as.character(e$message)
    ))
  })
}


#' 计算Coverage（真实值是否落在CI内）
calculate_coverage <- function(ci_results, beta_true_mat, config, 
                               key_times = NULL, alpha = 0.05) {
  
  if (is.null(key_times)) {
    key_times <- config$key_times
  }
  
  G <- config$G
  coverage_by_group <- numeric(G)
  ci_length_by_group <- numeric(G)
  
  for (g in 1:G) {
    if (!is.null(ci_results$beta_ci_functions[[g]])) {
      
      ci_info <- ci_results$beta_ci_functions[[g]]
      time_points <- ci_info$time
      ci_low <- ci_info$ci_low
      ci_up <- ci_info$ci_up
      
      # 在关键时间点评估
      eval_idx <- sapply(key_times, function(tt) which.min(abs(time_points - tt)))
      eval_idx <- eval_idx[eval_idx <= length(time_points) & eval_idx > 0]
      eval_idx <- unique(eval_idx)
      
      if (length(eval_idx) > 0) {
        # 真实值在这些时间点的值
        beta_true_interp <- approx(config$t_grid, beta_true_mat[, g], 
                                   xout = time_points[eval_idx], rule = 2)$y
        
        # 检查是否在CI内
        in_ci <- (beta_true_interp >= ci_low[eval_idx]) & 
          (beta_true_interp <= ci_up[eval_idx])
        coverage_by_group[g] <- mean(in_ci, na.rm = TRUE)
        
        # CI平均宽度
        ci_length_by_group[g] <- mean(ci_up[eval_idx] - ci_low[eval_idx], na.rm = TRUE)
      } else {
        coverage_by_group[g] <- NA
        ci_length_by_group[g] <- NA
      }
    } else {
      coverage_by_group[g] <- NA
      ci_length_by_group[g] <- NA
    }
  }
  
  # 计算总体统计量
  coverage_overall <- mean(coverage_by_group[1:4], na.rm = TRUE)  # 非零函数
  coverage_nonzero <- mean(coverage_by_group[1:4], na.rm = TRUE)
  coverage_zero <- mean(coverage_by_group[5:G], na.rm = TRUE)
  
  ci_length_overall <- mean(ci_length_by_group, na.rm = TRUE)
  ci_length_nonzero <- mean(ci_length_by_group[1:4], na.rm = TRUE)
  ci_length_zero <- mean(ci_length_by_group[5:G], na.rm = TRUE)
  
  return(list(
    coverage_by_group = coverage_by_group,
    coverage_overall = coverage_overall,
    coverage_nonzero = coverage_nonzero,
    coverage_zero = coverage_zero,
    ci_length_by_group = ci_length_by_group,
    ci_length_overall = ci_length_overall,
    ci_length_nonzero = ci_length_nonzero,
    ci_length_zero = ci_length_zero
  ))
}

# ===================================================================
# -----------多次重复仿真主函数-----------
# ===================================================================

run_repeated_simulations <- function(n_sim = 100, config = config, 
                                     parallel = TRUE, n_cores = NULL,
                                     verbose = TRUE, compute_ci = TRUE,
                                     n_bootstrap = 50) {
  
  cat("\n", rep("=", 70), "\n")
  cat(sprintf("   开始运行 %d 次重复仿真\n", n_sim))
  if (compute_ci) {
    cat(sprintf("   每次仿真包含 %d 次Bootstrap (计算置信区间)\n", n_bootstrap))
  }
  cat(rep("=", 70), "\n\n")
  
  start_time <- Sys.time()
  
  if (parallel) {
    # 并行计算
    if (is.null(n_cores)) {
      n_cores <- min(detectCores() - 1, n_sim)
    }
    
    cat(sprintf("使用并行计算，核心数: %d\n", n_cores))
    
    library(parallel)
    cl <- makeCluster(n_cores)
    
    # 导出必要的函数和对象到集群
    clusterExport(cl, varlist = ls(.GlobalEnv), envir = .GlobalEnv)
    
    # 加载必要的包
    clusterEvalQ(cl, {
      library(MASS)
      library(fda)
      library(fdapace)
      library(glmnet)
      library(AER)
      library(ncvreg)
    })
    
    # 运行仿真
    results_list <- parLapply(cl, 1:n_sim, function(i) {
      run_single_simulation(i, config, verbose = FALSE, 
                            compute_ci = compute_ci, n_bootstrap = n_bootstrap)
    })
    
    stopCluster(cl)
    
  } else {
    # 顺序计算
    cat("使用顺序计算\n")
    results_list <- vector("list", n_sim)
    
    for (i in 1:n_sim) {
      if (verbose && i %% 10 == 0) {
        cat(sprintf("完成: %d/%d (%.1f%%)\n", i, n_sim, 100*i/n_sim))
      }
      results_list[[i]] <- run_single_simulation(i, config, verbose = FALSE,
                                                 compute_ci = compute_ci,
                                                 n_bootstrap = n_bootstrap)
    }
  }
  
  end_time <- Sys.time()
  time_elapsed <- difftime(end_time, start_time, units = "mins")
  
  cat(sprintf("\n完成 %d 次仿真，用时: %.2f 分钟\n", n_sim, time_elapsed))
  
  # 汇总结果
  cat("\n汇总结果...\n")
  summary_results <- summarize_simulation_results(results_list, config, compute_ci)
  
  return(list(
    results_list = results_list,
    summary = summary_results,
    n_sim = n_sim,
    time_elapsed = time_elapsed,
    compute_ci = compute_ci
  ))
}


#' 汇总多次仿真结果
summarize_simulation_results <- function(results_list, config, compute_ci = TRUE) {
  
  n_sim <- length(results_list)
  
  # 统计成功率
  success_count <- sum(sapply(results_list, function(x) x$success))
  success_rate <- success_count / n_sim
  
  cat(sprintf("成功率: %.1f%% (%d/%d)\n", 100*success_rate, success_count, n_sim))
  
  # 只分析成功的仿真
  success_results <- results_list[sapply(results_list, function(x) x$success)]
  n_success <- length(success_results)
  
  if (n_success == 0) {
    cat("没有成功的仿真！\n")
    return(NULL)
  }
  
  # 提取MSE结果
  extract_mse <- function(method_name) {
    total_mse <- sapply(success_results, function(x) x[[method_name]]$total_mse)
    nonzero_mse <- sapply(success_results, function(x) x[[method_name]]$nonzero_mse)
    zero_mse <- sapply(success_results, function(x) x[[method_name]]$zero_mse)
    
    # 按组别的MSE
    mse_by_group <- t(sapply(success_results, function(x) x[[method_name]]$mse_by_group))
    
    list(
      total_mse = total_mse,
      nonzero_mse = nonzero_mse,
      zero_mse = zero_mse,
      mse_by_group = mse_by_group,
      
      # 统计量
      mean_total = mean(total_mse, na.rm = TRUE),
      sd_total = sd(total_mse, na.rm = TRUE),
      median_total = median(total_mse, na.rm = TRUE),
      
      mean_nonzero = mean(nonzero_mse, na.rm = TRUE),
      sd_nonzero = sd(nonzero_mse, na.rm = TRUE),
      
      mean_zero = mean(zero_mse, na.rm = TRUE),
      sd_zero = sd(zero_mse, na.rm = TRUE),
      
      mean_by_group = colMeans(mse_by_group, na.rm = TRUE),
      sd_by_group = apply(mse_by_group, 2, sd, na.rm = TRUE)
    )
  }
  
  # 提取Coverage和CI Length结果
  extract_coverage <- function(method_name) {
    if (!compute_ci) return(NULL)
    
    # 提取所有仿真的coverage结果
    coverage_overall_vec <- sapply(success_results, function(x) {
      if (!is.null(x[[method_name]]$coverage)) {
        x[[method_name]]$coverage$coverage_overall
      } else {
        NA
      }
    })
    
    coverage_nonzero_vec <- sapply(success_results, function(x) {
      if (!is.null(x[[method_name]]$coverage)) {
        x[[method_name]]$coverage$coverage_nonzero
      } else {
        NA
      }
    })
    
    coverage_zero_vec <- sapply(success_results, function(x) {
      if (!is.null(x[[method_name]]$coverage)) {
        x[[method_name]]$coverage$coverage_zero
      } else {
        NA
      }
    })
    
    # CI Length
    ci_length_overall_vec <- sapply(success_results, function(x) {
      if (!is.null(x[[method_name]]$coverage)) {
        x[[method_name]]$coverage$ci_length_overall
      } else {
        NA
      }
    })
    
    ci_length_nonzero_vec <- sapply(success_results, function(x) {
      if (!is.null(x[[method_name]]$coverage)) {
        x[[method_name]]$coverage$ci_length_nonzero
      } else {
        NA
      }
    })
    
    ci_length_zero_vec <- sapply(success_results, function(x) {
      if (!is.null(x[[method_name]]$coverage)) {
        x[[method_name]]$coverage$ci_length_zero
      } else {
        NA
      }
    })
    
    # 按组别的coverage和length（跨时间点的平均）
    coverage_by_group_mat <- t(sapply(success_results, function(x) {
      if (!is.null(x[[method_name]]$coverage)) {
        x[[method_name]]$coverage$coverage_by_group
      } else {
        rep(NA, config$G)
      }
    }))
    
    ci_length_by_group_mat <- t(sapply(success_results, function(x) {
      if (!is.null(x[[method_name]]$coverage)) {
        x[[method_name]]$coverage$ci_length_by_group
      } else {
        rep(NA, config$G)
      }
    }))
    
    list(
      # Coverage统计
      coverage_overall = coverage_overall_vec,
      mean_coverage_overall = mean(coverage_overall_vec, na.rm = TRUE),
      sd_coverage_overall = sd(coverage_overall_vec, na.rm = TRUE),
      
      coverage_nonzero = coverage_nonzero_vec,
      mean_coverage_nonzero = mean(coverage_nonzero_vec, na.rm = TRUE),
      sd_coverage_nonzero = sd(coverage_nonzero_vec, na.rm = TRUE),
      
      coverage_zero = coverage_zero_vec,
      mean_coverage_zero = mean(coverage_zero_vec, na.rm = TRUE),
      sd_coverage_zero = sd(coverage_zero_vec, na.rm = TRUE),
      
      coverage_by_group_mat = coverage_by_group_mat,
      mean_coverage_by_group = colMeans(coverage_by_group_mat, na.rm = TRUE),
      sd_coverage_by_group = apply(coverage_by_group_mat, 2, sd, na.rm = TRUE),
      
      # CI Length统计
      ci_length_overall = ci_length_overall_vec,
      mean_ci_length_overall = mean(ci_length_overall_vec, na.rm = TRUE),
      sd_ci_length_overall = sd(ci_length_overall_vec, na.rm = TRUE),
      
      ci_length_nonzero = ci_length_nonzero_vec,
      mean_ci_length_nonzero = mean(ci_length_nonzero_vec, na.rm = TRUE),
      sd_ci_length_nonzero = sd(ci_length_nonzero_vec, na.rm = TRUE),
      
      ci_length_zero = ci_length_zero_vec,
      mean_ci_length_zero = mean(ci_length_zero_vec, na.rm = TRUE),
      sd_ci_length_zero = sd(ci_length_zero_vec, na.rm = TRUE),
      
      ci_length_by_group_mat = ci_length_by_group_mat,
      mean_ci_length_by_group = colMeans(ci_length_by_group_mat, na.rm = TRUE),
      sd_ci_length_by_group = apply(ci_length_by_group_mat, 2, sd, na.rm = TRUE)
    )
  }
  
  FA_MFMR_summary <- extract_mse("FA_MFMR")
  MFMR_summary <- extract_mse("MFMR")
  MPCMR_summary <- extract_mse("MPCMR")
  PACE_2SFRI_summary <- extract_mse("PACE_2SFRI")
  
  FA_MFMR_coverage <- extract_coverage("FA_MFMR")
  MFMR_coverage <- extract_coverage("MFMR")
  MPCMR_coverage <- extract_coverage("MPCMR")
  PACE_2SFRI_coverage <- extract_coverage("PACE_2SFRI")
  
  # 因子数选择统计
  K_selected <- sapply(success_results, function(x) x$K_selected)
  K_true <- config$K_true
  
  # F统计量
  F_stats <- sapply(success_results, function(x) x$F_stat_mean)
  
  summary_list <- list(
    n_sim = n_sim,
    n_success = n_success,
    success_rate = success_rate,
    compute_ci = compute_ci,
    
    FA_MFMR = FA_MFMR_summary,
    MFMR = MFMR_summary,
    MPCMR = MPCMR_summary,
    PACE_2SFRI = PACE_2SFRI_summary,
    
    FA_MFMR_coverage = FA_MFMR_coverage,
    MFMR_coverage = MFMR_coverage,
    MPCMR_coverage = MPCMR_coverage,
    PACE_2SFRI_coverage = PACE_2SFRI_coverage,
    
    K_selected = list(
      values = K_selected,
      mean = mean(K_selected, na.rm = TRUE),
      sd = sd(K_selected, na.rm = TRUE),
      correct_rate = mean(K_selected == K_true, na.rm = TRUE)
    ),
    
    F_stat = list(
      values = F_stats,
      mean = mean(F_stats, na.rm = TRUE),
      sd = sd(F_stats, na.rm = TRUE)
    )
  )
  
  return(summary_list)
}


#' 打印汇总结果
print_summary <- function(summary) {
  
  cat("\n", rep("=", 70), "\n")
  cat("                   多次仿真结果汇总\n")
  cat(rep("=", 70), "\n")
  
  cat(sprintf("\n总仿真次数: %d\n", summary$n_sim))
  cat(sprintf("成功次数: %d (%.1f%%)\n", summary$n_success, 100*summary$success_rate))
  
  # 因子数选择
  cat(sprintf("\n因子数选择:\n"))
  cat(sprintf("  真实值: %d\n", config$K_true))
  cat(sprintf("  平均估计值: %.2f (±%.2f)\n", 
              summary$K_selected$mean, summary$K_selected$sd))
  cat(sprintf("  正确率: %.1f%%\n", 100*summary$K_selected$correct_rate))
  
  # F统计量
  cat(sprintf("\nF统计量 (FA-MFMR):\n"))
  cat(sprintf("  平均值: %.2f (±%.2f)\n", 
              summary$F_stat$mean, summary$F_stat$sd))
  
  # MSE对比
  cat("\n", rep("-", 70), "\n")
  cat("MSE 对比 (Mean ± SD)\n")
  cat(rep("-", 70), "\n")
  
  cat("\n【总体MSE】\n")
  cat(sprintf("  %-15s: %.6f ± %.6f\n", "FA-MFMR", 
              summary$FA_MFMR$mean_total, summary$FA_MFMR$sd_total))
  cat(sprintf("  %-15s: %.6f ± %.6f\n", "MFMR", 
              summary$MFMR$mean_total, summary$MFMR$sd_total))
  cat(sprintf("  %-15s: %.6f ± %.6f\n", "MPCMR", 
              summary$MPCMR$mean_total, summary$MPCMR$sd_total))
  cat(sprintf("  %-15s: %.6f ± %.6f\n", "PACE+2SFRI", 
              summary$PACE_2SFRI$mean_total, summary$PACE_2SFRI$sd_total))
  
  cat("\n【非零函数MSE (Beta1-4)】\n")
  cat(sprintf("  %-15s: %.6f ± %.6f\n", "FA-MFMR", 
              summary$FA_MFMR$mean_nonzero, summary$FA_MFMR$sd_nonzero))
  cat(sprintf("  %-15s: %.6f ± %.6f\n", "MFMR", 
              summary$MFMR$mean_nonzero, summary$MFMR$sd_nonzero))
  cat(sprintf("  %-15s: %.6f ± %.6f\n", "MPCMR", 
              summary$MPCMR$mean_nonzero, summary$MPCMR$sd_nonzero))
  cat(sprintf("  %-15s: %.6f ± %.6f\n", "PACE+2SFRI", 
              summary$PACE_2SFRI$mean_nonzero, summary$PACE_2SFRI$sd_nonzero))
  
  cat("\n【零函数MSE (Beta5-20)】\n")
  cat(sprintf("  %-15s: %.6f ± %.6f\n", "FA-MFMR", 
              summary$FA_MFMR$mean_zero, summary$FA_MFMR$sd_zero))
  cat(sprintf("  %-15s: %.6f ± %.6f\n", "MFMR", 
              summary$MFMR$mean_zero, summary$MFMR$sd_zero))
  cat(sprintf("  %-15s: %.6f ± %.6f\n", "MPCMR", 
              summary$MPCMR$mean_zero, summary$MPCMR$sd_zero))
  cat(sprintf("  %-15s: %.6f ± %.6f\n", "PACE+2SFRI", 
              summary$PACE_2SFRI$mean_zero, summary$PACE_2SFRI$sd_zero))
  
  # ========== 【修改】分组MSE - 显示前4个非零函数 ==========
  cat("\n【各Beta函数MSE (Mean)】\n")
  group_names <- c("Beta1(线性)", "Beta2(二次)", "Beta3(三角)", "Beta4(高斯)")
  
  for (g in 1:4) {
    cat(sprintf("\n  %s:\n", group_names[g]))
    cat(sprintf("    FA-MFMR    : %.6f\n", summary$FA_MFMR$mean_by_group[g]))
    cat(sprintf("    MFMR       : %.6f\n", summary$MFMR$mean_by_group[g]))
    cat(sprintf("    MPCMR      : %.6f\n", summary$MPCMR$mean_by_group[g]))
    cat(sprintf("    PACE+2SFRI : %.6f\n", summary$PACE_2SFRI$mean_by_group[g]))
  }
  
  # 显示零函数的平均MSE
  cat(sprintf("\n  零函数(Beta5-20)平均:\n"))
  cat(sprintf("    FA-MFMR    : %.6f\n", mean(summary$FA_MFMR$mean_by_group[5:config$G])))
  cat(sprintf("    MFMR       : %.6f\n", mean(summary$MFMR$mean_by_group[5:config$G])))
  cat(sprintf("    MPCMR      : %.6f\n", mean(summary$MPCMR$mean_by_group[5:config$G])))
  cat(sprintf("    PACE+2SFRI : %.6f\n", mean(summary$PACE_2SFRI$mean_by_group[5:config$G])))
  
  # ========== 【新增】Coverage和CI Length对比 ==========
  if (summary$compute_ci) {
    cat("\n", rep("-", 70), "\n")
    cat("Coverage 对比 (Mean ± SD, 名义水平: 95%)\n")
    cat(rep("-", 70), "\n")
    
    cat("\n【总体Coverage】\n")
    cat(sprintf("  %-15s: %.3f ± %.3f\n", "FA-MFMR", 
                summary$FA_MFMR_coverage$mean_coverage_overall, 
                summary$FA_MFMR_coverage$sd_coverage_overall))
    cat(sprintf("  %-15s: %.3f ± %.3f\n", "MFMR", 
                summary$MFMR_coverage$mean_coverage_overall,
                summary$MFMR_coverage$sd_coverage_overall))
    cat(sprintf("  %-15s: %.3f ± %.3f\n", "MPCMR", 
                summary$MPCMR_coverage$mean_coverage_overall,
                summary$MPCMR_coverage$sd_coverage_overall))
    cat(sprintf("  %-15s: %.3f ± %.3f\n", "PACE+2SFRI", 
                summary$PACE_2SFRI_coverage$mean_coverage_overall,
                summary$PACE_2SFRI_coverage$sd_coverage_overall))
    
    cat("\n【非零函数Coverage (Beta1-4)】\n")
    cat(sprintf("  %-15s: %.3f ± %.3f\n", "FA-MFMR", 
                summary$FA_MFMR_coverage$mean_coverage_nonzero,
                summary$FA_MFMR_coverage$sd_coverage_nonzero))
    cat(sprintf("  %-15s: %.3f ± %.3f\n", "MFMR", 
                summary$MFMR_coverage$mean_coverage_nonzero,
                summary$MFMR_coverage$sd_coverage_nonzero))
    cat(sprintf("  %-15s: %.3f ± %.3f\n", "MPCMR", 
                summary$MPCMR_coverage$mean_coverage_nonzero,
                summary$MPCMR_coverage$sd_coverage_nonzero))
    cat(sprintf("  %-15s: %.3f ± %.3f\n", "PACE+2SFRI", 
                summary$PACE_2SFRI_coverage$mean_coverage_nonzero,
                summary$PACE_2SFRI_coverage$sd_coverage_nonzero))
    
    cat("\n【零函数Coverage (Beta5-20)】\n")
    cat(sprintf("  %-15s: %.3f ± %.3f\n", "FA-MFMR", 
                summary$FA_MFMR_coverage$mean_coverage_zero,
                summary$FA_MFMR_coverage$sd_coverage_zero))
    cat(sprintf("  %-15s: %.3f ± %.3f\n", "MFMR", 
                summary$MFMR_coverage$mean_coverage_zero,
                summary$MFMR_coverage$sd_coverage_zero))
    cat(sprintf("  %-15s: %.3f ± %.3f\n", "MPCMR", 
                summary$MPCMR_coverage$mean_coverage_zero,
                summary$MPCMR_coverage$sd_coverage_zero))
    cat(sprintf("  %-15s: %.3f ± %.3f\n", "PACE+2SFRI", 
                summary$PACE_2SFRI_coverage$mean_coverage_zero,
                summary$PACE_2SFRI_coverage$sd_coverage_zero))
    
    cat("\n", rep("-", 70), "\n")
    cat("CI Length 对比 (Mean ± SD)\n")
    cat(rep("-", 70), "\n")
    
    cat("\n【总体CI Length】\n")
    cat(sprintf("  %-15s: %.4f ± %.4f\n", "FA-MFMR", 
                summary$FA_MFMR_coverage$mean_ci_length_overall,
                summary$FA_MFMR_coverage$sd_ci_length_overall))
    cat(sprintf("  %-15s: %.4f ± %.4f\n", "MFMR", 
                summary$MFMR_coverage$mean_ci_length_overall,
                summary$MFMR_coverage$sd_ci_length_overall))
    cat(sprintf("  %-15s: %.4f ± %.4f\n", "MPCMR", 
                summary$MPCMR_coverage$mean_ci_length_overall,
                summary$MPCMR_coverage$sd_ci_length_overall))
    cat(sprintf("  %-15s: %.4f ± %.4f\n", "PACE+2SFRI", 
                summary$PACE_2SFRI_coverage$mean_ci_length_overall,
                summary$PACE_2SFRI_coverage$sd_ci_length_overall))
    
    cat("\n【非零函数CI Length (Beta1-4)】\n")
    cat(sprintf("  %-15s: %.4f ± %.4f\n", "FA-MFMR", 
                summary$FA_MFMR_coverage$mean_ci_length_nonzero,
                summary$FA_MFMR_coverage$sd_ci_length_nonzero))
    cat(sprintf("  %-15s: %.4f ± %.4f\n", "MFMR", 
                summary$MFMR_coverage$mean_ci_length_nonzero,
                summary$MFMR_coverage$sd_ci_length_nonzero))
    cat(sprintf("  %-15s: %.4f ± %.4f\n", "MPCMR", 
                summary$MPCMR_coverage$mean_ci_length_nonzero,
                summary$MPCMR_coverage$sd_ci_length_nonzero))
    cat(sprintf("  %-15s: %.4f ± %.4f\n", "PACE+2SFRI", 
                summary$PACE_2SFRI_coverage$mean_ci_length_nonzero,
                summary$PACE_2SFRI_coverage$sd_ci_length_nonzero))
    
    cat("\n【零函数CI Length (Beta5-20)】\n")
    cat(sprintf("  %-15s: %.4f ± %.4f\n", "FA-MFMR", 
                summary$FA_MFMR_coverage$mean_ci_length_zero,
                summary$FA_MFMR_coverage$sd_ci_length_zero))
    cat(sprintf("  %-15s: %.4f ± %.4f\n", "MFMR", 
                summary$MFMR_coverage$mean_ci_length_zero,
                summary$MFMR_coverage$sd_ci_length_zero))
    cat(sprintf("  %-15s: %.4f ± %.4f\n", "MPCMR", 
                summary$MPCMR_coverage$mean_ci_length_zero,
                summary$MPCMR_coverage$sd_ci_length_zero))
    cat(sprintf("  %-15s: %.4f ± %.4f\n", "PACE+2SFRI", 
                summary$PACE_2SFRI_coverage$mean_ci_length_zero,
                summary$PACE_2SFRI_coverage$sd_ci_length_zero))
    
    # ========== 【新增】按组别显示Coverage和CI Length ==========
    cat("\n", rep("-", 70), "\n")
    cat("各Beta函数详细结果 (Coverage | CI Length)\n")
    cat(rep("-", 70), "\n")
    
    group_names <- c("Beta1(线性)", "Beta2(二次)", "Beta3(三角)", "Beta4(高斯)")
    
    for (g in 1:4) {
      cat(sprintf("\n【%s】\n", group_names[g]))
      
      cat("  Coverage:\n")
      cat(sprintf("    FA-MFMR    : %.3f\n", 
                  summary$FA_MFMR_coverage$mean_coverage_by_group[g]))
      cat(sprintf("    MFMR       : %.3f\n", 
                  summary$MFMR_coverage$mean_coverage_by_group[g]))
      cat(sprintf("    MPCMR      : %.3f\n", 
                  summary$MPCMR_coverage$mean_coverage_by_group[g]))
      cat(sprintf("    PACE+2SFRI : %.3f\n", 
                  summary$PACE_2SFRI_coverage$mean_coverage_by_group[g]))
      
      cat("  CI Length:\n")
      cat(sprintf("    FA-MFMR    : %.4f\n", 
                  summary$FA_MFMR_coverage$mean_ci_length_by_group[g]))
      cat(sprintf("    MFMR       : %.4f\n", 
                  summary$MFMR_coverage$mean_ci_length_by_group[g]))
      cat(sprintf("    MPCMR      : %.4f\n", 
                  summary$MPCMR_coverage$mean_ci_length_by_group[g]))
      cat(sprintf("    PACE+2SFRI : %.4f\n", 
                  summary$PACE_2SFRI_coverage$mean_ci_length_by_group[g]))
    }
    
    # 显示零函数的平均Coverage和CI Length
    cat(sprintf("\n【零函数(Beta5-20)平均】\n"))
    cat("  Coverage:\n")
    cat(sprintf("    FA-MFMR    : %.3f\n", 
                mean(summary$FA_MFMR_coverage$mean_coverage_by_group[5:config$G])))
    cat(sprintf("    MFMR       : %.3f\n", 
                mean(summary$MFMR_coverage$mean_coverage_by_group[5:config$G])))
    cat(sprintf("    MPCMR      : %.3f\n", 
                mean(summary$MPCMR_coverage$mean_coverage_by_group[5:config$G])))
    cat(sprintf("    PACE+2SFRI : %.3f\n", 
                mean(summary$PACE_2SFRI_coverage$mean_coverage_by_group[5:config$G])))
    
    cat("  CI Length:\n")
    cat(sprintf("    FA-MFMR    : %.4f\n", 
                mean(summary$FA_MFMR_coverage$mean_ci_length_by_group[5:config$G])))
    cat(sprintf("    MFMR       : %.4f\n", 
                mean(summary$MFMR_coverage$mean_ci_length_by_group[5:config$G])))
    cat(sprintf("    MPCMR      : %.4f\n", 
                mean(summary$MPCMR_coverage$mean_ci_length_by_group[5:config$G])))
    cat(sprintf("    PACE+2SFRI : %.4f\n", 
                mean(summary$PACE_2SFRI_coverage$mean_ci_length_by_group[5:config$G])))
  }
  
  # 最佳方法
  methods <- c("FA-MFMR", "MFMR", "MPCMR", "PACE+2SFRI")
  total_mse_means <- c(summary$FA_MFMR$mean_total, 
                       summary$MFMR$mean_total,
                       summary$MPCMR$mean_total,
                       summary$PACE_2SFRI$mean_total)
  
  best_idx <- which.min(total_mse_means)
  
  cat("\n", rep("-", 70), "\n")
  cat(sprintf("【最佳方法】(基于平均总体MSE): %s (%.6f)\n", 
              methods[best_idx], total_mse_means[best_idx]))
  cat(rep("-", 70), "\n")
}


#' 绘制MSE箱线图
plot_mse_boxplot <- function(summary) {
  
  # 准备数据
  mse_data <- data.frame(
    Method = rep(c("FA-MFMR", "MFMR", "MPCMR", "PACE+2SFRI"), 
                 each = summary$n_success),
    Total_MSE = c(summary$FA_MFMR$total_mse,
                  summary$MFMR$total_mse,
                  summary$MPCMR$total_mse,
                  summary$PACE_2SFRI$total_mse),
    Nonzero_MSE = c(summary$FA_MFMR$nonzero_mse,
                    summary$MFMR$nonzero_mse,
                    summary$MPCMR$nonzero_mse,
                    summary$PACE_2SFRI$nonzero_mse),
    Zero_MSE = c(summary$FA_MFMR$zero_mse,
                 summary$MFMR$zero_mse,
                 summary$MPCMR$zero_mse,
                 summary$PACE_2SFRI$zero_mse)
  )
  
  # 三个箱线图
  par(mfrow = c(1, 3), mar = c(8, 4, 3, 2))
  
  # 总体MSE
  boxplot(Total_MSE ~ Method, data = mse_data,
          main = "总体MSE",
          ylab = "MSE",
          las = 2,
          col = c("lightblue", "lightgreen", "lightyellow", "lightpink"))
  
  # 非零函数MSE
  boxplot(Nonzero_MSE ~ Method, data = mse_data,
          main = "非零函数MSE (Beta1-4)",
          ylab = "MSE",
          las = 2,
          col = c("lightblue", "lightgreen", "lightyellow", "lightpink"))
  
  # 零函数MSE
  boxplot(Zero_MSE ~ Method, data = mse_data,
          main = "零函数MSE (Beta5-20)",
          ylab = "MSE",
          las = 2,
          col = c("lightblue", "lightgreen", "lightyellow", "lightpink"))
  
  par(mfrow = c(1, 1))
}


#' 保存结果到文件
#' 保存结果到文件（在原有文件中添加Coverage和CI Length）
save_simulation_results <- function(repeated_results, 
                                    output_dir = "simulation_results",
                                    config = NULL) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("✓ 创建输出目录: %s\n", output_dir))
  }
  
  results_list <- repeated_results$results_list
  summary <- repeated_results$summary
  n_sim <- repeated_results$n_sim
  compute_ci <- repeated_results$compute_ci
  
  # ========== 1. 提取逐次仿真的逐组数据（MSE + Coverage + CI Length） ==========
  cat("处理逐次仿真结果...\n")
  
  detailed_results_list <- list()
  
  for (sim_id in 1:n_sim) {
    result <- results_list[[sim_id]]
    
    if (!result$success) {
      next
    }
    
    # 提取四种方法的结果
    for (method_name in c("FA_MFMR", "MFMR", "MPCMR", "PACE_2SFRI")) {
      method_data <- result[[method_name]]
      
      # 按组别保存
      for (g in 1:length(method_data$mse_by_group)) {
        row_data <- data.frame(
          sim_id = sim_id,
          method = method_name,
          group = g,
          mse = method_data$mse_by_group[g],
          total_mse = method_data$total_mse,
          nonzero_mse = method_data$nonzero_mse,
          zero_mse = method_data$zero_mse,
          stringsAsFactors = FALSE
        )
        
        # 如果计算了CI，添加Coverage和CI Length列
        if (compute_ci && !is.null(method_data$coverage)) {
          coverage_data <- method_data$coverage
          row_data$coverage <- coverage_data$coverage_by_group[g]
          row_data$ci_length <- coverage_data$ci_length_by_group[g]
          row_data$coverage_overall <- coverage_data$coverage_overall
          row_data$coverage_nonzero <- coverage_data$coverage_nonzero
          row_data$coverage_zero <- coverage_data$coverage_zero
          row_data$ci_length_overall <- coverage_data$ci_length_overall
          row_data$ci_length_nonzero <- coverage_data$ci_length_nonzero
          row_data$ci_length_zero <- coverage_data$ci_length_zero
        } else {
          row_data$coverage <- NA
          row_data$ci_length <- NA
          row_data$coverage_overall <- NA
          row_data$coverage_nonzero <- NA
          row_data$coverage_zero <- NA
          row_data$ci_length_overall <- NA
          row_data$ci_length_nonzero <- NA
          row_data$ci_length_zero <- NA
        }
        
        detailed_results_list[[length(detailed_results_list) + 1]] <- row_data
      }
    }
  }
  
  # 转换为数据框
  detailed_df <- do.call(rbind, detailed_results_list)
  rownames(detailed_df) <- NULL
  
  # 保存为txt
  txt_file <- file.path(output_dir, "detailed_simulation_results.txt")
  write.table(detailed_df, file = txt_file, row.names = FALSE, 
              sep = "\t", quote = FALSE)
  cat(sprintf("✓ 逐次逐组详细结果: %s (共%d行)\n", txt_file, nrow(detailed_df)))
  
  # ========== 2. 汇总表（按方法汇总 - MSE + Coverage + CI Length） ==========
  cat("生成方法汇总表...\n")
  
  methods <- c("FA_MFMR", "MFMR", "MPCMR", "PACE_2SFRI")
  summary_table <- list()
  
  for (method in methods) {
    method_filter <- detailed_df$method == method
    method_data <- detailed_df[method_filter, ]
    
    summary_row <- data.frame(
      Method = method,
      Total_MSE_Mean = mean(method_data$total_mse, na.rm = TRUE),
      Total_MSE_SD = sd(method_data$total_mse, na.rm = TRUE),
      Total_MSE_Median = median(method_data$total_mse, na.rm = TRUE),
      Total_MSE_Min = min(method_data$total_mse, na.rm = TRUE),
      Total_MSE_Max = max(method_data$total_mse, na.rm = TRUE),
      Nonzero_MSE_Mean = mean(method_data$nonzero_mse, na.rm = TRUE),
      Nonzero_MSE_SD = sd(method_data$nonzero_mse, na.rm = TRUE),
      Zero_MSE_Mean = mean(method_data$zero_mse, na.rm = TRUE),
      Zero_MSE_SD = sd(method_data$zero_mse, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    # 如果计算了CI，添加Coverage和CI Length列
    if (compute_ci) {
      summary_row$Coverage_Overall_Mean = mean(method_data$coverage_overall, na.rm = TRUE)
      summary_row$Coverage_Overall_SD = sd(method_data$coverage_overall, na.rm = TRUE)
      summary_row$Coverage_Nonzero_Mean = mean(method_data$coverage_nonzero, na.rm = TRUE)
      summary_row$Coverage_Nonzero_SD = sd(method_data$coverage_nonzero, na.rm = TRUE)
      summary_row$Coverage_Zero_Mean = mean(method_data$coverage_zero, na.rm = TRUE)
      summary_row$Coverage_Zero_SD = sd(method_data$coverage_zero, na.rm = TRUE)
      summary_row$CI_Length_Overall_Mean = mean(method_data$ci_length_overall, na.rm = TRUE)
      summary_row$CI_Length_Overall_SD = sd(method_data$ci_length_overall, na.rm = TRUE)
      summary_row$CI_Length_Nonzero_Mean = mean(method_data$ci_length_nonzero, na.rm = TRUE)
      summary_row$CI_Length_Nonzero_SD = sd(method_data$ci_length_nonzero, na.rm = TRUE)
      summary_row$CI_Length_Zero_Mean = mean(method_data$ci_length_zero, na.rm = TRUE)
      summary_row$CI_Length_Zero_SD = sd(method_data$ci_length_zero, na.rm = TRUE)
    }
    
    summary_row$N_Simulations = length(unique(method_data$sim_id))
    
    summary_table[[length(summary_table) + 1]] <- summary_row
  }
  
  summary_df <- do.call(rbind, summary_table)
  rownames(summary_df) <- NULL
  
  summary_txt <- file.path(output_dir, "summary_by_method.txt")
  write.table(summary_df, file = summary_txt, row.names = FALSE, 
              sep = "\t", quote = FALSE)
  cat(sprintf("✓ 方法汇总表: %s\n", summary_txt))
  
  # ========== 3. 各组别详细对比表（MSE + Coverage + CI Length） ==========
  cat("生成各组别对比表...\n")
  
  group_summary_list <- list()
  n_groups <- max(detailed_df$group, na.rm = TRUE)
  
  for (g in 1:n_groups) {
    group_filter <- detailed_df$group == g
    group_data <- detailed_df[group_filter, ]
    
    for (method in methods) {
      method_filter <- group_data$method == method
      method_group_data <- group_data[method_filter, ]
      
      if (nrow(method_group_data) > 0) {
        row <- data.frame(
          Group = g,
          Method = method,
          MSE_Mean = mean(method_group_data$mse, na.rm = TRUE),
          MSE_SD = sd(method_group_data$mse, na.rm = TRUE),
          MSE_Min = min(method_group_data$mse, na.rm = TRUE),
          MSE_Max = max(method_group_data$mse, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
        
        # 如果计算了CI，添加Coverage和CI Length列
        if (compute_ci) {
          row$Coverage_Mean = mean(method_group_data$coverage, na.rm = TRUE)
          row$Coverage_SD = sd(method_group_data$coverage, na.rm = TRUE)
          row$CI_Length_Mean = mean(method_group_data$ci_length, na.rm = TRUE)
          row$CI_Length_SD = sd(method_group_data$ci_length, na.rm = TRUE)
        }
        
        row$N_Simulations = length(unique(method_group_data$sim_id))
        
        group_summary_list[[length(group_summary_list) + 1]] <- row
      }
    }
  }
  
  group_summary_df <- do.call(rbind, group_summary_list)
  rownames(group_summary_df) <- NULL
  
  group_txt <- file.path(output_dir, "summary_by_group.txt")
  write.table(group_summary_df, file = group_txt, row.names = FALSE, 
              sep = "\t", quote = FALSE)
  cat(sprintf("✓ 各组别对比表: %s\n", group_txt))
  
  # ========== 4. 生成文本报告 ==========
  cat("生成文本报告...\n")
  
  report_file <- file.path(output_dir, "simulation_report.txt")
  sink(report_file)
  
  line_sep <- paste(rep("=", 80), collapse = "")
  dash_sep <- paste(rep("-", 80), collapse = "")
  
  cat(line_sep, "\n")
  cat("FA-MFMR 模拟研究结果报告\n")
  cat(line_sep, "\n\n")
  
  cat(sprintf("总仿真次数: %d\n", n_sim))
  cat(sprintf("成功率: %.1f%% (%d/%d)\n", 
              100*summary$success_rate, summary$n_success, n_sim))
  cat(sprintf("是否计算置信区间: %s\n\n", ifelse(compute_ci, "是", "否")))
  
  # 因子数选择统计
  cat("【因子数选择统计】\n")
  cat(sprintf("真实K值: %d\n", config$K_true))
  cat(sprintf("平均估计K: %.2f (±%.2f)\n", 
              summary$K_selected$mean, summary$K_selected$sd))
  cat(sprintf("选择正确的比例: %.1f%%\n\n", 100*summary$K_selected$correct_rate))
  
  cat(line_sep, "\n")
  cat("【方法MSE对比】\n")
  cat(line_sep, "\n\n")
  
  # 打印汇总表
  for (i in 1:nrow(summary_df)) {
    cat(sprintf("方法: %s\n", summary_df$Method[i]))
    cat(sprintf("  总体MSE:     Mean = %.6f, SD = %.6f, Median = %.6f\n",
                summary_df$Total_MSE_Mean[i],
                summary_df$Total_MSE_SD[i],
                summary_df$Total_MSE_Median[i]))
    cat(sprintf("  非零函数MSE: Mean = %.6f, SD = %.6f\n",
                summary_df$Nonzero_MSE_Mean[i],
                summary_df$Nonzero_MSE_SD[i]))
    cat(sprintf("  零函数MSE:   Mean = %.6f, SD = %.6f\n",
                summary_df$Zero_MSE_Mean[i],
                summary_df$Zero_MSE_SD[i]))
    
    # 如果有Coverage和CI Length数据，也打印出来
    if (compute_ci) {
      cat("  Coverage (名义水平: 95%):\n")
      cat(sprintf("    总体:     %.3f ± %.3f\n",
                  summary_df$Coverage_Overall_Mean[i],
                  summary_df$Coverage_Overall_SD[i]))
      cat(sprintf("    非零函数: %.3f ± %.3f\n",
                  summary_df$Coverage_Nonzero_Mean[i],
                  summary_df$Coverage_Nonzero_SD[i]))
      cat(sprintf("    零函数:   %.3f ± %.3f\n",
                  summary_df$Coverage_Zero_Mean[i],
                  summary_df$Coverage_Zero_SD[i]))
      cat("  CI Length:\n")
      cat(sprintf("    总体:     %.4f ± %.4f\n",
                  summary_df$CI_Length_Overall_Mean[i],
                  summary_df$CI_Length_Overall_SD[i]))
      cat(sprintf("    非零函数: %.4f ± %.4f\n",
                  summary_df$CI_Length_Nonzero_Mean[i],
                  summary_df$CI_Length_Nonzero_SD[i]))
      cat(sprintf("    零函数:   %.4f ± %.4f\n",
                  summary_df$CI_Length_Zero_Mean[i],
                  summary_df$CI_Length_Zero_SD[i]))
    }
    
    cat(sprintf("  仿真次数: %d\n\n",
                summary_df$N_Simulations[i]))
  }
  
  cat("\n")
  
  # 各组别详细
  cat(line_sep, "\n")
  cat("【各组别详细对比 (前4个非零函数)】\n")
  cat(line_sep, "\n\n")
  
  group_names <- c("Beta1(线性)", "Beta2(二次)", "Beta3(三角)", "Beta4(高斯)")
  for (g in 1:4) {
    cat(sprintf("【%s】\n", group_names[g]))
    group_g <- group_summary_df[group_summary_df$Group == g, ]
    
    cat("  MSE:\n")
    for (j in 1:nrow(group_g)) {
      cat(sprintf("    %s: %.6f ± %.6f [%.6f, %.6f]\n",
                  group_g$Method[j],
                  group_g$MSE_Mean[j],
                  group_g$MSE_SD[j],
                  group_g$MSE_Min[j],
                  group_g$MSE_Max[j]))
    }
    
    if (compute_ci) {
      cat("  Coverage:\n")
      for (j in 1:nrow(group_g)) {
        cat(sprintf("    %s: %.3f ± %.3f\n",
                    group_g$Method[j],
                    group_g$Coverage_Mean[j],
                    group_g$Coverage_SD[j]))
      }
      cat("  CI Length:\n")
      for (j in 1:nrow(group_g)) {
        cat(sprintf("    %s: %.4f ± %.4f\n",
                    group_g$Method[j],
                    group_g$CI_Length_Mean[j],
                    group_g$CI_Length_SD[j]))
      }
    }
    cat("\n")
  }
  
  # 零函数平均统计
  cat("【零函数(Beta5-20)平均统计】\n")
  cat("  MSE:\n")
  for (method in methods) {
    zero_groups <- group_summary_df[group_summary_df$Group >= 5 & 
                                      group_summary_df$Method == method, ]
    cat(sprintf("    %s: %.6f ± %.6f\n",
                method,
                mean(zero_groups$MSE_Mean, na.rm = TRUE),
                mean(zero_groups$MSE_SD, na.rm = TRUE)))
  }
  
  if (compute_ci) {
    cat("  Coverage:\n")
    for (method in methods) {
      zero_groups <- group_summary_df[group_summary_df$Group >= 5 & 
                                        group_summary_df$Method == method, ]
      cat(sprintf("    %s: %.3f ± %.3f\n",
                  method,
                  mean(zero_groups$Coverage_Mean, na.rm = TRUE),
                  mean(zero_groups$Coverage_SD, na.rm = TRUE)))
    }
    cat("  CI Length:\n")
    for (method in methods) {
      zero_groups <- group_summary_df[group_summary_df$Group >= 5 & 
                                        group_summary_df$Method == method, ]
      cat(sprintf("    %s: %.4f ± %.4f\n",
                  method,
                  mean(zero_groups$CI_Length_Mean, na.rm = TRUE),
                  mean(zero_groups$CI_Length_SD, na.rm = TRUE)))
    }
  }
  cat("\n")
  
  # 方法排名
  cat(line_sep, "\n")
  cat("【方法排名 (基于Total_MSE_Mean)】\n")
  cat(line_sep, "\n\n")
  
  ranking <- summary_df[order(summary_df$Total_MSE_Mean), ]
  
  for (i in 1:nrow(ranking)) {
    cat(sprintf("第%d名: %s\n", i, ranking$Method[i]))
    cat(sprintf("       Total_MSE = %.6f ± %.6f\n",
                ranking$Total_MSE_Mean[i],
                ranking$Total_MSE_SD[i]))
    cat(sprintf("       Nonzero_MSE = %.6f ± %.6f\n",
                ranking$Nonzero_MSE_Mean[i],
                ranking$Nonzero_MSE_SD[i]))
    cat(sprintf("       Zero_MSE = %.6f ± %.6f\n",
                ranking$Zero_MSE_Mean[i],
                ranking$Zero_MSE_SD[i]))
    
    if (compute_ci) {
      cat(sprintf("       Coverage_Overall = %.3f ± %.3f\n",
                  ranking$Coverage_Overall_Mean[i],
                  ranking$Coverage_Overall_SD[i]))
      cat(sprintf("       CI_Length_Overall = %.4f ± %.4f\n",
                  ranking$CI_Length_Overall_Mean[i],
                  ranking$CI_Length_Overall_SD[i]))
    }
    cat("\n")
  }
  
  sink()
  cat(sprintf("✓ 文本报告: %s\n", report_file))
  
  # ========== 5. 返回信息 ==========
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("✓ 所有结果已保存到目录: ", output_dir, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  cat("📁 文件清单:\n")
  cat("  1. detailed_simulation_results.txt\n")
  cat("     - 逐次仿真逐组详细数据\n")
  if (compute_ci) {
    cat("     - 包含: sim_id, method, group, mse, total_mse, nonzero_mse, zero_mse,\n")
    cat("             coverage, ci_length, coverage_overall, coverage_nonzero, coverage_zero,\n")
    cat("             ci_length_overall, ci_length_nonzero, ci_length_zero\n\n")
  } else {
    cat("     - 包含: sim_id, method, group, mse, total_mse, nonzero_mse, zero_mse\n\n")
  }
  
  cat("  2. summary_by_method.txt\n")
  cat("     - 按方法汇总 (4行 × 4种方法)\n")
  if (compute_ci) {
    cat("     - 包含: MSE统计量 + Coverage统计量 + CI Length统计量\n\n")
  } else {
    cat("     - 包含: MSE统计量\n\n")
  }
  
  cat("  3. summary_by_group.txt\n")
  cat("     - 按组别汇总 (20组 × 4种方法 = 80行)\n")
  if (compute_ci) {
    cat("     - 包含: MSE统计量 + Coverage统计量 + CI Length统计量\n\n")
  } else {
    cat("     - 包含: MSE统计量\n\n")
  }
  
  cat("  4. simulation_report.txt\n")
  cat("     - 完整的文本报告\n")
  if (compute_ci) {
    cat("     - 包含: 方法排名、各组别对比（MSE + Coverage + CI Length）\n\n")
  } else {
    cat("     - 包含: 方法排名、各组别对比（MSE）\n\n")
  }
  
  # 返回数据框供进一步使用
  return(list(
    detailed_df = detailed_df,
    summary_df = summary_df,
    group_summary_df = group_summary_df,
    output_dir = output_dir,
    files = c(txt_file, summary_txt, group_txt, report_file)
  ))
}


# ===================================================================
# -----------执行多次仿真-----------
# ===================================================================

cat("\n=== 多次重复仿真 ===\n\n")
cat("快速测试: repeated_results <- run_repeated_simulations(n_sim = 10, parallel = FALSE)\n")
cat("正式仿真: repeated_results <- run_repeated_simulations(n_sim = 100, parallel = TRUE)\n\n")

# 选择执行模式：

# 【选项1】快速测试（10次，顺序执行，约5-10分钟）
# repeated_results <- run_repeated_simulations(n_sim = 2, config = config, parallel = FALSE)

# 【选项2】正式仿真（100次，并行执行，推荐）
repeated_results <- run_repeated_simulations(n_sim = 100, config = config, parallel = TRUE)

# 查看结果
print_summary(repeated_results$summary)

# 绘制箱线图
plot_mse_boxplot(repeated_results$summary)

# 保存结果（可选）
save_simulation_results(repeated_results)
