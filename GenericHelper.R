library(BART)
library(C50) 
library(car)
#library(caret)
library(cluster)
library(corrplot)
library("data.table")
library(dendextend)
library(dunn.test) 
#library(extraTrees)
library(fingerprint) 
if (!require("devtools")) install.packages("devtools")
library(dplyr) 
library(e1071)
library(FactoMineR)
library(factoextra)
library(FNN)
if (!require("GGally")) install.packages("GGally")
library(GGally)
library(ggpubr)
library(ggplot2) 
library(glmnet)
library(glmmTMB)
library(PMCMRplus)
library(moments)
library(openxlsx)
library(partykit)
library(pls)
library(purrr)
library(PMCMRplus)
library(R6)
library(randomForest)
library(ranger)
library(rcdk)  
library(readODS) 
#library(recipes)
library(Rmpfr)
library(rpart.plot) 
library(stats)
library("sqldf")
library(tidyverse)
library(tidyr)
library(vegan)
library(xgboost)

GenericHelper <- R6Class(
  "GenericHelper",
  public = list(
    correlations_and_regressions = function(data,
                                            outFolder,
                                            filePrefix,
                                            cor_threshold = 0.7,
                                            p_threshold = 0.05,
                                            min_observations = 3,
                                            plot_width = 1600,
                                            plot_height = 1600,
                                            plot_res = 600) {
      # Check and create output folder
      if (!dir.exists(outFolder)) {
        dir.create(outFolder, recursive = TRUE)
        cat("Created folder:", outFolder, "\n")
      }
      
      if (!is.data.frame(data)) {
        stop("Parameter 'data' must be a data.frame")
      }
      
      # Select only numeric columns
      numeric_data <- data[, sapply(data, is.numeric), drop = FALSE]
      
      if (ncol(numeric_data) < 2) {
        stop("Data.frame must have at least 2 numeric columns for correlation")
      }
      
      cat("Correlation analysis for",
          ncol(numeric_data),
          "numeric variables\n")
      cat("Thresholds: correlation >=",
          cor_threshold,
          ", p-value <",
          p_threshold,
          "\n")
      
      # Function for calculating correlation and p-values
      calculate_pairwise = function(x, y, name_x, name_y) {
        complete_cases <- complete.cases(x, y)
        x_clean <- x[complete_cases]
        y_clean <- y[complete_cases]
        
        if (length(x_clean) < min_observations) {
          return(list(
            correlation = NA,
            pvalue = NA,
            n = length(x_clean)
          ))
        }
        
        cor_test <- cor.test(x_clean, y_clean)
        return(list(
          correlation = cor_test$estimate,
          pvalue = cor_test$p.value,
          n = length(x_clean)
        ))
      }
      
      # Initialize matrices
      n_vars <- ncol(numeric_data)
      var_names <- colnames(numeric_data)
      cor_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars)
      pvalue_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars)
      n_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars)
      
      colnames(cor_matrix) <- rownames(cor_matrix) <- var_names
      colnames(pvalue_matrix) <- rownames(pvalue_matrix) <- var_names
      colnames(n_matrix) <- rownames(n_matrix) <- var_names
      
      # Fill matrices
      for (i in 1:n_vars) {
        for (j in 1:n_vars) {
          if (i != j) {
            result <- calculate_pairwise(numeric_data[[i]], numeric_data[[j]], var_names[i], var_names[j])
            cor_matrix[i, j] <- result$correlation
            pvalue_matrix[i, j] <- result$pvalue
            n_matrix[i, j] <- result$n
          } else {
            cor_matrix[i, j] <- 1.0
            pvalue_matrix[i, j] <- 0.0
            n_matrix[i, j] <- nrow(numeric_data)
          }
        }
      }
      
      # 1. CORRELATION MATRIX VISUALIZATION
      plot_file <- file.path(outFolder, paste0(filePrefix, "_correlation_plot.png"))
      
      png(plot_file,
          width = plot_width,
          height = plot_height,
          res = plot_res)
      
      # Use complete correlation matrix
      cor_matrix_complete <- cor(numeric_data, use = "complete.obs")
      
      corrplot(
        cor_matrix_complete,
        method = "color",
        type = "upper",
        order = "hclust",
        tl.cex = 0.8,
        tl.col = "black",
        addCoef.col = "black",
        number.cex = 0.6,
        mar = c(0, 0, 2, 0),
        # Margins: bottom, left, top, right
        main = paste("Correlation Matrix -", filePrefix)
      )
      
      dev.off()
      cat("Correlation plot saved to:", plot_file, "\n")
      
      # 2. ADDITIONAL VISUALIZATION - Heatmap with ggplot2
      plot_file_heatmap <- file.path(outFolder, paste0(filePrefix, "_correlation_heatmap.png"))
      
      # Prepare data for ggplot2 heatmap
      cor_melt <- as.data.frame(as.table(cor_matrix_complete))
      names(cor_melt) <- c("Var1", "Var2", "Correlation")
      
      ggplot_heatmap <- ggplot(cor_melt, aes(Var1, Var2, fill = Correlation)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(
          low = "blue",
          high = "red",
          mid = "white",
          midpoint = 0,
          limit = c(-1, 1),
          space = "Lab",
          name = "Correlation"
        ) +
        geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(
            angle = 45,
            vjust = 1,
            hjust = 1
          ),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5)
        ) +
        labs(
          title = paste("Correlation Matrix -", filePrefix),
          subtitle = paste(n_vars, "numeric variables")
        )
      
      ggsave(
        plot_file_heatmap,
        ggplot_heatmap,
        width = 12,
        height = 10,
        dpi = 300
      )
      cat("Heatmap plot saved to:", plot_file_heatmap, "\n")
      
      # 3. Regression analysis for significant correlations
      regression_results <- list()
      significant_pairs <- 0
      
      for (i in 1:(n_vars - 1)) {
        for (j in (i + 1):n_vars) {
          if (!is.na(cor_matrix[i, j]) &&
              !is.na(pvalue_matrix[i, j]) &&
              abs(cor_matrix[i, j]) >= cor_threshold &&
              pvalue_matrix[i, j] < p_threshold) {
            significant_pairs <- significant_pairs + 1
            
            x <- numeric_data[[i]]
            y <- numeric_data[[j]]
            complete_cases <- complete.cases(x, y)
            x_clean <- x[complete_cases]
            y_clean <- y[complete_cases]
            
            # Linear regression
            lm_model <- lm(y_clean ~ x_clean)
            lm_summary <- summary(lm_model)
            
            # Save results
            regression_results[[paste(var_names[i], var_names[j], sep = "_vs_")]] <- data.table(
              predictor = var_names[i],
              response = var_names[j],
              correlation = round(cor_matrix[i, j], 4),
              p_value_cor = sprintf("%.4f", pvalue_matrix[i, j]),
              r_squared = round(lm_summary$r.squared, 4),
              adj_r_squared = round(lm_summary$adj.r.squared, 4),
              intercept = round(coef(lm_model)[1], 4),
              slope = round(coef(lm_model)[2], 4),
              intercept_p = sprintf("%.4f", coef(lm_summary)[1, 4]),
              slope_p = sprintf("%.4f", coef(lm_summary)[2, 4]),
              f_statistic = round(lm_summary$fstatistic[1], 2),
              df = paste(lm_summary$fstatistic[2], lm_summary$fstatistic[3], sep = ", "),
              n_observations = length(x_clean),
              equation = paste0("y = ", round(coef(lm_model)[2], 4), " * x + ", round(coef(lm_model)[1], 4))
            )
            
            # 4. REGRESSION VISUALIZATION for strongest correlations
            if (abs(cor_matrix[i, j]) >= 0.8) {
              # Only for very strong correlations
              reg_plot_file <- file.path(
                outFolder,
                paste0(
                  filePrefix,
                  "_regression_",
                  var_names[i],
                  "_vs_",
                  var_names[j],
                  ".png"
                )
              )
              
              regression_plot <- ggplot(data.frame(x = x_clean, y = y_clean), aes(x = x, y = y)) +
                geom_point(alpha = 0.6, color = "steelblue") +
                geom_smooth(method = "lm",
                            se = TRUE,
                            color = "red") +
                labs(
                  title = paste("Regression:", var_names[j], "~", var_names[i]),
                  x = var_names[i],
                  y = var_names[j],
                  subtitle = paste0(
                    "r = ",
                    round(cor_matrix[i, j], 3),
                    ", R² = ",
                    round(lm_summary$r.squared, 3),
                    ", p = ",
                    format.pval(pvalue_matrix[i, j], digits = 3)
                  )
                ) +
                theme_minimal() +
                theme(plot.title = element_text(face = "bold"))
              
              ggsave(
                reg_plot_file,
                regression_plot,
                width = 8,
                height = 6,
                dpi = 300
              )
            }
          }
        }
      }
      
      # 5. Save tabular results
      # Correlation matrix
      cor_dt <- as.data.table(cor_matrix, keep.rownames = "Variable")
      cor_file <- file.path(outFolder, paste0(filePrefix, ".correl.csv"))
      fwrite(cor_dt, cor_file)
      
      # P-values
      pvalue_dt <- as.data.table(pvalue_matrix, keep.rownames = "Variable")
      pvalue_file <- file.path(outFolder, paste0(filePrefix, ".pvalue.csv"))
      fwrite(pvalue_dt, pvalue_file)
      
      # Regression results
      if (length(regression_results) > 0) {
        regression_dt <- rbindlist(regression_results, fill = TRUE)
        regression_file <- file.path(outFolder, paste0(filePrefix, ".regression.csv"))
        fwrite(regression_dt, regression_file)
      } else {
        regression_dt <- data.table()
        regression_file <- NA
      }
      
      # 6. ADDITIONAL VISUALIZATION - Significant correlations
      if (length(regression_results) > 0) {
        # Network plot for significant correlations
        sig_cor_data <- regression_dt %>%
          select(predictor, response, correlation) %>%
          mutate(abs_cor = abs(correlation))
        
        if (nrow(sig_cor_data) > 0) {
          network_plot <- ggplot(sig_cor_data, aes(x = predictor, y = response)) +
            geom_point(aes(size = abs_cor, color = correlation), alpha = 0.7) +
            scale_size_continuous(range = c(3, 10), name = "Absolute correlation") +
            scale_color_gradient2(
              low = "blue",
              high = "red",
              mid = "white",
              midpoint = 0,
              name = "Correlation"
            ) +
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5)
            ) +
            labs(
              title = paste("Significant Correlations -", filePrefix),
              subtitle = paste("Threshold: |r| >=", cor_threshold, ", p <", p_threshold),
              x = "Predictor",
              y = "Response"
            )
          
          network_file <- file.path(outFolder,
                                    paste0(filePrefix, "_significant_correlations.png"))
          ggsave(
            network_file,
            network_plot,
            width = 10,
            height = 8,
            dpi = 300
          )
          cat("Network plot of significant correlations saved to:",
              network_file,
              "\n")
        }
      }
      
      # Analysis summary
      cat("\n=== ANALYSIS SUMMARY ===\n")
      cat("Total variables:", n_vars, "\n")
      cat("Total possible correlations:", choose(n_vars, 2), "\n")
      cat(
        "Significant correlations (r >=",
        cor_threshold,
        ", p <",
        p_threshold,
        "):",
        significant_pairs,
        "\n"
      )
      cat("Files saved to folder:", outFolder, "\n")
      cat(" - Correlations:", basename(cor_file), "\n")
      cat(" - P-values:", basename(pvalue_file), "\n")
      if (!is.na(regression_file)) {
        cat(" - Regressions:", basename(regression_file), "\n")
      }
      cat(" - Visualizations:\n")
      cat("   *", basename(plot_file), "\n")
      cat("   *", basename(plot_file_heatmap), "\n")
      if (exists("network_file")) {
        cat("   *", basename(network_file), "\n")
      }
      
      # Return results
      return(
        list(
          correlation_matrix = cor_dt,
          pvalue_matrix = pvalue_dt,
          regression_results = regression_dt,
          files = c(
            correlation = cor_file,
            pvalue = pvalue_file,
            regression = regression_file,
            plot_correlation = plot_file,
            plot_heatmap = plot_file_heatmap,
            plot_network = if (exists("network_file"))
              network_file
            else
              NA
          ),
          summary = list(
            n_variables = n_vars,
            n_significant_pairs = significant_pairs,
            correlation_threshold = cor_threshold,
            pvalue_threshold = p_threshold
          )
        )
      )
    }
    ,
    minmax_normalize = function(x) {
      if (length(unique(x)) == 1)
        return(rep(0, length(x))) # Handle constant variables
      (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    },
    
    # Function to calculate variance of min-max normalized data
    normalized_variance = function(x) {
      normalized <- self$minmax_normalize(x)
      var(normalized, na.rm = TRUE)
    }
    ,
 
    correlation_variable_selection = function(df,
                                              variables,
                                              correlation_method = "pearson",
                                              threshold = 0.8,
                                              alpha = 0.05,
                                              use_method = "pairwise.complete.obs") {
      #' Perform correlation-based variable selection
      #'
      #' @param df Data frame containing the variables
      #' @param variables Character vector of variable names to analyze
      #' @param correlation_method Correlation method: "pearson", "spearman", or "kendall" (default = "pearson")
      #' @param threshold Correlation threshold (absolute value, default = 0.8)
      #' @param alpha Significance level for correlation test (default = 0.05)
      #' @param use_method Method for handling missing values (default = "pairwise.complete.obs")
      #'
      #' @return List containing:
      #'   - removed_vars: Character vector of removed variables
      #'   - remaining_vars: Character vector of remaining variables
      #'   - all_vars: Character vector of all input variables
      #'   - similar_vars: List showing which variables were similar to removed ones
      #'   - correlation_method: Method used for correlation calculation
      #'   - threshold: Threshold used for variable removal
      
      # Input validation
      if (!all(variables %in% names(df))) {
        stop("Some specified variables are not in the data frame")
      }
      
      if (length(variables) < 2) {
        stop("At least 2 variables are required for correlation analysis")
      }
      
      # Validate correlation method
      valid_methods <- c("pearson", "spearman", "kendall")
      if (!correlation_method %in% valid_methods) {
        stop("Invalid correlation_method. Choose from: ", 
             paste(valid_methods, collapse = ", "))
      }
      
      # Subset the data to only include specified variables
      data_subset <- df[, variables, drop = FALSE]
      
      # Remove non-numeric columns
      numeric_vars <- sapply(data_subset, is.numeric)
      if (sum(numeric_vars) < 2) {
        stop("At least 2 numeric variables are required")
      }
      
      # Warn about non-numeric variables
      if (!all(numeric_vars)) {
        non_numeric_vars <- variables[!numeric_vars]
        warning("The following non-numeric variables will be excluded: ",
                paste(non_numeric_vars, collapse = ", "))
      }
      
      data_subset <- data_subset[, numeric_vars, drop = FALSE]
      all_vars <- names(data_subset)
      remaining_vars <- all_vars
      removed_vars <- character(0)
      similar_vars <- list()
      
      # Function to min-max normalize a vector
      minmax_normalize <- function(x) {
        if (length(unique(x)) == 1 || all(is.na(x)))
          return(rep(0, length(x))) # Handle constant variables or all NA
        (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
      }
      
      # Function to calculate variance of min-max normalized data
      normalized_variance <- function(x) {
        normalized <- minmax_normalize(x)
        var(normalized, na.rm = TRUE)
      }
      
      # Create correlation matrix with specified method
      cor_matrix <- cor(data_subset, 
                        method = correlation_method, 
                        use = use_method)
      
      # Get correlation symbol for messages
      cor_symbol <- switch(correlation_method,
                           "pearson" = "r",
                           "spearman" = "ρ", 
                           "kendall" = "τ")
      
      message("Performing variable selection using ", correlation_method, " correlation (", cor_symbol, ")")
      message("Threshold: |", cor_symbol, "| >= ", threshold)
      message("Significance level: α = ", alpha)
      
      # Find pairs with high correlation
      high_cor_pairs <- which(abs(cor_matrix) >= threshold &
                                upper.tri(cor_matrix),
                              arr.ind = TRUE)
      
      if (nrow(high_cor_pairs) == 0) {
        message("No variable pairs found with correlation >= ", threshold)
        return(
          list(
            removed_vars = removed_vars,
            remaining_vars = remaining_vars,
            all_vars = all_vars,
            similar_vars = similar_vars,
            correlation_method = correlation_method,
            threshold = threshold,
            alpha = alpha
          )
        )
      }
      
      message("Found ", nrow(high_cor_pairs), " variable pairs with |", cor_symbol, "| >= ", threshold)
      
      # Process each high correlation pair
      for (i in 1:nrow(high_cor_pairs)) {
        row_idx <- high_cor_pairs[i, 1]
        col_idx <- high_cor_pairs[i, 2]
        
        var1 <- all_vars[row_idx]
        var2 <- all_vars[col_idx]
        
        # Skip if either variable has already been removed
        if (!var1 %in% remaining_vars || !var2 %in% remaining_vars) {
          next
        }
        
        # Get actual correlation value
        cor_value <- cor_matrix[row_idx, col_idx]
        
        # Test correlation significance with specified method
        cor_test <- tryCatch({
          cor.test(data_subset[[var1]], data_subset[[var2]], 
                   method = correlation_method, 
                   use = use_method)
        }, error = function(e) {
          warning("Could not calculate p-value for pair (", var1, ", ", var2, "): ", e$message)
          return(list(p.value = NA))
        })
        
        # Check if correlation is significant (or if p-value is NA, proceed with caution)
        if (is.na(cor_test$p.value) || cor_test$p.value <= alpha) {
          # Calculate variances of min-max normalized variables
          var1_variance <- normalized_variance(data_subset[[var1]])
          var2_variance <- normalized_variance(data_subset[[var2]])
          
          # Determine which variable to remove (remove the one with lower normalized variance)
          if (var1_variance <= var2_variance) {
            variable_to_remove <- var1
            similar_variable <- var2
            variance_ratio <- ifelse(var2_variance > 0, var1_variance / var2_variance, 0)
          } else {
            variable_to_remove <- var2
            similar_variable <- var1
            variance_ratio <- ifelse(var1_variance > 0, var2_variance / var1_variance, 0)
          }
          
          # Format p-value for display
          p_value_display <- ifelse(is.na(cor_test$p.value), "NA", round(cor_test$p.value, 4))
          
          removed_reason <- paste0(
            "Similar to '", similar_variable, 
            "' (", cor_symbol, " = ", round(cor_value, 3),
            ", p = ", p_value_display,
            ", variance ratio = ", round(variance_ratio, 3), ")"
          )
          
          # Add to results
          removed_vars <- c(removed_vars, variable_to_remove)
          remaining_vars <- setdiff(remaining_vars, variable_to_remove)
          similar_vars[[variable_to_remove]] <- list(
            similar_to = similar_variable,
            correlation = cor_value,
            p_value = cor_test$p.value,
            variance_ratio = variance_ratio,
            reason = removed_reason,
            correlation_method = correlation_method
          )
          
          message("Removed '", variable_to_remove, "' - ", removed_reason)
        }
      }
      
      # Summary message
      if (length(removed_vars) > 0) {
        message("\nVariable selection completed:")
        message("  - Removed ", length(removed_vars), " variables: ", paste(removed_vars, collapse = ", "))
        message("  - Retained ", length(remaining_vars), " variables: ", paste(remaining_vars, collapse = ", "))
        message("  - Reduction: ", length(all_vars), " → ", length(remaining_vars), 
                " (", round((1 - length(remaining_vars)/length(all_vars)) * 100, 1), "% reduction)")
      } else {
        message("No variables were removed based on the specified criteria")
      }
      
      # Return comprehensive results
      return(
        list(
          removed_vars = removed_vars,
          remaining_vars = remaining_vars,
          all_vars = all_vars,
          similar_vars = similar_vars,
          correlation_method = correlation_method,
          threshold = threshold,
          alpha = alpha,
          summary = list(
            initial_count = length(all_vars),
            final_count = length(remaining_vars),
            reduction_percent = round((1 - length(remaining_vars)/length(all_vars)) * 100, 1),
            removed_count = length(removed_vars)
          )
        )
      )
    }
    
    
    ,
 
    create_group_boxplots_overview = function(df,
                                              output_folder,
                                              group_var,
                                              target_vars,
                                              dpi = 600,
                                              width = NULL,
                                              height = NULL) {
      # Load required packages
      if (!require(ggplot2, quietly = TRUE)) {
        stop("ggplot2 package is required. Please install it using install.packages('ggplot2')")
      }
      if (!require(gridExtra, quietly = TRUE)) {
        stop(
          "gridExtra package is required. Please install it using install.packages('gridExtra')"
        )
      }
      if (!require(cowplot, quietly = TRUE)) {
        stop("cowplot package is required. Please install it using install.packages('cowplot')")
      }
      
      # Validate inputs
      if (!is.data.frame(df)) {
        stop("The first argument must be a data frame")
      }
      
      # Create output folder if it doesn't exist
      if (!dir.exists(output_folder)) {
        cat("Creating output folder:", output_folder, "\n")
        dir.create(output_folder,
                   recursive = TRUE,
                   showWarnings = FALSE)
      }
      
      # Check if folder was created successfully
      if (!dir.exists(output_folder)) {
        stop("Failed to create output folder: ", output_folder)
      }
      
      if (!group_var %in% names(df)) {
        stop("Group variable '", group_var, "' not found in data frame")
      }
      
      missing_vars <- target_vars[!target_vars %in% names(df)]
      if (length(missing_vars) > 0) {
        stop(
          "The following target variables not found in data frame: ",
          paste(missing_vars, collapse = ", ")
        )
      }
      
      # Validate resolution parameters
      if (!is.null(width) && (!is.numeric(width) || width <= 0)) {
        stop("Width must be a positive numeric value")
      }
      if (!is.null(height) &&
          (!is.numeric(height) || height <= 0)) {
        stop("Height must be a positive numeric value")
      }
      if (!is.numeric(dpi) || dpi <= 0) {
        stop("DPI must be a positive numeric value")
      }
      
      # Convert group variable to factor for better plotting
      df[[group_var]] <- as.factor(df[[group_var]])
      group_levels <- levels(df[[group_var]])
      
      n_plots <- length(target_vars)
      cat("Creating", n_plots, "boxplots...\n")
      
      # Optimize grid layout for large numbers of plots
      if (n_plots <= 4) {
        n_cols <- 2
        n_rows <- ceiling(n_plots / 2)
      } else if (n_plots <= 9) {
        n_cols <- 3
        n_rows <- ceiling(n_plots / 3)
      } else if (n_plots <= 16) {
        n_cols <- 4
        n_rows <- ceiling(n_plots / 4)
      } else if (n_plots <= 25) {
        n_cols <- 5
        n_rows <- ceiling(n_plots / 5)
      } else {
        n_cols <- 6
        n_rows <- ceiling(n_plots / 6)
      }
      
      # Increase base dimensions by 10%
      base_width_per_col <- 12 * 1.1  # 13.2 inches per column
      base_height_per_row <- 10 * 1.1 # 11 inches per row
      extra_legend_height <- 2 * 1.1  # 2.2 inches for legend
      
      # Calculate dimensions optimized for visibility (increased by 10%)
      if (is.null(width)) {
        width <- base_width_per_col * n_cols
      }
      
      if (is.null(height)) {
        height <- (base_height_per_row * n_rows) + extra_legend_height
      }
      
      cat("Grid layout:", n_rows, "rows ×", n_cols, "columns\n")
      cat("Image dimensions:",
          round(width, 1),
          "×",
          round(height, 1),
          "inches\n")
      cat("Base dimensions increased by 10%\n")
      
      # Create individual boxplots for each target variable
      plot_list <- list()
      
      for (i in seq_along(target_vars)) {
        target_var <- target_vars[i]
        
        p <- ggplot(df, aes(x = .data[[group_var]], y = .data[[target_var]])) +
          geom_boxplot(
            aes(fill = .data[[group_var]]),
            alpha = 0.8,
            outlier.shape = 16,
            outlier.alpha = 0.5,
            size = 0.4,
            # Decreased box line thickness
            width = 0.9,
            # Increased box width
            fatten = 1.5          # Slightly decreased median line thickness
          ) +
          geom_jitter(
            width = 0.2,
            alpha = 0.3,
            size = 0.5,
            height = 0
          ) +
          labs(
            title = target_var,
            x = NULL,
            # Remove x-axis label
            y = NULL,
            fill = NULL  # Remove legend title from individual plots
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(
              size = 15,
              face = "bold",
              hjust = 0.5,
              margin = margin(b = 12)
            ),
            # Increased 10%
            axis.text.x = element_blank(),
            # Remove x-axis text
            axis.ticks.x = element_blank(),
            # Remove x-axis ticks
            axis.text.y = element_text(size = 12),
            # Increased 10%
            axis.title.x = element_blank(),
            # Remove x-axis title
            legend.position = "none",
            # Remove legend from individual plots
            panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(
              color = "grey70",
              fill = NA,
              linewidth = 0.5
            ),
            plot.margin = margin(16, 16, 16, 16, "pt"),
            # Increased 10%
            panel.background = element_rect(fill = "white", color = NA)
          )
        
        # Use viridis colors if available, otherwise use default
        if (require(viridis, quietly = TRUE)) {
          p <- p + scale_fill_viridis_d(begin = 0.2, end = 0.8)
        } else {
          p <- p + scale_fill_brewer(palette = "Set2")
        }
        
        plot_list[[i]] <- p
      }
      
      # Create a shared legend from the first plot
      legend_plot <- ggplot(df, aes(x = .data[[group_var]], y = .data[[target_vars[1]]])) +
        geom_boxplot(aes(fill = .data[[group_var]]), alpha = 0.8) +
        labs(fill = group_var) +
        theme(
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(size = 13, face = "bold"),
          # Increased 10%
          legend.text = element_text(size = 12),
          # Increased 10%
          legend.key.size = unit(0.9, "cm"),
          # Increased 10%
          legend.margin = margin(t = 11, b = 11) # Increased 10%
        )
      
      if (require(viridis, quietly = TRUE)) {
        legend_plot <- legend_plot + scale_fill_viridis_d(begin = 0.2, end = 0.8)
      } else {
        legend_plot <- legend_plot + scale_fill_brewer(palette = "Set2")
      }
      
      # Extract the legend
      legend_grob <- cowplot::get_legend(legend_plot)
      
      # Create text grob for group levels information
      levels_grob <- grid::textGrob(
        paste("Levels:", toString(group_levels)),
        gp = grid::gpar(fontsize = 12, fontface = "italic"),
        # Increased 10%
        just = "center"
      )
      
      # Arrange everything: plots, then legend, then levels info
      combined_plot <- gridExtra::grid.arrange(
        # Main plots
        gridExtra::arrangeGrob(
          grobs = plot_list,
          ncol = n_cols,
          nrow = n_rows,
          top = grid::textGrob(
            paste("Boxplot Overview:", n_plots, "Variables by", group_var),
            gp = grid::gpar(fontsize = 18, fontface = "bold") # Increased 10%
          )
        ),
        # Legend
        legend_grob,
        # Levels information
        levels_grob,
        nrow = 3,
        heights = c(0.88, 0.06, 0.06)  # Adjusted heights for larger plots
      )
      
      # Create output filename
      output_file <- file.path(output_folder, paste0(group_var, "_boxplots_overview.png"))
      
      # Save the plot with limitsize = FALSE to allow large dimensions
      ggsave(
        filename = output_file,
        plot = combined_plot,
        width = width,
        height = height,
        dpi = dpi,
        bg = "white",
        limitsize = FALSE
      )
      
      # Print confirmation message
      cat("\n✅ Boxplot overview successfully created!\n")
      cat("📁 Output file:", output_file, "\n")
      cat("📐 Final resolution:",
          round(width, 1),
          "×",
          round(height, 1),
          "inches\n")
      cat("🖼️  DPI:", dpi, "\n")
      cat("📊 Number of plots:", n_plots, "\n")
      cat("🔢 Grid layout:", n_rows, "rows ×", n_cols, "columns\n")
      cat("🎨 Group variable:", group_var, "\n")
      cat("🔤 Group levels:", toString(group_levels), "\n")
      cat("💾 File size:", round(file.size(output_file) / 1024, 1), "KB\n")
      cat("📈 All dimensions increased by 10%\n")
      
      return(invisible(combined_plot))
    }
    
    ,
    
    create_variable_boxplots_overview = function(df,
                                                 output_file,
                                                 target_vars,
                                                 dpi = 600,
                                                 width = NULL,
                                                 height = NULL) {
      # Load required packages
      if (!require(ggplot2, quietly = TRUE)) {
        stop("ggplot2 package is required. Please install it using install.packages('ggplot2')")
      }
      if (!require(gridExtra, quietly = TRUE)) {
        stop(
          "gridExtra package is required. Please install it using install.packages('gridExtra')"
        )
      }
      
      # Validate inputs
      if (!is.data.frame(df)) {
        stop("The first argument must be a data frame")
      }
      
      # Validate output file
      if (!is.character(output_file) || length(output_file) != 1) {
        stop("output_file must be a single character string")
      }
      
      # Ensure output file has .png extension
      if (!grepl("\\.png$", output_file, ignore.case = TRUE)) {
        stop("output_file must have a .png extension")
      }
      
      # Extract directory from output file path
      output_folder <- dirname(output_file)
      
      # Create output folder if it doesn't exist
      if (!dir.exists(output_folder)) {
        cat("Creating output folder:", output_folder, "\n")
        dir.create(output_folder,
                   recursive = TRUE,
                   showWarnings = FALSE)
      }
      
      # Check if folder was created successfully
      if (!dir.exists(output_folder)) {
        stop("Failed to create output folder: ", output_folder)
      }
      
      missing_vars <- target_vars[!target_vars %in% names(df)]
      if (length(missing_vars) > 0) {
        stop(
          "The following target variables not found in data frame: ",
          paste(missing_vars, collapse = ", ")
        )
      }
      
      # Validate resolution parameters
      if (!is.null(width) && (!is.numeric(width) || width <= 0)) {
        stop("Width must be a positive numeric value")
      }
      if (!is.null(height) &&
          (!is.numeric(height) || height <= 0)) {
        stop("Height must be a positive numeric value")
      }
      if (!is.numeric(dpi) || dpi <= 0) {
        stop("DPI must be a positive numeric value")
      }
      
      n_plots <- length(target_vars)
      cat("Creating", n_plots, "boxplots...\n")
      
      # Optimize grid layout for large numbers of plots
      if (n_plots <= 4) {
        n_cols <- 2
        n_rows <- ceiling(n_plots / 2)
      } else if (n_plots <= 9) {
        n_cols <- 3
        n_rows <- ceiling(n_plots / 3)
      } else if (n_plots <= 16) {
        n_cols <- 4
        n_rows <- ceiling(n_plots / 4)
      } else if (n_plots <= 25) {
        n_cols <- 5
        n_rows <- ceiling(n_plots / 5)
      } else {
        n_cols <- 6
        n_rows <- ceiling(n_plots / 6)
      }
      
      # Increase base dimensions by 10%
      base_width_per_col <- 12 * 1.1  # 13.2 inches per column
      base_height_per_row <- 10 * 1.1 # 11 inches per row
      
      # Calculate dimensions optimized for visibility (increased by 10%)
      if (is.null(width)) {
        width <- base_width_per_col * n_cols
      }
      
      if (is.null(height)) {
        height <- base_height_per_row * n_rows
      }
      
      cat("Grid layout:", n_rows, "rows ×", n_cols, "columns\n")
      cat("Image dimensions:",
          round(width, 1),
          "×",
          round(height, 1),
          "inches\n")
      cat("Base dimensions increased by 10%\n")
      
      # Create individual boxplots for each target variable
      plot_list <- list()
      
      for (i in seq_along(target_vars)) {
        target_var <- target_vars[i]
        
        # Create a dummy constant for x-axis
        p <- ggplot(df, aes(x = factor(1), y = .data[[target_var]])) +
          geom_boxplot(
            fill = "steelblue",
            alpha = 0.8,
            outlier.shape = 16,
            outlier.alpha = 0.5,
            size = 0.4,
            # Decreased box line thickness
            width = 0.7,
            # Box width
            fatten = 1.5          # Median line thickness
          ) +
          geom_jitter(
            width = 0.2,
            alpha = 0.3,
            size = 0.5,
            height = 0,
            color = "steelblue"
          ) +
          labs(title = target_var, x = NULL, # Remove x-axis label
               y = NULL) +
          theme_minimal() +
          theme(
            plot.title = element_text(
              size = 15,
              face = "bold",
              hjust = 0.5,
              margin = ggplot2::margin(b = 12)
            ),
            axis.text.x = element_blank(),
            # Remove x-axis text
            axis.ticks.x = element_blank(),
            # Remove x-axis ticks
            axis.text.y = element_text(size = 12),
            axis.title.x = element_blank(),
            # Remove x-axis title
            panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(
              color = "grey70",
              fill = NA,
              linewidth = 0.5
            ),
            plot.margin = ggplot2::margin(16, 16, 16, 16, "pt"),
            panel.background = element_rect(fill = "white", color = NA)
          )
        
        plot_list[[i]] <- p
      }
      
      # Arrange all plots in a grid
      combined_plot <- gridExtra::grid.arrange(
        grobs = plot_list,
        ncol = n_cols,
        nrow = n_rows,
        top = grid::textGrob(
          paste("Boxplot Overview:", n_plots, "Variables"),
          gp = grid::gpar(fontsize = 18, fontface = "bold")
        )
      )
      
      # Save the plot with limitsize = FALSE to allow large dimensions
      ggplot2::ggsave(
        filename = output_file,
        plot = combined_plot,
        width = width,
        height = height,
        dpi = dpi,
        bg = "white",
        limitsize = FALSE
      )
      
      # Print confirmation message
      cat("\n✅ Boxplot overview successfully created!\n")
      cat("📁 Output file:", output_file, "\n")
      cat("📐 Final resolution:",
          round(width, 1),
          "×",
          round(height, 1),
          "inches\n")
      cat("🖼️  DPI:", dpi, "\n")
      cat("📊 Number of plots:", n_plots, "\n")
      cat("🔢 Grid layout:", n_rows, "rows ×", n_cols, "columns\n")
      cat("💾 File size:", round(file.size(output_file) / 1024, 1), "KB\n")
      
      return(invisible(combined_plot))
    }
    
    ,
    
    create_individual_boxplots = function(df,
                                          output_folder,
                                          group_var = NULL,
                                          target_vars,
                                          dpi = 600,
                                          width = 10,
                                          height = 8) {
      # Load required packages
      if (!require(ggplot2, quietly = TRUE)) {
        stop("ggplot2 package is required. Please install it using install.packages('ggplot2')")
      }
      
      # Validate inputs
      if (!is.data.frame(df)) {
        stop("The first argument must be a data frame")
      }
      
      # Create output folder if it doesn't exist
      if (!dir.exists(output_folder)) {
        cat("Creating output folder:", output_folder, "\n")
        dir.create(output_folder,
                   recursive = TRUE,
                   showWarnings = FALSE)
      }
      
      # Check if folder was created successfully
      if (!dir.exists(output_folder)) {
        stop("Failed to create output folder: ", output_folder)
      }
      
      # Validate target variables
      missing_vars <- target_vars[!target_vars %in% names(df)]
      if (length(missing_vars) > 0) {
        stop(
          "The following target variables not found in data frame: ",
          paste(missing_vars, collapse = ", ")
        )
      }
      
      # Validate group variable if provided
      if (!is.null(group_var)) {
        if (!group_var %in% names(df)) {
          stop("Group variable '", group_var, "' not found in data frame")
        }
        # Convert group variable to factor for better plotting
        df[[group_var]] <- as.factor(df[[group_var]])
        group_levels <- levels(df[[group_var]])
      }
      
      # Validate resolution parameters
      if (!is.numeric(width) || width <= 0) {
        stop("Width must be a positive numeric value")
      }
      if (!is.numeric(height) || height <= 0) {
        stop("Height must be a positive numeric value")
      }
      if (!is.numeric(dpi) || dpi <= 0) {
        stop("DPI must be a positive numeric value")
      }
      
      n_plots <- length(target_vars)
      cat("Creating", n_plots, "individual boxplot files...\n")
      
      # Create individual boxplots for each target variable
      created_files <- character(n_plots)
      
      for (i in seq_along(target_vars)) {
        target_var <- target_vars[i]
        
        if (!is.null(group_var)) {
          # Create grouped boxplot
          p <- ggplot(df, aes(x = .data[[group_var]], y = .data[[target_var]])) +
            geom_boxplot(
              aes(fill = .data[[group_var]]),
              alpha = 0.8,
              outlier.shape = 16,
              outlier.alpha = 0.5,
              size = 0.4,
              width = 0.8,
              fatten = 1.5
            ) +
            geom_jitter(
              width = 0.2,
              alpha = 0.6,
              # Increased alpha for better visibility
              size = 1.5,
              # Increased dot size (was 0.5)
              height = 0
            ) +
            labs(
              title = paste("Boxplot of", target_var, "by", group_var),
              x = group_var,
              y = target_var,
              fill = group_var
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(
                size = 16,
                face = "bold",
                hjust = 0.5,
                margin = margin(b = 12)
              ),
              axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                size = 12
              ),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 13, face = "bold"),
              legend.position = "right",
              panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(
                color = "grey70",
                fill = NA,
                linewidth = 0.5
              ),
              plot.margin = margin(16, 16, 16, 16, "pt"),
              panel.background = element_rect(fill = "white", color = NA)
            )
          
          # Use viridis colors if available, otherwise use default
          if (require(viridis, quietly = TRUE)) {
            p <- p + scale_fill_viridis_d(begin = 0.2, end = 0.8)
          } else {
            p <- p + scale_fill_brewer(palette = "Set2")
          }
          
          # Create filename for grouped plot
          filename <- paste0(target_var, "_by_", group_var, ".png")
          
        } else {
          # Create ungrouped boxplot
          p <- ggplot(df, aes(x = factor(1), y = .data[[target_var]])) +
            geom_boxplot(
              fill = "steelblue",
              alpha = 0.8,
              outlier.shape = 16,
              outlier.alpha = 0.5,
              size = 0.4,
              width = 0.6,
              fatten = 1.5
            ) +
            geom_jitter(
              width = 0.2,
              alpha = 0.6,
              # Increased alpha for better visibility
              size = 1.5,
              # Increased dot size (was 0.5)
              height = 0,
              color = "steelblue"
            ) +
            labs(
              title = paste("Boxplot of", target_var),
              x = NULL,
              y = target_var
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(
                size = 16,
                face = "bold",
                hjust = 0.5,
                margin = margin(b = 12)
              ),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 13, face = "bold"),
              panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(
                color = "grey70",
                fill = NA,
                linewidth = 0.5
              ),
              plot.margin = margin(16, 16, 16, 16, "pt"),
              panel.background = element_rect(fill = "white", color = NA)
            )
          
          # Create filename for ungrouped plot
          filename <- paste0(target_var, "_boxplot.png")
        }
        
        # Create full output path
        output_file <- file.path(output_folder, filename)
        
        # Save the individual plot
        ggsave(
          filename = output_file,
          plot = p,
          width = width,
          height = height,
          dpi = dpi,
          bg = "white"
        )
        
        created_files[i] <- output_file
        cat("Created:", filename, "\n")
      }
      
      # Print summary message
      cat("\n✅ All individual boxplots successfully created!\n")
      cat("📁 Output folder:", output_folder, "\n")
      cat("📐 Plot dimensions:", width, "×", height, "inches\n")
      cat("🖼️  DPI:", dpi, "\n")
      cat("📊 Number of plots created:", n_plots, "\n")
      cat("🔵 Dot size: Increased to 1.5 for better visibility\n")
      if (!is.null(group_var)) {
        cat("🎨 Group variable:", group_var, "\n")
      } else {
        cat("🎨 Plot type: Ungrouped boxplots\n")
      }
      
      return(invisible(created_files))
    }
    
    ,
    create_individual_histograms = function(df,
                                            output_folder,
                                            group_var = NULL,
                                            target_vars,
                                            dpi = 600,
                                            width = 10,
                                            height = 8,
                                            bins = 30) {
      # Load required packages
      if (!require(ggplot2, quietly = TRUE)) {
        stop("ggplot2 package is required. Please install it using install.packages('ggplot2')")
      }
      
      # Validate inputs
      if (!is.data.frame(df)) {
        stop("The first argument must be a data frame")
      }
      
      # Create output folder if it doesn't exist
      if (!dir.exists(output_folder)) {
        cat("Creating output folder:", output_folder, "\n")
        dir.create(output_folder,
                   recursive = TRUE,
                   showWarnings = FALSE)
      }
      
      # Check if folder was created successfully
      if (!dir.exists(output_folder)) {
        stop("Failed to create output folder: ", output_folder)
      }
      
      # Validate target variables
      missing_vars <- target_vars[!target_vars %in% names(df)]
      if (length(missing_vars) > 0) {
        stop(
          "The following target variables not found in data frame: ",
          paste(missing_vars, collapse = ", ")
        )
      }
      
      # Validate group variable if provided
      if (!is.null(group_var)) {
        if (!group_var %in% names(df)) {
          stop("Group variable '", group_var, "' not found in data frame")
        }
        # Convert group variable to factor for better plotting
        df[[group_var]] <- as.factor(df[[group_var]])
        group_levels <- levels(df[[group_var]])
      }
      
      # Validate resolution parameters
      if (!is.numeric(width) || width <= 0) {
        stop("Width must be a positive numeric value")
      }
      if (!is.numeric(height) || height <= 0) {
        stop("Height must be a positive numeric value")
      }
      if (!is.numeric(dpi) || dpi <= 0) {
        stop("DPI must be a positive numeric value")
      }
      if (!is.numeric(bins) || bins <= 0) {
        stop("Bins must be a positive numeric value")
      }
      
      n_plots <- length(target_vars)
      cat("Creating", n_plots, "individual histogram files...\n")
      
      # Create individual histograms for each target variable
      created_files <- character(n_plots)
      
      for (i in seq_along(target_vars)) {
        target_var <- target_vars[i]
        
        if (!is.null(group_var)) {
          # Create grouped histogram
          p <- ggplot(df, aes(x = .data[[target_var]], fill = .data[[group_var]])) +
            geom_histogram(
              alpha = 0.7,
              position = "identity",
              bins = bins,
              color = "white",
              size = 0.2
            ) +
            labs(
              title = paste("Histogram of", target_var, "by", group_var),
              x = target_var,
              y = "Frequency",
              fill = group_var
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(
                size = 16,
                face = "bold",
                hjust = 0.5,
                margin = margin(b = 12)
              ),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 13, face = "bold"),
              legend.position = "right",
              panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(
                color = "grey70",
                fill = NA,
                linewidth = 0.5
              ),
              plot.margin = margin(16, 16, 16, 16, "pt"),
              panel.background = element_rect(fill = "white", color = NA)
            )
          
          # Use viridis colors if available, otherwise use default
          if (require(viridis, quietly = TRUE)) {
            p <- p + scale_fill_viridis_d(begin = 0.2, end = 0.8)
          } else {
            p <- p + scale_fill_brewer(palette = "Set2")
          }
          
          # Create filename for grouped histogram
          filename <- paste0(target_var, "_histogram_by_", group_var, ".png")
          
        } else {
          # Create ungrouped histogram
          p <- ggplot(df, aes(x = .data[[target_var]])) +
            geom_histogram(
              fill = "steelblue",
              alpha = 0.8,
              bins = bins,
              color = "white",
              size = 0.3
            ) +
            labs(
              title = paste("Histogram of", target_var),
              x = target_var,
              y = "Frequency"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(
                size = 16,
                face = "bold",
                hjust = 0.5,
                margin = margin(b = 12)
              ),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 13, face = "bold"),
              panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(
                color = "grey70",
                fill = NA,
                linewidth = 0.5
              ),
              plot.margin = margin(16, 16, 16, 16, "pt"),
              panel.background = element_rect(fill = "white", color = NA)
            )
          
          # Create filename for ungrouped histogram
          filename <- paste0(target_var, "_histogram.png")
        }
        
        # Create full output path
        output_file <- file.path(output_folder, filename)
        
        # Save the individual histogram
        ggsave(
          filename = output_file,
          plot = p,
          width = width,
          height = height,
          dpi = dpi,
          bg = "white"
        )
        
        created_files[i] <- output_file
        cat("Created:", filename, "\n")
      }
      
      # Print summary message
      cat("\n✅ All individual histograms successfully created!\n")
      cat("📁 Output folder:", output_folder, "\n")
      cat("📐 Plot dimensions:", width, "×", height, "inches\n")
      cat("🖼️  DPI:", dpi, "\n")
      cat("📊 Number of histograms created:", n_plots, "\n")
      cat("📦 Number of bins:", bins, "\n")
      if (!is.null(group_var)) {
        cat("🎨 Group variable:", group_var, "\n")
      } else {
        cat("🎨 Plot type: Ungrouped histograms\n")
      }
      
      return(invisible(created_files))
    }
    
    ,
    create_histograms_overview = function(df,
                                          output_file,
                                          target_vars,
                                          dpi = 600,
                                          width = NULL,
                                          height = NULL,
                                          bins = 30) {
      # Load required packages
      if (!require(ggplot2, quietly = TRUE)) {
        stop("ggplot2 package is required. Please install it using install.packages('ggplot2')")
      }
      if (!require(gridExtra, quietly = TRUE)) {
        stop(
          "gridExtra package is required. Please install it using install.packages('gridExtra')"
        )
      }
      
      # Validate inputs
      if (!is.data.frame(df)) {
        stop("The first argument must be a data frame")
      }
      
      # Validate output file path
      if (!is.character(output_file) || length(output_file) != 1) {
        stop("Output file must be a single character string")
      }
      
      # Ensure output file has .png extension
      if (!grepl("\\.png$", output_file, ignore.case = TRUE)) {
        output_file <- paste0(tools::file_path_sans_ext(output_file), ".png")
        cat("Added .png extension to output file:", output_file, "\n")
      }
      
      # Create output directory if it doesn't exist
      output_dir <- dirname(output_file)
      if (!dir.exists(output_dir) && output_dir != ".") {
        cat("Creating output directory:", output_dir, "\n")
        dir.create(output_dir,
                   recursive = TRUE,
                   showWarnings = FALSE)
      }
      
      # Check if directory was created successfully (if needed)
      if (output_dir != "." && !dir.exists(output_dir)) {
        stop("Failed to create output directory: ", output_dir)
      }
      
      # Validate target variables
      missing_vars <- target_vars[!target_vars %in% names(df)]
      if (length(missing_vars) > 0) {
        stop(
          "The following target variables not found in data frame: ",
          paste(missing_vars, collapse = ", ")
        )
      }
      
      # Validate resolution parameters
      if (!is.null(width) && (!is.numeric(width) || width <= 0)) {
        stop("Width must be a positive numeric value")
      }
      if (!is.null(height) && (!is.numeric(height) || height <= 0)) {
        stop("Height must be a positive numeric value")
      }
      if (!is.numeric(dpi) || dpi <= 0) {
        stop("DPI must be a positive numeric value")
      }
      if (!is.numeric(bins) || bins <= 0) {
        stop("Bins must be a positive numeric value")
      }
      
      n_plots <- length(target_vars)
      cat("Creating histogram overview with", n_plots, "variables...\n")
      
      # Optimize grid layout for large numbers of plots
      if (n_plots <= 4) {
        n_cols <- 2
        n_rows <- ceiling(n_plots / 2)
      } else if (n_plots <= 9) {
        n_cols <- 3
        n_rows <- ceiling(n_plots / 3)
      } else if (n_plots <= 16) {
        n_cols <- 4
        n_rows <- ceiling(n_plots / 4)
      } else if (n_plots <= 25) {
        n_cols <- 5
        n_rows <- ceiling(n_plots / 5)
      } else {
        n_cols <- 6
        n_rows <- ceiling(n_plots / 6)
      }
      
      # Calculate dimensions optimized for visibility
      base_width_per_col <- 5 * 1.1
      base_height_per_row <- 4 * 1.1
      
      if (is.null(width)) {
        width <- base_width_per_col * n_cols
      }
      
      if (is.null(height)) {
        height <- base_height_per_row * n_rows
      }
      
      cat("Grid layout:", n_rows, "rows ×", n_cols, "columns\n")
      cat("Image dimensions:",
          round(width, 1),
          "×",
          round(height, 1),
          "inches\n")
      
      # Create individual histograms for each target variable
      plot_list <- list()
      
      for (i in seq_along(target_vars)) {
        target_var <- target_vars[i]
        
        # Create ungrouped histogram for overview
        p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[target_var]])) +
          ggplot2::geom_histogram(
            fill = "steelblue",
            alpha = 0.8,
            bins = bins,
            color = "white",
            size = 0.2
          ) +
          ggplot2::labs(title = target_var, x = NULL, y = NULL) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            plot.title = ggplot2::element_text(
              size = 10,
              face = "bold",
              hjust = 0.5,
              margin = ggplot2::margin(b = 5, unit = "pt")
            ),
            axis.text = ggplot2::element_text(size = 8),
            axis.title = ggplot2::element_text(size = 9),
            panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.2),
            panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(
              color = "grey80",
              fill = NA,
              linewidth = 0.3
            ),
            plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "pt"),
            panel.background = ggplot2::element_rect(fill = "white", color = NA)
          )
        
        plot_list[[i]] <- p
      }
      
      # Arrange all plots in grid
      combined_plot <- gridExtra::grid.arrange(
        grobs = plot_list,
        ncol = n_cols,
        nrow = n_rows,
        top = paste("Histogram Overview:", n_plots, "Variables")
      )
      
      # Save the overview plot
      ggplot2::ggsave(
        filename = output_file,
        plot = combined_plot,
        width = width,
        height = height,
        dpi = dpi,
        bg = "white",
        limitsize = FALSE
      )
      
      # Print confirmation message
      cat("\n✅ Histogram overview successfully created!\n")
      cat("📁 Output file:", output_file, "\n")
      cat("📐 Final resolution:", round(width, 1), "×", round(height, 1), "inches\n")
      cat("🖼️  DPI:", dpi, "\n")
      cat("📊 Number of histograms:", n_plots, "\n")
      cat("🔢 Grid layout:", n_rows, "rows ×", n_cols, "columns\n")
      cat("📦 Number of bins:", bins, "\n")
      
      return(invisible(combined_plot))
    }
    
,

create_scatter_matrix = function(data, var_names, output_path,  
                                 correlation_method = "pearson",
                                 resolution = 300, width = 8, height = 6,
                                 axis_text_size = 5,        # For ALL axis text
                                 axis_title_size = 9,       # Axis titles
                                 strip_text_size = 5,       # Facet labels
                                 dot_size = 0.5,           # Dot size
                                 cor_text_size = 2.5,      # Correlation text
                                 cor_text_color = "black", # Correlation text color
                                 cor_digits = 2,           # Decimal places for correlation
                                 cor_use = "complete.obs") { # Handling of missing data
  
  # Load required libraries
  required_packages <- c("GGally", "ggplot2")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop("Package ", pkg, " is required. Please install it.")
    }
  }
  
  # Validate correlation method
  valid_methods <- c("pearson", "spearman", "kendall")
  if (!correlation_method %in% valid_methods) {
    stop("Invalid correlation_method. Choose from: ", 
         paste(valid_methods, collapse = ", "))
  }
  
  # Validate inputs
  if (!all(var_names %in% names(data))) {
    missing_vars <- setdiff(var_names, names(data))
    stop("The following variables were not found in the data: ", 
         paste(missing_vars, collapse = ", "))
  }
  
  # Custom correlation function with specified method
  custom_cor <- function(data, mapping, 
                         method = correlation_method,
                         use = cor_use,
                         digits = cor_digits,
                         color = cor_text_color,
                         size = cor_text_size,
                         ...) {
    
    # Extract data using GGally's helper function
    x_data <- GGally::eval_data_col(data, mapping$x)
    y_data <- GGally::eval_data_col(data, mapping$y)
    
    # Calculate correlation
    cor_result <- cor(x_data, y_data, 
                      method = method, 
                      use = use)
    
    # Format correlation text
    cor_text <- format(round(cor_result, digits), nsmall = digits)
    
    # Create label with correlation coefficient and method abbreviation
    method_abbr <- switch(method,
                          "pearson" = "r",
                          "spearman" = "ρ",
                          "kendall" = "τ")
    
    label <- paste0(method_abbr, " = ", cor_text)
    
    # Create the plot
    ggplot2::ggplot(data = data, mapping = mapping) +
      ggplot2::annotate("text", 
                        x = mean(range(x_data, na.rm = TRUE)), 
                        y = mean(range(y_data, na.rm = TRUE)), 
                        label = label,
                        color = color, 
                        size = size,
                        fontface = "bold") +
      ggplot2::theme_void()
  }
  
  # Customize scatterplot points
  custom_points <- function(data, mapping, ...) {
    ggplot2::ggplot(data = data, mapping = mapping) + 
      ggplot2::geom_point(..., alpha = 0.6) +
      ggplot2::theme_minimal()
  }
  
  # Create the pairs plot with consistent tick sizes
  p <- GGally::ggpairs(
    data = data[, var_names, drop = FALSE],
    progress = FALSE,
    # Customize the lower triangle (scatterplots)
    lower = list(
      continuous = GGally::wrap(custom_points, size = dot_size)
    ),
    # Customize upper triangle (correlations)
    upper = list(
      continuous = GGally::wrap(custom_cor, 
                                method = correlation_method,
                                use = cor_use,
                                digits = cor_digits,
                                color = cor_text_color,
                                size = cor_text_size)
    ),
    # Customize diagonal (density plots)
    diag = list(
      continuous = GGally::wrap("densityDiag", alpha = 0.7, color = "black", fill = "gray80")
    ),
    # Axis labels appearance
    axisLabels = "show"
  ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      # Single control for ALL axis text (left, right, bottom, top ticks)
      axis.text = ggplot2::element_text(size = axis_text_size),
      
      # Axis titles (x and y axis labels)
      axis.title = ggplot2::element_text(size = axis_title_size),
      
      # Facet/strip labels - FIXED: use ggplot2::margin() with proper namespace
      strip.text = ggplot2::element_text(
        size = strip_text_size,
        face = "bold",
        margin = ggplot2::margin(t = 2, b = 2)  # Fixed this line
      ),
      
      # Panel grid settings
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.1),
      
      # Ensure consistent appearance across all panels
      panel.spacing = grid::unit(0.5, "lines"),
      
      # Panel background
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  
  # Add title with correlation method info
  p <- p + ggplot2::ggtitle(paste("Scatter Plot Matrix -", 
                                  tools::toTitleCase(correlation_method), 
                                  "Correlation"))
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Save the plot with inch dimensions
  grDevices::png(
    filename = output_path,
    width = width, 
    height = height, 
    units = "in",
    res = resolution
  )
  print(p)
  grDevices::dev.off()
  
  # Informative message
  message("Scatter matrix saved to: ", output_path)
  message("Dimensions: ", width, " x ", height, " inches")
  message("Resolution: ", resolution, " dpi")
  message("Correlation method: ", correlation_method)
  message("All axis ticks (left/right/bottom/top): ", axis_text_size, "pt")
  message("Axis titles: ", axis_title_size, "pt")
  message("Facet labels: ", strip_text_size, "pt")
  message("Correlation decimal places: ", cor_digits)
  
  return(invisible(p))
}


,
create_correlation_plot = function(data, var_names, output_file, 
                                    width_inches = 8, height_inches = 8, 
                                    dpi = 300) {
  # Load required library
  if (!require(corrplot)) {
    stop("corrplot package is required. Please install it using install.packages('corrplot')")
  }
  
  # Validate inputs
  if (!all(var_names %in% names(data))) {
    missing_vars <- setdiff(var_names, names(data))
    stop("The following variables were not found in the data: ", 
         paste(missing_vars, collapse = ", "))
  }
  
  # Subset data to selected variables
  plot_data <- data[, var_names, drop = FALSE]
  
  # Calculate correlation matrix
  cor_matrix <- cor(plot_data, use = "complete.obs")
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Convert inches to pixels for PNG
  width_pixels <- width_inches * dpi
  height_pixels <- height_inches * dpi
  
  # Create the plot
  grDevices::png(
    filename = output_file,
    width = width_pixels,
    height = height_pixels,
    res = dpi
  )
  
  # Create correlation plot
  corrplot::corrplot(cor_matrix, 
                     method = "color", 
                     type = "upper", 
                     order = "hclust", 
                     tl.cex = 0.7, 
                     tl.col = "black",
                     addCoef.col = "black", 
                     number.cex = 0.6)
  
  grDevices::dev.off()
  
  message("Correlation plot saved to: ", output_file)
  message("Dimensions: ", width_inches, " x ", height_inches, " inches")
  message("Resolution: ", dpi, " dpi")
  message("Variables included: ", paste(var_names, collapse = ", "))
  message("Correlation matrix dimensions: ", nrow(cor_matrix), " x ", ncol(cor_matrix))
  
  return(invisible(cor_matrix))
}
  ,
create_correlation_plot_sig = function(data, var_names, output_file, 
                                       correlation_method = "pearson",
                                       width_inches = 8, height_inches = 8, 
                                       dpi = 300, alpha = 0.05,
                                       order_method = "hclust",
                                       insig_method = "pch",
                                       pch_color = "green",
                                       pch_size = 3,
                                       pch_symbol = 4) {
  
  # Load required libraries
  if (!require(corrplot, quietly = TRUE)) {
    stop("corrplot package is required. Please install it using install.packages('corrplot')")
  }
  
  # Validate correlation method
  valid_methods <- c("pearson", "spearman", "kendall")
  if (!correlation_method %in% valid_methods) {
    stop("Invalid correlation_method. Choose from: ", 
         paste(valid_methods, collapse = ", "))
  }
  
  # Validate order method
  valid_order_methods <- c("original", "AOE", "FPC", "hclust", "alphabet")
  if (!order_method %in% valid_order_methods) {
    stop("Invalid order_method. Choose from: ", 
         paste(valid_order_methods, collapse = ", "))
  }
  
  # Validate inputs
  if (!all(var_names %in% names(data))) {
    missing_vars <- setdiff(var_names, names(data))
    stop("The following variables were not found in the data: ", 
         paste(missing_vars, collapse = ", "))
  }
  
  # Subset data to selected variables
  plot_data <- data[, var_names, drop = FALSE]
  
  # Remove non-numeric variables with warning
  numeric_vars <- sapply(plot_data, is.numeric)
  if (!all(numeric_vars)) {
    non_numeric <- var_names[!numeric_vars]
    warning("The following variables are not numeric and will be excluded: ",
            paste(non_numeric, collapse = ", "))
    var_names <- var_names[numeric_vars]
    plot_data <- plot_data[, numeric_vars, drop = FALSE]
  }
  
  if (length(var_names) < 2) {
    stop("At least 2 numeric variables are required for correlation analysis")
  }
  
  # Calculate correlation matrix with specified method
  cor_matrix <- cor(plot_data, 
                    method = correlation_method, 
                    use = "complete.obs")
  
  # Function to calculate correlation p-values with specified method
  cor.mtest <- function(mat, method = correlation_method, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], method = method, ...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  
  # Calculate p-value matrix with specified method
  p_matrix <- cor.mtest(plot_data, method = correlation_method)
  
  # Create significance indicator matrix (TRUE = significant, FALSE = not significant)
  sig_matrix <- p_matrix <= alpha
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Convert inches to pixels for PNG
  width_pixels <- width_inches * dpi
  height_pixels <- height_inches * dpi
  
  # Get correlation method symbol for title
  cor_symbol <- switch(correlation_method,
                       "pearson" = "r",
                       "spearman" = "ρ", 
                       "kendall" = "τ")
  
  # Create the plot
  grDevices::png(
    filename = output_file,
    width = width_pixels,
    height = height_pixels,
    res = dpi
  )
  
  # Set up layout for title
  layout_matrix <- matrix(c(1, 2), nrow = 2, ncol = 1)
  layout(layout_matrix, heights = c(0.1, 0.9))
  
  # Add title
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(0.5, 0.5, 
       paste("Correlation Matrix (", tools::toTitleCase(correlation_method), " ", cor_symbol, ")", sep = ""),
       cex = 1.5, font = 2)
  
  # Create correlation plot with significance-based coloring
  par(mar = c(0, 0, 1, 0))
  corrplot::corrplot(cor_matrix, 
                     method = "color", 
                     type = "upper", 
                     order = order_method, 
                     tl.cex = 0.7, 
                     tl.col = "black",
                     addCoef.col = "black", 
                     number.cex = 0.6,
                     p.mat = p_matrix,           # Add p-value matrix
                     sig.level = alpha,          # Significance level
                     insig = insig_method,       # Method for non-significant correlations
                     pch.col = pch_color,        # Color for non-significant symbols
                     pch.cex = pch_size,         # Size of symbols
                     pch = pch_symbol)           # Symbol type
  
  grDevices::dev.off()
  
  # Calculate summary statistics
  total_correlations <- sum(!is.na(cor_matrix[upper.tri(cor_matrix)]))
  significant_correlations <- sum(sig_matrix[upper.tri(sig_matrix)], na.rm = TRUE)
  
  # Get descriptive statistics for correlations
  cor_values <- cor_matrix[upper.tri(cor_matrix)]
  cor_values <- cor_values[!is.na(cor_values)]
  
  message("Correlation plot saved to: ", output_file)
  message("Dimensions: ", width_inches, " x ", height_inches, " inches")
  message("Resolution: ", dpi, " dpi")
  message("Correlation method: ", correlation_method, " (", cor_symbol, ")")
  message("Variables included: ", paste(var_names, collapse = ", "))
  message("Correlation matrix dimensions: ", nrow(cor_matrix), " x ", ncol(cor_matrix))
  message("Significance level: α = ", alpha)
  message("Significant correlations: ", significant_correlations, " out of ", total_correlations, 
          " (", round(significant_correlations/total_correlations * 100, 1), "%)")
  message("Correlation statistics:")
  message("  Range: [", round(min(cor_values), 3), ", ", round(max(cor_values), 3), "]")
  message("  Mean |r|: ", round(mean(abs(cor_values)), 3))
  message("  Median: ", round(median(cor_values), 3))
  
  return(invisible(list(
    correlation_matrix = cor_matrix,
    p_value_matrix = p_matrix,
    significance_matrix = sig_matrix,
    method = correlation_method,
    summary = list(
      total_correlations = total_correlations,
      significant_correlations = significant_correlations,
      percent_significant = round(significant_correlations/total_correlations * 100, 1),
      correlation_range = range(cor_values),
      mean_absolute_correlation = mean(abs(cor_values))
    )
  )))
} 


,


regression_models = function(data, var_names, threshold = 0.7, cor_method = "pearson") {
  # Load required libraries
  if (!require(dplyr)) {
    stop("dplyr package is required. Please install it using install.packages('dplyr')")
  }
  
  # Validate inputs
  if (!all(var_names %in% names(data))) {
    missing_vars <- setdiff(var_names, names(data))
    stop("The following variables were not found in the data: ", 
         paste(missing_vars, collapse = ", "))
  }
  
  if (length(var_names) < 2) {
    stop("At least 2 variables are required for correlation analysis.")
  }
  
  # Validate correlation method
  valid_methods <- c("pearson", "kendall", "spearman")
  if (!cor_method %in% valid_methods) {
    stop("Invalid correlation method. Choose from: ", paste(valid_methods, collapse = ", "))
  }
  
  # Subset data to selected variables
  plot_data <- data[, var_names, drop = FALSE]
  
  # Calculate correlation matrix with specified method
  cor_matrix <- cor(plot_data, use = "complete.obs", method = cor_method)
  
  # Initialize results data frame
  results <- data.frame(
    var1 = character(),
    var2 = character(),
    method = character(),
    r = numeric(),
    p_value_r = numeric(),
    a = numeric(),
    b = numeric(),
    p_value_intercept = numeric(),
    p_value_slope = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through all variable pairs
  for (i in 1:(length(var_names) - 1)) {
    for (j in (i + 1):length(var_names)) {
      var1 <- var_names[i]
      var2 <- var_names[j]
      
      # Remove rows with NA in either variable
      clean_data <- plot_data[complete.cases(plot_data[c(var1, var2)]), ]
      
      if (nrow(clean_data) < 3) {
        warning(sprintf("Insufficient data for pair %s ~ %s. Skipping.", var1, var2))
        next
      }
      
      # Get correlation and p-value with specified method
      cor_test <- cor.test(clean_data[[var1]], clean_data[[var2]], method = cor_method)
      r_value <- cor_test$estimate
      p_value_r <- cor_test$p.value
      
      # Check if absolute correlation meets threshold
      if (abs(r_value) >= threshold) {
        # Perform linear regression using data frame indexing to handle special characters
        lm_model <- lm(clean_data[[var2]] ~ clean_data[[var1]])
        lm_summary <- summary(lm_model)
        
        # Extract coefficients and p-values
        a <- coef(lm_model)[1]  # Intercept
        b <- coef(lm_model)[2]  # Slope
        p_value_intercept <- coef(lm_summary)[1, 4]  # p-value for intercept
        p_value_slope <- coef(lm_summary)[2, 4]      # p-value for slope
        
        # Add to results
        results <- rbind(results, data.frame(
          var1 = var1,
          var2 = var2,
          method = cor_method,
          r = round(r_value, 4),
          p_value_r = round(p_value_r, 6),
          a = round(a, 4),
          b = round(b, 4),
          p_value_intercept = round(p_value_intercept, 6),
          p_value_slope = round(p_value_slope, 6),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Return results sorted by absolute correlation (descending)
  if (nrow(results) > 0) {
    results <- results[order(-abs(results$r)), ]
    rownames(results) <- NULL
  }
  
  # Print summary
  message("Analysis complete!")
  message("Variables analyzed: ", paste(var_names, collapse = ", "))
  message("Correlation method: ", cor_method)
  message("Correlation threshold: |r| >= ", threshold)
  message("Strong correlations found: ", nrow(results))
  
  if (nrow(results) > 0) {
    message("\nStrong correlation pairs:")
    for (i in 1:nrow(results)) {
      message(sprintf("  %s ~ %s: r = %.3f", 
                      results$var1[i], results$var2[i], results$r[i]))
    }
  } else {
    message("No correlations meeting the threshold were found.")
  }
  
  return(results)
}


,
create_regression_plot = function(data, var_x, var_y, 
                                  output_file, width_inches = 8, height_inches = 6, 
                                  dpi = 300) {
  # Load required libraries
  if (!require(ggplot2)) {
    stop("ggplot2 package is required. Please install it using install.packages('ggplot2')")
  }
  
  # Validate inputs
  if (!var_x %in% names(data)) {
    stop("Variable '", var_x, "' was not found in the data")
  }
  if (!var_y %in% names(data)) {
    stop("Variable '", var_y, "' was not found in the data")
  }
  
  # Remove rows with NA values in the selected variables
  plot_data <- data[complete.cases(data[c(var_x, var_y)]), ]
  
  if (nrow(plot_data) == 0) {
    stop("No complete cases found for the specified variables")
  }
  
  if (nrow(plot_data) < 2) {
    stop("At least 2 complete cases are required for regression")
  }
  
  # Calculate linear regression
  lm_model <- lm(plot_data[[var_y]] ~ plot_data[[var_x]])
  lm_summary <- summary(lm_model)
  
  # Extract coefficients
  intercept <- coef(lm_model)[1]
  slope <- coef(lm_model)[2]
  
  # Extract statistics
  r_squared <- lm_summary$r.squared
  p_value <- coef(lm_summary)[2, 4]
  correlation <- cor(plot_data[[var_x]], plot_data[[var_y]])
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir) && output_dir != ".") {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create the regression line data
  x_range <- range(plot_data[[var_x]], na.rm = TRUE)
  x_seq <- seq(x_range[1], x_range[2], length.out = 100)
  y_pred <- intercept + slope * x_seq
  line_data <- data.frame(x = x_seq, y = y_pred)
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = .data[[var_x]], y = .data[[var_y]])) +
    geom_point(alpha = 0.6, color = "blue", size = 2) +
    geom_line(data = line_data, aes(x = x, y = y), 
              color = "red", linewidth = 1.2, linetype = "solid") +
    labs(
      title = paste("Scatter Plot with Regression Line:", var_y, "vs", var_x),
      x = var_x,
      y = var_y,
      subtitle = paste("y =", round(intercept, 4), "+", round(slope, 4), "* x",
                       "| R² =", round(r_squared, 4), "| p =", round(p_value, 6))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  # Save the plot
  ggsave(
    filename = output_file,
    plot = p,
    width = width_inches,
    height = height_inches,
    dpi = dpi
  )
  
  # Print summary information
  message("Regression plot saved to: ", output_file)
  message("Dimensions: ", width_inches, " x ", height_inches, " inches")
  message("Resolution: ", dpi, " dpi")
  message("Variables: ", var_y, " ~ ", var_x)
  message("Regression equation: y = ", round(intercept, 4), " + ", round(slope, 4), " * x")
  message("R-squared: ", round(r_squared, 4))
  message("P-value: ", round(p_value, 6))
  message("Correlation (r): ", round(correlation, 4))
  message("Data points plotted: ", nrow(plot_data))
  
  return(invisible(list(
    plot = p,
    intercept = intercept,
    slope = slope,
    r_squared = r_squared,
    p_value = p_value,
    correlation = correlation,
    n_observations = nrow(plot_data),
    lm_model = lm_model
  )))
}
,
shapiro_test_comprehensive = function(data, variables, alpha = 0.05) {
  # Load required packages
  if (!require(dplyr)) {
    stop("dplyr package is required. Please install it using install.packages('dplyr')")
  }
  
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }
  
  if (missing(variables)) {
    variables <- names(data)[sapply(data, is.numeric)]
    message("No variables specified. Using all numeric variables: ", 
            paste(variables, collapse = ", "))
  }
  
  # Check if variables exist
  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars) > 0) {
    stop("The following variables were not found in the data: ", 
         paste(missing_vars, collapse = ", "))
  }
  
  # Check if variables are numeric
  non_numeric_vars <- variables[!sapply(data[variables], is.numeric)]
  if (length(non_numeric_vars) > 0) {
    warning("The following variables are not numeric and will be skipped: ", 
            paste(non_numeric_vars, collapse = ", "))
    variables <- setdiff(variables, non_numeric_vars)
  }
  
  if (length(variables) == 0) {
    stop("No numeric variables to test")
  }
  
  # Perform tests and create results
  results_list <- lapply(variables, function(var) {
    clean_data <- na.omit(data[[var]])
    n <- length(clean_data)
    
    if (n < 3) {
      return(data.frame(
        variable = var,
        statistic = NA,
        p_value = NA,
        normality = "Insufficient data",
        mean = round(mean(clean_data), 4),
        sd = round(sd(clean_data), 4),
        skewness = ifelse(n >= 3, round(moments::skewness(clean_data), 4), NA),
        kurtosis = ifelse(n >= 3, round(moments::kurtosis(clean_data), 4), NA),
        stringsAsFactors = FALSE
      ))
    }
    
    if (n > 5000) {
      warning("Variable '", var, "' has more than 5000 observations. ",
              "Shapiro-Wilk test may not be reliable for large samples.")
    }
    
    shapiro_result <- shapiro.test(clean_data)
    
    # Calculate additional statistics
    if (require(moments, quietly = TRUE)) {
      skew_val <- moments::skewness(clean_data)
      kurt_val <- moments::kurtosis(clean_data)
    } else {
      skew_val <- NA
      kurt_val <- NA
      message("Install 'moments' package for skewness and kurtosis calculations")
    }
    
    data.frame(
      variable = var,
      W = round(shapiro_result$statistic, 4),
      p_value = round(shapiro_result$p.value, 6),
      normality = ifelse(shapiro_result$p.value >= alpha, "Normal", "Non-normal"),
      mean = round(mean(clean_data), 4),
      sd = round(sd(clean_data), 4),
      skewness = round(skew_val, 4),
      kurtosis = round(kurt_val, 4),
      stringsAsFactors = FALSE
    )
  })
  
  # Combine all results
  results <- do.call(rbind, results_list)
  
  # Sort by p-value
  results <- results[order(results$p_value, na.last = TRUE), ]
  rownames(results) <- NULL
  
  # Print summary
  cat("=== Shapiro-Wilk Normality Test Results ===\n")
  cat("Significance level: α =", alpha, "\n")
  cat("Variables tested:", length(variables), "\n")
  cat("Normal distributions:", sum(results$normality == "Normal", na.rm = TRUE), "\n")
  cat("Non-normal distributions:", sum(results$normality == "Non-normal", na.rm = TRUE), "\n")
  cat("Insufficient data:", sum(results$normality == "Insufficient data", na.rm = TRUE), "\n")
  
  return(results)
}
,
kruskal_wallis_comprehensive = function(data, iv_vars, dv_vars, 
                                         alpha = 0.05, posthoc = FALSE,
                                         effect_size = FALSE) {
  # Load required packages
  if (!require(dplyr, quietly = TRUE)) {
    stop("dplyr package is required")
  }
  
  if (posthoc && !require(PMCMRplus, quietly = TRUE)) {
    warning("PMCMRplus package not available for post-hoc tests. Install with: install.packages('PMCMRplus')")
    posthoc <- FALSE
  }
  
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }
  
  if (missing(iv_vars) || missing(dv_vars)) {
    stop("Both IV (independent) and DV (dependent) variables must be specified")
  }
  
  # Check if variables exist
  missing_iv <- setdiff(iv_vars, names(data))
  missing_dv <- setdiff(dv_vars, names(data))
  
  if (length(missing_iv) > 0) {
    stop("The following IV variables were not found: ", paste(missing_iv, collapse = ", "))
  }
  if (length(missing_dv) > 0) {
    stop("The following DV variables were not found: ", paste(missing_dv, collapse = ", "))
  }
  
  # Convert IV variables to factors if needed
  for (iv in iv_vars) {
    if (!is.factor(data[[iv]])) {
      data[[iv]] <- as.factor(data[[iv]])
      message("Converted IV variable '", iv, "' to factor")
    }
  }
  
  # Initialize empty data frame with proper column structure
  results <- data.frame(
    iv_variable = character(),
    dv_variable = character(),
    chi_squared = numeric(),
    df = integer(),
    p_value = numeric(),
    n_groups = integer(),
    n_obs = integer(),
    significance = character(),
    stringsAsFactors = FALSE
  )
  
  # Add optional columns if requested
  if (effect_size) {
    results$epsilon_squared <- numeric()
  }
  
  if (posthoc) {
    results$posthoc_significant_pairs <- character()
  }
  
  # Perform Kruskal-Wallis test for each combination
  for (iv in iv_vars) {
    for (dv in dv_vars) {
      # Remove NA values for current combination - handle variable names with spaces
      clean_data <- data[complete.cases(data[c(iv, dv)]), ]
      
      if (nrow(clean_data) == 0) {
        warning("No complete cases for combination: ", iv, " ~ ", dv)
        next
      }
      
      # Check if we have at least 2 groups with data
      group_counts <- table(clean_data[[iv]])
      valid_groups <- sum(group_counts > 0)
      
      if (valid_groups < 2) {
        warning("Need at least 2 groups with data for ", iv, " ~ ", dv, ". Skipping.")
        next
      }
      
      # Perform Kruskal-Wallis test - handle variables with spaces using [[ ]]
      kruskal_test <- tryCatch({
        kruskal.test(x = clean_data[[dv]], g = clean_data[[iv]])
      }, error = function(e) {
        warning("Kruskal-Wallis test failed for ", iv, " ~ ", dv, ": ", e$message)
        return(NULL)
      })
      
      if (is.null(kruskal_test)) {
        next
      }
      
      # Calculate effect size (epsilon squared)
      epsilon_squared <- NA
      if (effect_size) {
        epsilon_squared <- (kruskal_test$statistic - kruskal_test$parameter + 1) / 
          (nrow(clean_data) - valid_groups)
      }
      
      # Determine significance
      sig_star <- ifelse(kruskal_test$p.value < 0.001, "***",
                         ifelse(kruskal_test$p.value < 0.01, "**",
                                ifelse(kruskal_test$p.value < 0.05, "*", "ns")))
      
      # Initialize posthoc results
      posthoc_sig_pairs <- NA_character_
      
      # Perform post-hoc tests if requested and test is significant
      if (posthoc && kruskal_test$p.value < alpha && valid_groups > 2) {
        tryCatch({
          posthoc_results <- PMCMRplus::kwAllPairsNemenyiTest(
            x = clean_data[[dv]], 
            g = clean_data[[iv]]  # Handle variables with spaces
          )
          
          # Extract significant pairs from posthoc results
          posthoc_matrix <- posthoc_results$p.value
          sig_pairs <- which(posthoc_matrix < alpha, arr.ind = TRUE)
          
          if (nrow(sig_pairs) > 0) {
            pair_names <- apply(sig_pairs, 1, function(idx) {
              row_name <- rownames(posthoc_matrix)[idx[1]]
              col_name <- colnames(posthoc_matrix)[idx[2]]
              p_val <- posthoc_matrix[idx[1], idx[2]]
              paste0(row_name, "-", col_name, " (p=", round(p_val, 4), ")")
            })
            posthoc_sig_pairs <- paste(pair_names, collapse = "; ")
          }
        }, error = function(e) {
          warning("Post-hoc test failed for ", iv, " ~ ", dv, ": ", e$message)
        })
      }
      
      # Create new row as a list
      new_row <- list(
        iv_variable = iv,
        dv_variable = dv,
        chi_squared = round(kruskal_test$statistic, 4),
        df = as.integer(kruskal_test$parameter),
        p_value = round(kruskal_test$p.value, 6),
        n_groups = as.integer(valid_groups),
        n_obs = as.integer(nrow(clean_data)),
        significance = sig_star
      )
      
      # Add optional columns
      if (effect_size) {
        new_row$epsilon_squared <- round(epsilon_squared, 4)
      }
      
      if (posthoc) {
        new_row$posthoc_significant_pairs <- posthoc_sig_pairs
      }
      
      # Convert to data frame and bind rows
      new_row_df <- as.data.frame(new_row, stringsAsFactors = FALSE)
      results <- dplyr::bind_rows(results, new_row_df)
    }
  }
  
  # Check if we have any results
  if (nrow(results) == 0) {
    stop("No valid results were produced. Check your data and variable combinations.")
  }
  
  # Sort by p-value
  results <- results[order(results$p_value), ]
  rownames(results) <- NULL
  
  # Print summary
  cat("=== Kruskal-Wallis Analysis Summary ===\n")
  cat("Total combinations tested:", nrow(results), "\n")
  cat("Significant results (p <", alpha, "):", sum(results$p_value < alpha), "\n")
  cat("Significance codes: *** p < 0.001, ** p < 0.01, * p < 0.05, ns not significant\n")
  
  if (effect_size) {
    cat("\nEffect size interpretation (epsilon squared):\n")
    cat("0.01 < 0.06 (small), 0.06 < 0.14 (medium), ≥ 0.14 (large)\n")
  }
  
  # Return the data frame
  return(results)
}
,
dunn_tests = function(data, dv_names, iv_names, 
                      method = "bonferroni", 
                      alpha = 0.05) {
  # Initialize empty list to store results
  all_results <- list()
  result_counter <- 1
  
  # Helper function to extract Dunn's test results
  extract_dunn_results <- function(dunn_output, dv, iv, kw_pval) {
    comparisons <- strsplit(dunn_output$comparisons, " - ")
    
    # Handle case where no comparisons are returned
    if(length(comparisons) == 0) {
      return(data.frame(
        dependent_var = character(0),
        independent_var = character(0),
        group1 = character(0),
        group2 = character(0),
        kw_statistic = numeric(0),
        kw_p_value = numeric(0),
        z_value = numeric(0),
        p_value = numeric(0),
        p_adjusted = numeric(0),
        significant = logical(0),
        stringsAsFactors = FALSE
      ))
    }
    
    data.frame(
      dependent_var = dv,
      independent_var = iv,
      group1 = sapply(comparisons, function(x) x[1]),
      group2 = sapply(comparisons, function(x) x[2]),
      kw_statistic = kw_pval$statistic,
      kw_p_value = kw_pval$p.value,
      z_value = dunn_output$Z,
      p_value = dunn_output$P,
      p_adjusted = dunn_output$P.adjusted,
      significant = dunn_output$P.adjusted < alpha,
      stringsAsFactors = FALSE
    )
  }
  
  # Loop through each dependent and independent variable combination
  for(dv in dv_names) {
    for(iv in iv_names) {
      
      # Check if the IV exists and has at least 2 groups
      if(!iv %in% names(data)) next
      
      group_counts <- table(data[[iv]])
      if(length(group_counts) < 2) next
      
      cat("Analyzing", dv, "by", iv, "\n")
      
      # Perform Kruskal-Wallis test
      kw_formula <- as.formula(paste(dv, "~", iv))
      kw_test <- tryCatch({
        kruskal.test(kw_formula, data = data)
      }, error = function(e) {
        cat("Error in Kruskal-Wallis for", dv, "~", iv, ":", e$message, "\n")
        return(NULL)
      })
      
      if(is.null(kw_test)) next
      
      # Store Kruskal-Wallis result (even if not significant)
      kw_result <- data.frame(
        dependent_var = dv,
        independent_var = iv,
        group1 = "Overall",
        group2 = "Overall", 
        kw_statistic = kw_test$statistic,
        kw_p_value = kw_test$p.value,
        z_value = NA_real_,
        p_value = NA_real_,
        p_adjusted = NA_real_,
        significant = kw_test$p.value < alpha,
        stringsAsFactors = FALSE
      )
      
      all_results[[result_counter]] <- kw_result
      result_counter <- result_counter + 1
      
      # If Kruskal-Wallis is significant, perform Dunn's test
      if(kw_test$p.value < alpha) {
        dunn_result <- tryCatch({
          dunn_output <- dunn.test(data[[dv]], data[[iv]], 
                                   method = method, kw = FALSE)
          extract_dunn_results(dunn_output, dv, iv, kw_test)
        }, error = function(e) {
          cat("Error in Dunn's test for", dv, "~", iv, ":", e$message, "\n")
          return(NULL)
        })
        
        if(!is.null(dunn_result) && nrow(dunn_result) > 0) {
          all_results[[result_counter]] <- dunn_result
          result_counter <- result_counter + 1
        }
      }
    }
  }
  
  # Combine all results
  if(length(all_results) == 0) {
    return(data.frame(
      dependent_var = character(0),
      independent_var = character(0),
      group1 = character(0),
      group2 = character(0),
      kw_statistic = numeric(0),
      kw_p_value = numeric(0),
      z_value = numeric(0),
      p_value = numeric(0),
      p_adjusted = numeric(0),
      significant = logical(0),
      stringsAsFactors = FALSE
    ))
  }
  
  bind_rows(all_results)
}
,
permanova = function(data, dv_vars, iv_vars, permutations = 9999) {
  # 1. PERMANOVA for multivariate effects
  lipid_matrix <- data[, dv_vars]
  complete_cases <- complete.cases(lipid_matrix)
  lipid_clean <- lipid_matrix[complete_cases, ]
  
  # Fix: Use row.names = NULL to avoid warnings
  meta_clean <- as.data.frame(data[complete_cases, iv_vars], row.names = NULL)
  
  if (nrow(lipid_clean) > 0) {
    dist_matrix <- vegdist(lipid_clean, method = "bray")
    permanova_formula <- as.formula(paste("dist_matrix ~", paste(iv_vars, collapse = " + ")))
    permanova_result <- adonis2(permanova_formula,
                                data = meta_clean,
                                permutations = permutations)
    print(permanova_result)
    # Convert to data.frame
    result_df <- as.data.frame(permanova_result)
    
    # Add variable names as a column
    result_df$Variable <- rownames(result_df)
    rownames(result_df) <- NULL
    
    # Reorder columns to have Variable first
    result_df <- result_df[, c("Variable", setdiff(names(result_df), "Variable"))]
    
    return(result_df)
  } else {
    # Return empty data.frame with expected structure if no complete cases
    empty_df <- data.frame(
      Variable = character(),
      Df = numeric(),
      SumOfSqs = numeric(),
      F = numeric(),
      R2 = numeric(),
      `Pr(>F)` = numeric(),
      stringsAsFactors = FALSE
    )
    return(empty_df)
  }
}
,
hierarchical_clustering_plot = function(data_frame, 
                                        var_names, 
                                        data_label_var = NULL,
                                        distance_method = "euclidean",
                                        output_file = "hclust_plot.png",
                                        width = 10,
                                        height = 8,
                                        dpi = 300,
                                        linkage_method = "ward.D2",
                                        num_clusters = NULL,
                                        color_palette = NULL,
                                        label_font_size = 0.7,
                                        horizontal = FALSE) {
  
  # Load required packages
  if (!requireNamespace("dendextend", quietly = TRUE)) {
    install.packages("dendextend")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  if (!requireNamespace("factoextra", quietly = TRUE)) {
    install.packages("factoextra")
  }
  
  library(dendextend)
  library(ggplot2)
  library(factoextra)
  
  # Input validation
  if (!is.data.frame(data_frame)) {
    stop("data_frame must be a data frame")
  }
  
  if (!all(var_names %in% names(data_frame))) {
    missing_vars <- setdiff(var_names, names(data_frame))
    stop("The following variables are not in the data frame: ", 
         paste(missing_vars, collapse = ", "))
  }
  
  if (!is.null(data_label_var) && !data_label_var %in% names(data_frame)) {
    stop("data_label_var '", data_label_var, "' not found in data frame")
  }
  
  if (!is.numeric(label_font_size) || label_font_size <= 0) {
    stop("label_font_size must be a positive number")
  }
  
  if (!is.logical(horizontal)) {
    stop("horizontal must be a logical value (TRUE or FALSE)")
  }
  
  # Prepare data for clustering
  clustering_data <- data_frame[, var_names, drop = FALSE]
  
  # Remove rows with missing values
  complete_cases <- complete.cases(clustering_data)
  clustering_data <- clustering_data[complete_cases, ]
  
  if (nrow(clustering_data) == 0) {
    stop("No complete cases available after removing missing values")
  }
  
  # Standardize the data
  scaled_data <- scale(clustering_data)
  
  # Calculate distance matrix
  dist_matrix <- dist(scaled_data, method = distance_method)
  
  # Perform hierarchical clustering
  hc <- hclust(dist_matrix, method = linkage_method)
  
  # Determine optimal number of clusters if not specified
  if (is.null(num_clusters)) {
    # Use silhouette method to suggest optimal clusters
    sil_width <- fviz_nbclust(scaled_data, FUN = hcut, method = "silhouette", 
                              k.max = min(10, nrow(scaled_data) - 1))
    num_clusters <- which.max(sil_width$data$y)
    message("Optimal number of clusters suggested by silhouette method: ", num_clusters)
  }
  
  # Set color palette
  if (is.null(color_palette)) {
    color_palette <- c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07", "#FF0066", 
                       "#00CC99", "#9966FF", "#FF9933", "#66FF33", "#FF6666")
  }
  
  # Create labels for dendrogram
  if (!is.null(data_label_var)) {
    labels <- as.character(data_frame[complete_cases, data_label_var])
  } else {
    labels <- rownames(clustering_data)
    if (is.null(labels)) {
      labels <- paste0("Obs", 1:nrow(clustering_data))
    }
  }
  
  # Create and customize dendrogram
  dend <- as.dendrogram(hc)
  
  # Color branches and labels
  dend <- color_branches(dend, k = num_clusters, palette = color_palette)
  dend <- color_labels(dend, k = num_clusters, palette = color_palette)
  
  # Set labels and font size
  labels(dend) <- labels[order.dendrogram(dend)]
  
  # Start PNG device
  png(filename = output_file, width = width, height = height, 
      units = "in", res = dpi, type = "cairo")
  
  # Adjust margins based on orientation and label font size
  if (horizontal) {
    # Horizontal plot - more space on left for labels
    left_margin <- max(6, 4 + (label_font_size * 3))
    par(mar = c(4, left_margin, 4, 2) + 0.1)
  } else {
    # Vertical plot - more space on bottom for labels
    bottom_margin <- max(8, 6 + (label_font_size * 2))
    par(mar = c(bottom_margin, 4, 4, 2) + 0.1)
  }
  
  # Plot title
  plot_title <- paste("Hierarchical Clustering\n",
                      "Method:", linkage_method, 
                      "| Distance:", distance_method,
                      "| Variables:", length(var_names))
  
  # Plot dendrogram with orientation
  if (horizontal) {
    plot(dend, 
         main = plot_title,
         xlab = "Height",
         ylab = "",
         cex = label_font_size,
         horiz = TRUE,
         axes = TRUE)
    
    # Add cluster rectangles for horizontal plot
    rect_data <- rect.dendrogram(dend, k = num_clusters, 
                                 border = head(color_palette, num_clusters),
                                 horiz = TRUE)
  } else {
    plot(dend, 
         main = plot_title,
         ylab = "Height",
         xlab = "",
         cex = label_font_size,
         horiz = FALSE)
    
    # Add cluster rectangles for vertical plot
    rect.hclust(hc, k = num_clusters, border = head(color_palette, num_clusters))
  }
  
  # Add legend with adjusted position based on orientation
  legend_font_size <- max(0.6, label_font_size * 0.9)
  
  if (horizontal) {
    legend("bottomright", 
           legend = paste("Cluster", 1:num_clusters),
           fill = head(color_palette, num_clusters),
           cex = legend_font_size,
           box.lty = 0)
  } else {
    legend("topright", 
           legend = paste("Cluster", 1:num_clusters),
           fill = head(color_palette, num_clusters),
           cex = legend_font_size,
           box.lty = 0)
  }
  
  # Close PNG device
  dev.off()
  
  message("Clustering plot saved to: ", output_file)
  message("Clustering details:")
  message("  - Number of observations: ", nrow(clustering_data))
  message("  - Number of variables: ", length(var_names))
  message("  - Distance method: ", distance_method)
  message("  - Linkage method: ", linkage_method)
  message("  - Number of clusters: ", num_clusters)
  message("  - Label font size: ", label_font_size)
  message("  - Orientation: ", ifelse(horizontal, "Horizontal", "Vertical"))
  
  # Return cluster assignments
  cluster_assignments <- cutree(hc, k = num_clusters)
  
  # Create a summary data frame
  result <- list(
    hclust_object = hc,
    cluster_assignments = cluster_assignments,
    dendrogram = dend,
    num_clusters = num_clusters,
    label_font_size = label_font_size,
    horizontal = horizontal,
    output_file = output_file
  )
  
  return(invisible(result))
}
,
binarize_variables_by_mean = function(data_frame, 
                                       var_names, 
                                       suffix = "_bin", 
                                       na_action = "keep_na") {
  # Input validation
  if (!is.data.frame(data_frame)) {
    stop("data_frame must be a data frame")
  }
  
  if (!is.character(var_names)) {
    stop("var_names must be a character vector")
  }
  
  if (!all(var_names %in% names(data_frame))) {
    missing_vars <- setdiff(var_names, names(data_frame))
    stop("The following variables are not in the data frame: ", 
         paste(missing_vars, collapse = ", "))
  }
  
  valid_na_actions <- c("keep_na", "to_zero", "to_one")
  if (!na_action %in% valid_na_actions) {
    stop("na_action must be one of: ", paste(valid_na_actions, collapse = ", "))
  }
  
  # Create a copy of the original data frame
  result_df <- data_frame
  
  # Initialize a data frame to store means
  means_summary <- data.frame(
    variable = character(),
    mean_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Binarize each variable
  for (var in var_names) {
    # Calculate mean (remove NAs for calculation)
    var_mean <- mean(data_frame[[var]], na.rm = TRUE)
    
    # Store mean in summary
    means_summary <- rbind(means_summary, 
                           data.frame(variable = var, mean_value = var_mean))
    
    # Create binarized variable
    binarized_var <- ifelse(data_frame[[var]] < var_mean, 0, 1)
    
    # Handle NA values based on na_action
    if (na_action == "to_zero") {
      binarized_var[is.na(binarized_var)] <- 0
    } else if (na_action == "to_one") {
      binarized_var[is.na(binarized_var)] <- 1
    }
    # If "keep_na", NAs remain as NA
    
    # Add binarized variable to result data frame
    new_var_name <- paste0(var, suffix)
    result_df[[new_var_name]] <- binarized_var
  }
  
  # Print means summary
  cat("=== MEANS OF VARIABLES ===\n")
  for (i in 1:nrow(means_summary)) {
    cat(sprintf("%-20s: %10.4f\n", 
                means_summary$variable[i], 
                means_summary$mean_value[i]))
  }
  cat("==========================\n\n")
  
  # Print binarization summary
  cat("=== BINARIZATION SUMMARY ===\n")
  cat("Variables binarized:", length(var_names), "\n")
  cat("Suffix used:", suffix, "\n")
  cat("NA handling:", na_action, "\n")
  
  # Show sample of binarized values for first variable
  if (length(var_names) > 0) {
    first_bin_var <- paste0(var_names[1], suffix)
    cat("Sample binarized values for", first_bin_var, ":\n")
    sample_values <- head(result_df[[first_bin_var]])
    print(sample_values)
  }
  cat("============================\n")
  
  # Add means summary as attribute for further use
  attr(result_df, "means_summary") <- means_summary
  
  return(result_df)
}
,
generate_c50_tree = function(data_frame, 
                              input_var_names, 
                              target_var_name, 
                              output_png_file = "c50_tree.png",
                              width_inches = 12, 
                              height_inches = 8, 
                              dpi = 600) {
  
  # Load packages
  required_packages <- c("C50", "rpart", "rpart.plot")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  
  # Input validation with detailed error messages
  if (!target_var_name %in% names(data_frame)) {
    stop("Target variable '", target_var_name, "' not found. Available variables: ", 
         paste(names(data_frame), collapse = ", "))
  }
  
  missing_inputs <- setdiff(input_var_names, names(data_frame))
  if (length(missing_inputs) > 0) {
    stop("Missing input variables: ", paste(missing_inputs, collapse = ", "),
         "\nAvailable variables: ", paste(names(data_frame), collapse = ", "))
  }
  
  # Handle target variable
  data_frame[[target_var_name]] <- as.factor(data_frame[[target_var_name]])
  message("Target variable '", target_var_name, "' converted to factor")
  cat("Target levels:", paste(levels(data_frame[[target_var_name]]), collapse = ", "), "\n")
  
  # Create a function to safely handle any variable name
  make_formula_safe <- function(target, predictors) {
    # Wrap all variable names in backticks to handle special characters
    safe_target <- paste0("`", gsub("`", "\\\\`", target), "`")
    safe_predictors <- sapply(predictors, function(x) paste0("`", gsub("`", "\\\\`", x), "`"))
    
    formula_str <- paste(safe_target, "~", paste(safe_predictors, collapse = " + "))
    return(as.formula(formula_str))
  }
  
  # Create safe formula
  formula <- make_formula_safe(target_var_name, input_var_names)
  cat("Using formula:", deparse(formula), "\n")
  
  # Build models
  c50_model <- C50::C5.0.formula(formula, data = data_frame)
  rpart_model <- rpart::rpart(formula, data = data_frame, method = "class",
                              control = rpart.control(minsplit = 2, cp = 0.005))
  
  # Create visualization with enhanced settings
  png(output_png_file, width = width_inches, height = height_inches, units = "in", res = dpi)
  
  # Try to create a clean plot even with complex variable names
  rpart.plot(rpart_model,
             main = paste("Decision Tree Analysis\nTarget:", target_var_name),
             type = 5,  # More compact style
             extra = 101,  # Show n and percentage
             box.palette = list("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),  # ColorBrewer palette
             shadow.col = "gray",
             nn = FALSE,  # Hide node numbers for cleaner look
             fallen.leaves = FALSE,
             branch = 0.3,
             roundint = FALSE,
             tweak = 1.05,  # Text size adjustment
             space = 0.1,   # More space between nodes
             gap = 0.5)     # Gap between boxes
  
  dev.off()
  
  # Comprehensive output
  cat(paste0("\n", rep("=", 60), "\n"))
  cat("C5.0 DECISION TREE ANALYSIS RESULTS\n")
  cat(paste0("\n", rep("=", 60), "\n"))
  
  # Model performance
  preds <- predict(c50_model, data_frame, type = "class")
  acc <- mean(preds == data_frame[[target_var_name]])
  cat(sprintf("Training Accuracy: %.1f%%\n", acc * 100))
  
  # Variable importance
  imp <- C50::C5imp(c50_model)
  cat("\nTop 10 Most Important Variables:\n")
  print(head(imp[order(-imp$Overall), , drop = FALSE], 10))
  
  # Variables with special characters
  special_vars <- input_var_names[grepl("[[:space:]]|[^a-zA-Z0-9._]", input_var_names)]
  if (length(special_vars) > 0) {
    cat("\nNote: Handled", length(special_vars), "variables with special characters\n")
  }
  
  cat("Output file:", output_png_file, "\n")
  cat(rep("=", 60) , "\n")
  
  return(list(
    c50_model = c50_model,
    rpart_model = rpart_model,
    accuracy = acc,
    importance = imp,
    formula = deparse(formula),
    special_variables = special_vars,
    output_file = output_png_file
  ))
}
,
# Install and load required packages

# Function to get MACCS fingerprint as binary string
get_fingerprint_maccs_binary = function(smiles) {
  # Parse the molecule from SMILES
  molecule <- parse.smiles(smiles)[[1]]
  
  # Generate MACCS fingerprint (166 bits)
  fp <- get.fingerprint(molecule, type = 'maccs')
  
  # Convert to binary representation
  binary_vector <- as.logical(fp@bits)
  
  # Pad to 166 bits if necessary
  if (length(binary_vector) < 166) {
    binary_vector <- c(binary_vector, rep(FALSE, 166 - length(binary_vector)))
  }
  
  # Convert to binary string (1s and 0s)
  binary_string <- paste(as.integer(binary_vector), collapse = "")
  
  return(list(
    binary_string = binary_string,
    binary_vector = binary_vector,
    length = length(binary_vector)
  ))
}
,
random_forest = function(data_frame, 
                         input_var_names, 
                         target_var_name, 
                         output_png_file = "rf_robust.png",
                         width_inches = 12, 
                         height_inches = 8, 
                         dpi = 600,
                         ntree = 500) {
  
  library(randomForest)
  
  # Input validation
  if (!target_var_name %in% names(data_frame)) {
    stop("Target variable '", target_var_name, "' not found")
  }
  
  missing_vars <- setdiff(input_var_names, names(data_frame))
  if (length(missing_vars) > 0) {
    stop("Missing variables: ", paste(missing_vars, collapse = ", "))
  }
  
  # Convert target to factor
  data_frame[[target_var_name]] <- as.factor(data_frame[[target_var_name]])
  cat("Target levels:", paste(levels(data_frame[[target_var_name]]), collapse = ", "), "\n")
  
  # Use data frame interface directly - most reliable for special characters
  cat("Building Random Forest using data frame interface...\n")
  
  # Select only the columns we need
  model_data <- data_frame[, c(target_var_name, input_var_names)]
  
  # Build model using data frame interface (no formula parsing issues)
  rf_model <- randomForest::randomForest(
    model_data[, -1],  # All columns except the first (target)
    model_data[, 1],   # First column is target
    ntree = ntree,
    importance = TRUE,
    do.trace = ifelse(ntree > 100, floor(ntree/10), ntree)
  )
  
  # Create plot
  png(output_png_file, width = width_inches, height = height_inches, units = "in", res = dpi)
  par(mfrow = c(1, 2), mar = c(7, 4, 4, 2) + 0.1)  # Extra bottom margin for labels
  
  # Error plot
  plot(rf_model, main = "Error Rate vs Trees", lwd = 2)
  grid()
  
  # Importance plot
  if (!is.null(rf_model$importance)) {
    imp <- importance(rf_model)
    
    # Get top variables
    n_show <- min(10, nrow(imp))
    imp_sorted <- imp[order(-imp[, "MeanDecreaseAccuracy"]), ]
    imp_top <- head(imp_sorted, n_show)
    
    # Plot with proper labels
    barplot(imp_top[, "MeanDecreaseAccuracy"],
            names.arg = rownames(imp_top),
            las = 2,  # Vertical labels
            cex.names = 0.7,
            main = paste("Top", n_show, "Variable Importance"),
            ylab = "Mean Decrease Accuracy",
            col = rainbow(n_show, alpha = 0.7))
  }
  
  dev.off()
  
  # Results
  cat(sprintf("\nOOB Error: %.3f (Accuracy: %.1f%%)\n", 
              rf_model$err.rate[ntree, 1],
              (1 - rf_model$err.rate[ntree, 1]) * 100))
  cat("Plot saved to:", output_png_file, "\n")
  
  return(rf_model)
}
,
filter_variables_usingXY = function(data, y_vars, x_vars, cor_threshold, pvalue_threshold) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  if (cor_threshold < 0 | cor_threshold > 1) {
    stop("cor_threshold must be between 0 and 1")
  }
  if (pvalue_threshold < 0 | pvalue_threshold > 1) {
    stop("pvalue_threshold must be between 0 and 1")
  }
  
  # Check if all specified variables exist in the data
  all_vars <- c(y_vars, x_vars)
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("The following variables are not found in the data: ", 
         paste(missing_vars, collapse = ", "))
  }
  
  # Initialize lists to track variables
  all_variables <- x_vars
  removed_zero_sd <- character()
  removed_low_correlation <- character()
  remaining_variables <- character()
  
  # Step 1: Remove x variables with zero standard deviation
  for (x_var in x_vars) {
    if (is.numeric(data[[x_var]])) {
      sd_val <- sd(data[[x_var]], na.rm = TRUE)
      if (!is.na(sd_val) && sd_val == 0) {
        removed_zero_sd <- c(removed_zero_sd, x_var)
      } else {
        remaining_variables <- c(remaining_variables, x_var)
      }
    } else {
      # Keep non-numeric variables for now, they'll be removed in correlation check
      remaining_variables <- c(remaining_variables, x_var)
    }
  }
  
  # Step 2: Check correlation and significance for remaining variables
  final_remaining <- character()
  
  for (x_var in remaining_variables) {
    # Skip non-numeric variables for correlation check
    if (!is.numeric(data[[x_var]])) {
      removed_low_correlation <- c(removed_low_correlation, x_var)
      next
    }
    
    variable_passed <- FALSE
    
    # Check correlation with each y variable
    for (y_var in y_vars) {
      if (!is.numeric(data[[y_var]])) {
        next  # Skip non-numeric y variables for correlation
      }
      
      # Remove rows with NA in either variable for correlation test
      complete_cases <- complete.cases(data[[x_var]], data[[y_var]])
      x_vals <- data[[x_var]][complete_cases]
      y_vals <- data[[y_var]][complete_cases]
      
      if (length(x_vals) < 3) {  # Need at least 3 observations for correlation test
        next
      }
      
      # Calculate correlation and p-value
      cor_test <- cor.test(x_vals, y_vals, method = "pearson")
      abs_cor <- abs(cor_test$estimate)
      p_value <- cor_test$p.value
      
      # Check if correlation meets threshold and is significant
      if (abs_cor >= cor_threshold && p_value <= pvalue_threshold) {
        variable_passed <- TRUE
        break  # No need to check other y variables if one passes
      }
    }
    
    if (variable_passed) {
      final_remaining <- c(final_remaining, x_var)
    } else {
      removed_low_correlation <- c(removed_low_correlation, x_var)
    }
  }
  
  # Return results
  list(
    all_variables = all_variables,
    removed_variables = list(
      zero_sd = removed_zero_sd,
      low_correlation = removed_low_correlation
    ),
    remaining_variables = final_remaining,
    summary = list(
      total_variables = length(all_variables),
      removed_zero_sd = length(removed_zero_sd),
      removed_low_correlation = length(removed_low_correlation),
      remaining = length(final_remaining)
    )
  )
}
,
KNN_label=function(train_data, test_data, target_var, predictor_vars, k = 1, l = 0) {
  # Convert to data frames if they are data tables
  if (is.data.table(train_data)) {
    train_data <- as.data.frame(train_data)
  }
  if (is.data.table(test_data)) {
    test_data <- as.data.frame(test_data)
  }
  
  # Validate inputs
  if (!target_var %in% names(train_data)) {
    stop("Target variable '", target_var, "' not found in training data")
  }
  
  if (!target_var %in% names(test_data)) {
    stop("Target variable '", target_var, "' not found in test data")
  }
  
  if (!all(predictor_vars %in% names(train_data))) {
    missing_vars <- setdiff(predictor_vars, names(train_data))
    stop("Predictor variable(s) not found in training data: ", paste(missing_vars, collapse = ", "))
  }
  
  if (!all(predictor_vars %in% names(test_data))) {
    missing_vars <- setdiff(predictor_vars, names(test_data))
    stop("Predictor variable(s) not found in test data: ", paste(missing_vars, collapse = ", "))
  }
  
  if (k <= 0 || k > nrow(train_data)) {
    stop("k must be a positive integer less than or equal to the number of training observations")
  }
  
  # Extract features and target from training data
  train_features <- train_data[, predictor_vars, drop = FALSE]
  train_labels <- train_data[[target_var]]
  
  # Extract features from test data (target is only for evaluation)
  test_features <- test_data[, predictor_vars, drop = FALSE]
  test_actual <- test_data[[target_var]]
  
  # Ensure target is a factor for classification
  if (!is.factor(train_labels)) {
    warning("Converting training target variable to factor")
    train_labels <- as.factor(train_labels)
  }
  
  if (!is.factor(test_actual)) {
    test_actual <- as.factor(test_actual)
  }
  
  # Apply KNN with fixed prob = TRUE and use.all = TRUE
  predictions <- knn(
    train = train_features,
    test = test_features,
    cl = train_labels,
    k = k,
    l = l,
    prob = TRUE,      # Fixed to TRUE
    use.all = TRUE    # Fixed to TRUE
  )
  
  # Extract probabilities
  probabilities <- attr(predictions, "prob")
  
  # Calculate accuracy
  accuracy <- mean(predictions == test_actual)
  
  # Create confusion matrix
  conf_matrix <- table(Predicted = predictions, Actual = test_actual)
  
  # Return comprehensive results
  result <- list(
    predictions = predictions,
    probabilities = probabilities,
    test_actual = test_actual,
    accuracy = accuracy,
    confusion_matrix = conf_matrix,
    target_variable = target_var,
    predictor_variables = predictor_vars,
    k = k,
    l = l,
    classes = levels(train_labels),
    n_train = nrow(train_data),
    n_test = nrow(test_data)
  )
  
  class(result) <- "knn_result"
  return(result)
}
,
KNN_reg = function(train_data, test_data, target_var, predictor_vars, k = 5) {
  # Convert to data frames if they are data tables
  if (is.data.table(train_data)) {
    train_data <- as.data.frame(train_data)
  }
  if (is.data.table(test_data)) {
    test_data <- as.data.frame(test_data)
  }
  
  # Validate inputs
  if (!target_var %in% names(train_data)) {
    stop("Target variable '", target_var, "' not found in training data")
  }
  
  if (!target_var %in% names(test_data)) {
    stop("Target variable '", target_var, "' not found in test data")
  }
  
  if (!all(predictor_vars %in% names(train_data))) {
    missing_vars <- setdiff(predictor_vars, names(train_data))
    stop("Predictor variable(s) not found in training data: ", paste(missing_vars, collapse = ", "))
  }
  
  if (!all(predictor_vars %in% names(test_data))) {
    missing_vars <- setdiff(predictor_vars, names(test_data))
    stop("Predictor variable(s) not found in test data: ", paste(missing_vars, collapse = ", "))
  }
  
  if (k <= 0 || k > nrow(train_data)) {
    stop("k must be a positive integer less than or equal to the number of training observations")
  }
  
  # Extract features and target from training data
  train_features <- train_data[, predictor_vars, drop = FALSE]
  train_target <- train_data[[target_var]]
  
  # Extract features from test data (target is only for evaluation)
  test_features <- test_data[, predictor_vars, drop = FALSE]
  test_actual <- test_data[[target_var]]
  
  # Ensure target is numeric for regression
  if (!is.numeric(train_target)) {
    warning("Converting training target variable to numeric")
    train_target <- as.numeric(train_target)
  }
  
  if (!is.numeric(test_actual)) {
    test_actual <- as.numeric(test_actual)
  }
  
  # Apply KNN Regression
  knn_result <- knn.reg(
    train = train_features,
    test = test_features,
    y = train_target,
    k = k
  )
  
  predictions <- knn_result$pred
  
  # Calculate regression metrics
  residuals <- test_actual - predictions
  mse <- mean(residuals^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(residuals))
  r_squared <- 1 - (sum(residuals^2) / sum((test_actual - mean(test_actual))^2))
  
  # Return comprehensive results
  result <- list(
    predictions = predictions,
    test_actual = test_actual,
    residuals = residuals,
    mse = mse,
    rmse = rmse,
    mae = mae,
    r_squared = r_squared,
    target_variable = target_var,
    predictor_variables = predictor_vars,
    k = k,
    n_train = nrow(train_data),
    n_test = nrow(test_data)
  )
  
  class(result) <- "knn_regression_result"
  return(result)
}
,
generate_correlation_plots_scatter = function(data, variables, 
                                              correlation_method = "pearson",
                                              pvalue = 0.05, 
                                              threshold = 0.5, 
                                              output_folder = "./correlation_plots",
                                              width = 8, height = 6, dpi = 300,
                                              smooth_method = "lm") {   
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("The 'data' argument must be a data frame")
  }
  
  # Validate correlation method
  valid_methods <- c("pearson", "spearman", "kendall")
  if (!correlation_method %in% valid_methods) {
    stop("Invalid correlation_method. Choose from: ", 
         paste(valid_methods, collapse = ", "))
  }
  
  if (!all(variables %in% names(data))) {
    missing_vars <- variables[!variables %in% names(data)]
    stop("The following variables are not in the data frame: ", 
         paste(missing_vars, collapse = ", "))
  }
  
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    message("Created output folder: ", output_folder)
  }
  
  # Subset data to selected variables
  data_subset <- data[, variables, drop = FALSE]
  
  # Remove non-numeric variables with warning
  numeric_vars <- sapply(data_subset, is.numeric)
  if (!all(numeric_vars)) {
    non_numeric <- variables[!numeric_vars]
    warning("The following variables are not numeric and will be excluded: ",
            paste(non_numeric, collapse = ", "))
    variables <- variables[numeric_vars]
    data_subset <- data_subset[, numeric_vars, drop = FALSE]
  }
  
  if (length(variables) < 2) {
    stop("At least 2 numeric variables are required for correlation analysis")
  }
  
  # Calculate correlation matrix with specified method
  cor_matrix <- cor(data_subset, 
                    method = correlation_method, 
                    use = "pairwise.complete.obs")
  
  # Function to calculate correlation p-values with specified method
  cor_pvalue <- function(x, y, method = correlation_method) {
    cor_test <- cor.test(x, y, method = method, use = "pairwise.complete.obs")
    return(cor_test$p.value)
  }
  
  # Create p-value matrix
  pvalue_matrix <- matrix(NA, nrow = length(variables), ncol = length(variables))
  rownames(pvalue_matrix) <- colnames(pvalue_matrix) <- variables
  
  for (i in 1:length(variables)) {
    for (j in 1:length(variables)) {
      if (i != j) {
        pvalue_matrix[i, j] <- cor_pvalue(data_subset[[variables[i]]], 
                                          data_subset[[variables[j]]],
                                          method = correlation_method)
      }
    }
  }
  
  # Create results data frame
  results <- data.frame()
  
  # Generate plots for significant correlations
  plot_count <- 0
  
  # Get correlation symbol for plot labels
  cor_symbol <- switch(correlation_method,
                       "pearson" = "r",
                       "spearman" = "ρ", 
                       "kendall" = "τ")
  
  for (i in 1:(length(variables)-1)) {
    for (j in (i+1):length(variables)) {
      var1 <- variables[i]
      var2 <- variables[j]
      cor_value <- cor_matrix[i, j]
      p_val <- pvalue_matrix[i, j]
      
      # Check if correlation meets criteria
      if (!is.na(p_val) && p_val < pvalue && abs(cor_value) >= threshold) {
        
        # Create plot - ALWAYS use linear regression for consistency
        p <- ggplot2::ggplot(data_subset, ggplot2::aes(x = .data[[var1]], y = .data[[var2]])) +
          ggplot2::geom_point(alpha = 0.6, size = 2) +
          ggplot2::geom_smooth(method = smooth_method,  # Use the specified smooth method
                               se = TRUE, color = "red", linetype = "solid") +
          ggplot2::labs(
            title = paste("Correlation:", var1, "vs", var2),
            subtitle = paste0(cor_symbol, " = ", round(cor_value, 3), 
                              ", p = ", format.pval(p_val, digits = 3),
                              " (", correlation_method, ")"),
            x = var1,
            y = var2
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 14),
            plot.subtitle = ggplot2::element_text(size = 11, color = "darkred")
          )
        
        # Create filename with method indication
        filename <- file.path(output_folder, 
                              paste0("correlation_", correlation_method, "_", 
                                     var1, "_", var2, ".png"))
        
        # Save plot
        ggplot2::ggsave(
          filename = filename,
          plot = p,
          width = width,
          height = height,
          dpi = dpi,
          device = "png"
        )
        
        plot_count <- plot_count + 1
        
        # Add to results
        results <- rbind(results, data.frame(
          Variable1 = var1,
          Variable2 = var2,
          Correlation = cor_value,
          P_value = p_val,
          Method = correlation_method,
          Plot_file = basename(filename),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Print summary
  message("Correlation analysis completed using '", correlation_method, "' method")
  message("Smoothing method used in plots: '", smooth_method, "'")
  if (plot_count == 0) {
    message("No correlations met the specified criteria (p < ", pvalue, 
            " and |", cor_symbol, "| >= ", threshold, ")")
  } else {
    message("Generated ", plot_count, " scatter plots in '", output_folder, "'")
    message("Correlation threshold: |", cor_symbol, "| >= ", threshold)
    message("Significance level: p < ", pvalue)
  }
  
  # Return results invisibly
  invisible(results)
}
,

random_forest_model = function(train_data, test_data, target_var, input_vars, 
                               ntree = 500, mtry = NULL, regression = TRUE) {
  
  # Create formula
  formula <- as.formula(paste(target_var, "~", paste(input_vars, collapse = " + ")))
  
  # Set mtry if not provided
  if (is.null(mtry)) {
    mtry <- if (regression) max(floor(length(input_vars)/3), 1) else floor(sqrt(length(input_vars)))
  }
  
  # Train Random Forest model
  if (regression) {
    model <- randomForest(
      formula,
      data = train_data,
      ntree = ntree,
      mtry = mtry,
      importance = TRUE,
      na.action = na.omit
    )
  } else {
    model <- randomForest(
      formula,
      data = train_data,
      ntree = ntree,
      mtry = mtry,
      importance = TRUE,
      na.action = na.omit
    )
  }
  
  # Make predictions on test set
  predictions <- predict(model, newdata = test_data)
  
  # Get actual values
  actual <- test_data[[target_var]]
  
  # Calculate comprehensive performance metrics
  if (regression) {
    # Remove NA values for calculation
    valid_idx <- !is.na(actual) & !is.na(predictions)
    actual_clean <- actual[valid_idx]
    predicted_clean <- predictions[valid_idx]
    
    # Regression metrics
    residuals <- actual_clean - predicted_clean
    mae <- mean(abs(residuals))
    rmse <- sqrt(mean(residuals^2))
    r_squared <- cor(actual_clean, predicted_clean)^2
    mape <- mean(abs(residuals / actual_clean)) * 100  # Mean Absolute Percentage Error
    mean_error <- mean(residuals)  # Bias
    median_ae <- median(abs(residuals))
    
    metrics <- list(
      MAE = mae,
      RMSE = rmse,
      R_squared = r_squared,
      MAPE = mape,
      Mean_Error = mean_error,
      Median_AE = median_ae,
      Bias_Direction = ifelse(mean_error > 0, "Under-predicting", "Over-predicting"),
      Predictions = data.frame(
        Actual = actual, 
        Predicted = predictions,
        Residuals = actual - predictions
      ),
      Correlation_Test = cor.test(actual_clean, predicted_clean),
      Summary = data.frame(
        Metric = c("RMSE", "MAE", "R²", "MAPE", "Mean_Error", "Median_AE"),
        Value = c(rmse, mae, r_squared, mape, mean_error, median_ae)
      )
    )
    
  } else {
    # Classification metrics
    confusion_matrix <- confusionMatrix(predictions, actual)
    
    # Calculate AUC for binary classification
    auc_value <- NA
    if (length(unique(actual)) == 2) {
      tryCatch({
        roc_obj <- roc(actual, as.numeric(predictions))
        auc_value <- auc(roc_obj)
      }, error = function(e) {
        auc_value <- NA
      })
    }
    
    metrics <- list(
      Confusion_Matrix = confusion_matrix,
      Accuracy = confusion_matrix$overall["Accuracy"],
      Kappa = confusion_matrix$overall["Kappa"],
      Sensitivity = confusion_matrix$byClass["Sensitivity"],
      Specificity = confusion_matrix$byClass["Specificity"],
      F1_Score = confusion_matrix$byClass["F1"],
      AUC = auc_value,
      Predictions = data.frame(Actual = actual, Predicted = predictions),
      Class_Proportions = list(
        Actual = prop.table(table(actual)),
        Predicted = prop.table(table(predictions))
      )
    )
  }
  
  # Get variable importance
  var_importance <- importance(model)
  
  # Return results
  result <- list(
    model = model,
    predictions = predictions,
    metrics = metrics,
    variable_importance = var_importance,
    parameters = list(
      ntree = ntree,
      mtry = mtry,
      regression = regression,
      n_predictors = length(input_vars)
    ),
    data_info = list(
      train_size = nrow(train_data),
      test_size = nrow(test_data),
      target_variable = target_var
    )
  )
  
  return(result)
}
,
xgboost_model = function(train_data, test_data, target_var, input_vars, 
                         nrounds = 100, params = NULL, regression = TRUE, 
                         cv_folds = 5) {
  
  # Prepare data for XGBoost
  prepare_xgb_data <- function(data, target_var, input_vars) {
    # Convert to matrix (XGBoost requires numeric matrix)
    x_data <- as.matrix(data[, input_vars])
    
    if (regression) {
      y_data <- data[[target_var]]
    } else {
      # Convert factor to numeric (0-indexed)
      y_data <- as.numeric(data[[target_var]]) - 1
    }
    
    return(list(x = x_data, y = y_data))
  }
  
  # Prepare train and test data
  train_prep <- prepare_xgb_data(train_data, target_var, input_vars)
  test_prep <- prepare_xgb_data(test_data, target_var, input_vars)
  
  # Set default parameters if not provided
  if (is.null(params)) {
    if (regression) {
      params <- list(
        objective = "reg:squarederror",
        max_depth = 6,
        eta = 0.1,
        gamma = 0,
        colsample_bytree = 0.8,
        min_child_weight = 1,
        subsample = 0.8
      )
    } else {
      params <- list(
        objective = "binary:logistic",
        max_depth = 6,
        eta = 0.1,
        gamma = 0,
        colsample_bytree = 0.8,
        min_child_weight = 1,
        subsample = 0.8
      )
    }
  }
  
  # Create DMatrix objects (XGBoost native format)
  dtrain <- xgb.DMatrix(data = train_prep$x, label = train_prep$y)
  dtest <- xgb.DMatrix(data = test_prep$x, label = test_prep$y)
  
  # Use cross-validation to find optimal nrounds
  cv_result <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    nfold = cv_folds,
    metrics = if(regression) "rmse" else "error",
    early_stopping_rounds = 10,
    verbose = FALSE
  )
  
  # Get optimal number of rounds
  best_nrounds <- if (!is.null(cv_result$best_iteration)) {
    cv_result$best_iteration
  } else {
    nrounds
  }
  
  # Train final model
  model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_nrounds,
    watchlist = list(train = dtrain, test = dtest),
    verbose = 0
  )
  
  # Make predictions
  if (regression) {
    predictions <- predict(model, dtest)
  } else {
    pred_probs <- predict(model, dtest)
    predictions <- ifelse(pred_probs > 0.5, 1, 0)
    # Convert back to original factor levels
    predictions <- factor(predictions, levels = c(0, 1), 
                          labels = levels(train_data[[target_var]]))
  }
  
  # Get actual values
  actual <- test_data[[target_var]]
  
  # Calculate comprehensive performance metrics
  if (regression) {
    # Remove NA values for calculation
    valid_idx <- !is.na(actual) & !is.na(predictions)
    actual_clean <- actual[valid_idx]
    predicted_clean <- predictions[valid_idx]
    
    # Regression metrics
    residuals <- actual_clean - predicted_clean
    mae <- mean(abs(residuals))
    rmse <- sqrt(mean(residuals^2))
    r_squared <- cor(actual_clean, predicted_clean)^2
    mape <- mean(abs(residuals / actual_clean)) * 100
    mean_error <- mean(residuals)
    median_ae <- median(abs(residuals))
    
    metrics <- list(
      MAE = mae,
      RMSE = rmse,
      R_squared = r_squared,
      MAPE = mape,
      Mean_Error = mean_error,
      Median_AE = median_ae,
      Bias_Direction = ifelse(mean_error > 0, "Under-predicting", "Over-predicting"),
      Predictions = data.frame(
        Actual = actual, 
        Predicted = predictions,
        Residuals = actual - predictions
      ),
      CV_Results = cv_result,
      Best_Iteration = best_nrounds
    )
    
  } else {
    # Classification metrics
    confusion_matrix <- confusionMatrix(predictions, actual)
    
    # Calculate AUC for binary classification
    auc_value <- NA
    if (length(unique(actual)) == 2) {
      tryCatch({
        roc_obj <- roc(actual, as.numeric(predictions))
        auc_value <- auc(roc_obj)
      }, error = function(e) {
        auc_value <- NA
      })
    }
    
    metrics <- list(
      Confusion_Matrix = confusion_matrix,
      Accuracy = confusion_matrix$overall["Accuracy"],
      Kappa = confusion_matrix$overall["Kappa"],
      Sensitivity = confusion_matrix$byClass["Sensitivity"],
      Specificity = confusion_matrix$byClass["Specificity"],
      F1_Score = confusion_matrix$byClass["F1"],
      AUC = auc_value,
      Predictions = data.frame(Actual = actual, Predicted = predictions),
      CV_Results = cv_result,
      Best_Iteration = best_nrounds
    )
  }
  
  # Get variable importance
  importance_matrix <- xgb.importance(
    feature_names = input_vars, 
    model = model
  )
  
  # Return results
  result <- list(
    model = model,
    predictions = predictions,
    metrics = metrics,
    variable_importance = importance_matrix,
    parameters = list(
      nrounds = best_nrounds,
      params = params,
      regression = regression,
      n_predictors = length(input_vars)
    ),
    data_info = list(
      train_size = nrow(train_data),
      test_size = nrow(test_data),
      target_variable = target_var
    )
  )
  
  return(result)
}
,

dnn_model = function(train_data, test_data, target_var, input_vars, 
                     layers = c(64, 32), activation = "relu", dropout = 0.2,
                     learning_rate = 0.001, epochs = 100, batch_size = 32,
                     validation_split = 0.2, regression = TRUE, verbose = 0) {
  get_variable_importance <- function(model, feature_names) {
    tryCatch({
      # Get weights from first layer
      weights <- get_weights(model)[[1]]
      # Calculate importance as mean absolute weight
      importance <- apply(abs(weights), 1, mean)
      names(importance) <- feature_names
      return(sort(importance, decreasing = TRUE))
    }, error = function(e) {
      return(NULL)
    })
  }
  # Prepare data for neural network
  prepare_nn_data <- function(data, target_var, input_vars, regression) {
    # Scale features (important for neural networks)
    x_data <- as.matrix(data[, input_vars])
    
    # Scale to 0-1 range
    x_scaled <- scale(x_data)
    
    if (regression) {
      y_data <- data[[target_var]]
      # Scale target for regression
      y_mean <- mean(y_data, na.rm = TRUE)
      y_sd <- sd(y_data, na.rm = TRUE)
      y_scaled <- (y_data - y_mean) / y_sd
      return(list(x = x_scaled, y = y_scaled, y_stats = c(y_mean, y_sd)))
    } else {
      # For classification, convert to appropriate format
      y_data <- data[[target_var]]
      if (length(unique(y_data)) == 2) {
        # Binary classification - convert to 0/1
        y_encoded <- as.numeric(y_data) - 1
        return(list(x = x_scaled, y = y_encoded, y_levels = levels(y_data)))
      } else {
        # Multi-class classification - one-hot encoding
        y_encoded <- to_categorical(as.numeric(y_data) - 1)
        return(list(x = x_scaled, y = y_encoded, y_levels = levels(y_data)))
      }
    }
  }
  
  # Prepare train and test data
  train_prep <- prepare_nn_data(train_data, target_var, input_vars, regression)
  test_prep <- prepare_nn_data(test_data, target_var, input_vars, regression)
  
  # Get input shape
  input_shape <- ncol(train_prep$x)
  
  # Build model
  model <- keras_model_sequential()
  
  # Input layer
  model %>% 
    layer_dense(units = layers[1], activation = activation, input_shape = input_shape)
  
  # Add dropout if specified
  if (dropout > 0) {
    model %>% layer_dropout(rate = dropout)
  }
  
  # Additional hidden layers
  if (length(layers) > 1) {
    for (i in 2:length(layers)) {
      model %>% layer_dense(units = layers[i], activation = activation)
      if (dropout > 0) {
        model %>% layer_dropout(rate = dropout)
      }
    }
  }
  
  # Output layer
  if (regression) {
    model %>% layer_dense(units = 1, activation = "linear")
  } else {
    if (length(unique(train_data[[target_var]])) == 2) {
      # Binary classification
      model %>% layer_dense(units = 1, activation = "sigmoid")
    } else {
      # Multi-class classification
      model %>% layer_dense(units = length(unique(train_data[[target_var]])), 
                            activation = "softmax")
    }
  }
  
  # Compile model
  if (regression) {
    model %>% compile(
      loss = "mse",
      optimizer = optimizer_adam(learning_rate = learning_rate),
      metrics = c("mae")
    )
  } else {
    model %>% compile(
      loss = if (length(unique(train_data[[target_var]])) == 2) "binary_crossentropy" else "categorical_crossentropy",
      optimizer = optimizer_adam(learning_rate = learning_rate),
      metrics = c("accuracy")
    )
  }
  
  # Set up early stopping
  early_stop <- callback_early_stopping(
    monitor = "val_loss",
    patience = 10,
    restore_best_weights = TRUE
  )
  
  # Train model
  history <- model %>% fit(
    x = train_prep$x,
    y = train_prep$y,
    epochs = epochs,
    batch_size = batch_size,
    validation_split = validation_split,
    callbacks = list(early_stop),
    verbose = verbose
  )
  
  # Make predictions
  raw_predictions <- predict(model, test_prep$x)
  
  # Process predictions based on problem type
  if (regression) {
    # Rescale predictions back to original scale
    predictions <- raw_predictions * test_prep$y_stats[2] + test_prep$y_stats[1]
  } else {
    if (length(unique(train_data[[target_var]])) == 2) {
      # Binary classification
      predictions <- ifelse(raw_predictions > 0.5, 1, 0)
      predictions <- factor(predictions, levels = c(0, 1), 
                            labels = test_prep$y_levels)
    } else {
      # Multi-class classification
      predictions <- factor(max.col(raw_predictions), 
                            levels = 1:length(test_prep$y_levels),
                            labels = test_prep$y_levels)
    }
  }
  
  # Get actual values
  actual <- test_data[[target_var]]
  
  # Calculate comprehensive performance metrics
  if (regression) {
    # Remove NA values for calculation
    valid_idx <- !is.na(actual) & !is.na(predictions)
    actual_clean <- actual[valid_idx]
    predicted_clean <- predictions[valid_idx]
    
    # Regression metrics
    residuals <- actual_clean - predicted_clean
    mae <- mean(abs(residuals))
    rmse <- sqrt(mean(residuals^2))
    r_squared <- cor(actual_clean, predicted_clean)^2
    mape <- mean(abs(residuals / actual_clean)) * 100
    mean_error <- mean(residuals)
    median_ae <- median(abs(residuals))
    
    metrics <- list(
      MAE = mae,
      RMSE = rmse,
      R_squared = r_squared,
      MAPE = mape,
      Mean_Error = mean_error,
      Median_AE = median_ae,
      Bias_Direction = ifelse(mean_error > 0, "Under-predicting", "Over-predicting"),
      Predictions = data.frame(
        Actual = actual, 
        Predicted = predictions,
        Residuals = actual - predictions
      ),
      Training_History = history
    )
    
  } else {
    # Classification metrics
    confusion_matrix <- confusionMatrix(predictions, actual)
    
    # Calculate AUC for binary classification
    auc_value <- NA
    if (length(unique(actual)) == 2) {
      tryCatch({
        roc_obj <- roc(actual, as.numeric(predictions))
        auc_value <- auc(roc_obj)
      }, error = function(e) {
        auc_value <- NA
      })
    }
    
    metrics <- list(
      Confusion_Matrix = confusion_matrix,
      Accuracy = confusion_matrix$overall["Accuracy"],
      Kappa = confusion_matrix$overall["Kappa"],
      Sensitivity = confusion_matrix$byClass["Sensitivity"],
      Specificity = confusion_matrix$byClass["Specificity"],
      F1_Score = confusion_matrix$byClass["F1"],
      AUC = auc_value,
      Predictions = data.frame(Actual = actual, Predicted = predictions),
      Training_History = history
    )
  }
  
  # Get variable importance (approximate using connection weights)
  variable_importance <- get_variable_importance(model, input_vars)
  
  # Return results
  result <- list(
    model = model,
    predictions = predictions,
    metrics = metrics,
    variable_importance = variable_importance,
    parameters = list(
      layers = layers,
      activation = activation,
      dropout = dropout,
      learning_rate = learning_rate,
      epochs = epochs,
      batch_size = batch_size,
      validation_split = validation_split,
      regression = regression,
      n_predictors = length(input_vars)
    ),
    data_info = list(
      train_size = nrow(train_data),
      test_size = nrow(test_data),
      target_variable = target_var
    ),
    scaling_info = list(
      input_means = attr(train_prep$x, "scaled:center"),
      input_sds = attr(train_prep$x, "scaled:scale")
    )
  )
  
  return(result)
}
,
knn_model = function(train_data, test_data, target_var, input_vars, 
                     k = 5, scale_data = TRUE, regression = TRUE, 
                     tune_k = FALSE, tune_length = 10) {
  
  # Prepare data for KNN
  prepare_knn_data <- function(data, target_var, input_vars, scale_data) {
    x_data <- as.matrix(data[, input_vars])
    y_data <- data[[target_var]]
    
    if (scale_data) {
      # Scale features (important for KNN)
      x_scaled <- scale(x_data)
      return(list(x = x_scaled, y = y_data, 
                  center = attr(x_scaled, "scaled:center"),
                  scale = attr(x_scaled, "scaled:scale")))
    } else {
      return(list(x = x_data, y = y_data))
    }
  }
  
  # Prepare train and test data
  train_prep <- prepare_knn_data(train_data, target_var, input_vars, scale_data)
  test_prep <- prepare_knn_data(test_data, target_var, input_vars, scale_data)
  
  # If scaling was applied to train, apply same scaling to test
  if (scale_data && !is.null(train_prep$center)) {
    test_prep$x <- scale(test_prep$x, 
                         center = train_prep$center, 
                         scale = train_prep$scale)
  }
  
  # Tune k if requested
  if (tune_k) {
    if (regression) {
      # For regression, use custom tuning
      k_values <- seq(1, min(tune_length, nrow(train_data) %/% 2), by = 2)
      cv_results <- data.frame()
      
      for (k_val in k_values) {
        # Simple cross-validation
        set.seed(123)
        cv_folds <- createFolds(train_prep$y, k = 5)
        cv_errors <- c()
        
        for (fold in cv_folds) {
          train_fold_x <- train_prep$x[-fold, ]
          train_fold_y <- train_prep$y[-fold]
          test_fold_x <- train_prep$x[fold, ]
          test_fold_y <- train_prep$y[fold]
          
          pred <- FNN::knn.reg(train = train_fold_x, 
                               test = test_fold_x, 
                               y = train_fold_y, 
                               k = k_val)$pred
          cv_errors <- c(cv_errors, sqrt(mean((pred - test_fold_y)^2)))
        }
        
        cv_results <- rbind(cv_results, 
                            data.frame(k = k_val, RMSE = mean(cv_errors)))
      }
      
      best_k <- cv_results$k[which.min(cv_results$RMSE)]
      
    } else {
      # For classification, use caret's built-in tuning
      ctrl <- trainControl(method = "cv", number = 5)
      tune_grid <- expand.grid(k = seq(1, min(tune_length, nrow(train_data) %/% 2), by = 2))
      
      knn_tune <- train(x = train_prep$x,
                        y = train_prep$y,
                        method = "knn",
                        trControl = ctrl,
                        tuneGrid = tune_grid,
                        preProcess = if(scale_data) c("center", "scale") else NULL)
      
      best_k <- knn_tune$bestTune$k
      cv_results <- knn_tune$results
    }
  } else {
    best_k <- k
    cv_results <- NULL
  }
  
  # Train final KNN model and make predictions
  if (regression) {
    knn_pred <- FNN::knn.reg(train = train_prep$x, 
                             test = test_prep$x, 
                             y = train_prep$y, 
                             k = best_k)
    predictions <- knn_pred$pred
  } else {
    knn_pred <- FNN::knn(train = train_prep$x, 
                         test = test_prep$x, 
                         cl = train_prep$y, 
                         k = best_k,
                         prob = TRUE)
    predictions <- knn_pred
  }
  
  # Get actual values
  actual <- test_data[[target_var]]
  
  # Calculate comprehensive performance metrics
  if (regression) {
    # Remove NA values for calculation
    valid_idx <- !is.na(actual) & !is.na(predictions)
    actual_clean <- actual[valid_idx]
    predicted_clean <- predictions[valid_idx]
    
    # Regression metrics
    residuals <- actual_clean - predicted_clean
    mae <- mean(abs(residuals))
    rmse <- sqrt(mean(residuals^2))
    r_squared <- cor(actual_clean, predicted_clean)^2
    mape <- mean(abs(residuals / actual_clean)) * 100
    mean_error <- mean(residuals)
    median_ae <- median(abs(residuals))
    
    metrics <- list(
      MAE = mae,
      RMSE = rmse,
      R_squared = r_squared,
      MAPE = mape,
      Mean_Error = mean_error,
      Median_AE = median_ae,
      Bias_Direction = ifelse(mean_error > 0, "Under-predicting", "Over-predicting"),
      Predictions = data.frame(
        Actual = actual, 
        Predicted = predictions,
        Residuals = actual - predictions
      )
    )
    
  } else {
    # Classification metrics
    confusion_matrix <- confusionMatrix(predictions, actual)
    
    # Calculate AUC for binary classification
    auc_value <- NA
    if (length(unique(actual)) == 2) {
      tryCatch({
        prob_attr <- attr(predictions, "prob")
        # Adjust probabilities for the winning class
        adjusted_probs <- ifelse(predictions == levels(actual)[1], 
                                 1 - prob_attr, prob_attr)
        roc_obj <- roc(actual, adjusted_probs)
        auc_value <- auc(roc_obj)
      }, error = function(e) {
        auc_value <- NA
      })
    }
    
    metrics <- list(
      Confusion_Matrix = confusion_matrix,
      Accuracy = confusion_matrix$overall["Accuracy"],
      Kappa = confusion_matrix$overall["Kappa"],
      Sensitivity = confusion_matrix$byClass["Sensitivity"],
      Specificity = confusion_matrix$byClass["Specificity"],
      F1_Score = confusion_matrix$byClass["F1"],
      AUC = auc_value,
      Predictions = data.frame(Actual = actual, Predicted = predictions)
    )
  }
  
  # Store scaling information
  scaling_info <- if (scale_data) {
    list(center = train_prep$center, scale = train_prep$scale)
  } else {
    NULL
  }
  
  # Return results
  result <- list(
    model = list(
      train_x = train_prep$x,
      train_y = train_prep$y,
      k = best_k,
      scaling_info = scaling_info,
      input_vars = input_vars,
      regression = regression
    ),
    predictions = predictions,
    metrics = metrics,
    tuning_results = cv_results,
    parameters = list(
      k = best_k,
      scale_data = scale_data,
      regression = regression,
      tuned = tune_k,
      n_predictors = length(input_vars)
    ),
    data_info = list(
      train_size = nrow(train_data),
      test_size = nrow(test_data),
      target_variable = target_var
    )
  )
  
  return(result)
}
,
extra_trees_ranger = function(train_data, test_data, target_var, input_vars, 
                              ntree = 500, mtry = NULL, min_node_size = 5,
                              splitrule = "extratrees", regression = TRUE, 
                              num_random_splits = 1, tune_mtry = FALSE) {
  
  # Create formula
  formula <- as.formula(paste(target_var, "~", paste(input_vars, collapse = " + ")))
  
  # Set mtry if not provided
  if (is.null(mtry)) {
    mtry <- if (regression) max(floor(length(input_vars)/3), 1) else floor(sqrt(length(input_vars)))
  }
  
  # Tune mtry if requested
  if (tune_mtry) {
    mtry_values <- unique(round(seq(2, length(input_vars), length.out = 5)))
    cv_results <- data.frame()
    
    for (mtry_val in mtry_values) {
      set.seed(123)
      cv_folds <- createFolds(train_data[[target_var]], k = 5)
      cv_errors <- c()
      
      for (fold in cv_folds) {
        train_fold <- train_data[-fold, ]
        test_fold <- train_data[fold, ]
        
        model_cv <- ranger(
          formula = formula,
          data = train_fold,
          num.trees = 100,
          mtry = mtry_val,
          min.node.size = min_node_size,
          splitrule = splitrule,
          num.random.splits = num_random_splits,
          importance = "impurity"
        )
        
        pred <- predict(model_cv, test_fold)$predictions
        if (regression) {
          cv_errors <- c(cv_errors, sqrt(mean((pred - test_fold[[target_var]])^2)))
        } else {
          cv_errors <- c(cv_errors, 1 - mean(pred == test_fold[[target_var]]))
        }
      }
      
      cv_results <- rbind(cv_results, 
                          data.frame(mtry = mtry_val, Error = mean(cv_errors)))
    }
    
    best_mtry <- cv_results$mtry[which.min(cv_results$Error)]
  } else {
    best_mtry <- mtry
    cv_results <- NULL
  }
  
  # Train final model
  model <- ranger(
    formula = formula,
    data = train_data,
    num.trees = ntree,
    mtry = best_mtry,
    min.node.size = min_node_size,
    splitrule = splitrule,
    num.random.splits = num_random_splits,
    importance = "impurity",
    seed = 123
  )
  
  # Make predictions
  predictions <- predict(model, test_data)$predictions
  
  # Get actual values
  actual <- test_data[[target_var]]
  
  # Calculate comprehensive performance metrics
  if (regression) {
    valid_idx <- !is.na(actual) & !is.na(predictions)
    actual_clean <- actual[valid_idx]
    predicted_clean <- predictions[valid_idx]
    
    residuals <- actual_clean - predicted_clean
    mae <- mean(abs(residuals))
    rmse <- sqrt(mean(residuals^2))
    r_squared <- cor(actual_clean, predicted_clean)^2
    mape <- mean(abs(residuals / actual_clean)) * 100
    mean_error <- mean(residuals)
    median_ae <- median(abs(residuals))
    
    metrics <- list(
      MAE = mae,
      RMSE = rmse,
      R_squared = r_squared,
      MAPE = mape,
      Mean_Error = mean_error,
      Median_AE = median_ae,
      Bias_Direction = ifelse(mean_error > 0, "Under-predicting", "Over-predicting"),
      Predictions = data.frame(
        Actual = actual, 
        Predicted = predictions,
        Residuals = actual - predictions
      )
    )
    
  } else {
    confusion_matrix <- confusionMatrix(predictions, actual)
    
    auc_value <- NA
    if (length(unique(actual)) == 2) {
      tryCatch({
        pred_probs <- predict(model, test_data, type = "response")$predictions
        if (is.matrix(pred_probs)) {
          positive_probs <- pred_probs[, 2]
        } else {
          positive_probs <- pred_probs
        }
        roc_obj <- roc(actual, positive_probs)
        auc_value <- auc(roc_obj)
      }, error = function(e) {
        auc_value <- NA
      })
    }
    
    metrics <- list(
      Confusion_Matrix = confusion_matrix,
      Accuracy = confusion_matrix$overall["Accuracy"],
      Kappa = confusion_matrix$overall["Kappa"],
      Sensitivity = confusion_matrix$byClass["Sensitivity"],
      Specificity = confusion_matrix$byClass["Specificity"],
      F1_Score = confusion_matrix$byClass["F1"],
      AUC = auc_value,
      Predictions = data.frame(Actual = actual, Predicted = predictions)
    )
  }
  
  # Get variable importance
  var_importance <- model$variable.importance
  
  # Return results
  result <- list(
    model = model,
    predictions = predictions,
    metrics = metrics,
    variable_importance = var_importance,
    tuning_results = cv_results,
    parameters = list(
      ntree = ntree,
      mtry = best_mtry,
      min_node_size = min_node_size,
      splitrule = splitrule,
      num_random_splits = num_random_splits,
      regression = regression,
      tuned = tune_mtry,
      n_predictors = length(input_vars)
    ),
    data_info = list(
      train_size = nrow(train_data),
      test_size = nrow(test_data),
      target_variable = target_var
    )
  )
  
  return(result)
}
,
svm_model = function(train_data, test_data, target_var, input_vars, 
                      kernel = "radial", cost = 1, gamma = "auto",
                      scale = TRUE, regression = TRUE) {
  
  library(e1071)
  library(caret)
  
  # Create formula
  formula <- as.formula(paste(target_var, "~", paste(input_vars, collapse = " + ")))
  
  # Set automatic gamma if not specified
  if (gamma == "auto") {
    n_features <- length(input_vars)
    gamma <- 1 / n_features  # Common default for automatic gamma
  }
  
  # Train SVM model
  if (regression) {
    model <- svm(
      formula,
      data = train_data,
      type = "eps-regression",
      kernel = kernel,
      cost = cost,
      gamma = gamma,
      scale = scale
    )
  } else {
    model <- svm(
      formula,
      data = train_data,
      type = "C-classification",
      kernel = kernel,
      cost = cost,
      gamma = gamma,
      probability = TRUE,
      scale = scale
    )
  }
  
  # Make predictions
  predictions <- predict(model, test_data)
  
  # Get actual values
  actual <- test_data[[target_var]]
  
  # Calculate comprehensive performance metrics
  if (regression) {
    valid_idx <- !is.na(actual) & !is.na(predictions)
    actual_clean <- actual[valid_idx]
    predicted_clean <- predictions[valid_idx]
    
    residuals <- actual_clean - predicted_clean
    mae <- mean(abs(residuals))
    rmse <- sqrt(mean(residuals^2))
    r_squared <- cor(actual_clean, predicted_clean)^2
    mape <- mean(abs(residuals / actual_clean)) * 100
    mean_error <- mean(residuals)
    
    metrics <- list(
      MAE = mae,
      RMSE = rmse,
      R_squared = r_squared,
      MAPE = mape,
      Mean_Error = mean_error,
      Bias_Direction = ifelse(mean_error > 0, "Under-predicting", "Over-predicting"),
      Predictions = data.frame(
        Actual = actual, 
        Predicted = predictions,
        Residuals = actual - predictions
      )
    )
    
  } else {
    confusion_matrix <- confusionMatrix(predictions, actual)
    
    auc_value <- NA
    if (length(unique(actual)) == 2) {
      tryCatch({
        pred_probs <- attr(predict(model, test_data, probability = TRUE), "probabilities")[, 2]
        roc_obj <- roc(actual, pred_probs)
        auc_value <- auc(roc_obj)
      }, error = function(e) {
        auc_value <- NA
      })
    }
    
    metrics <- list(
      Confusion_Matrix = confusion_matrix,
      Accuracy = confusion_matrix$overall["Accuracy"],
      Kappa = confusion_matrix$overall["Kappa"],
      Sensitivity = confusion_matrix$byClass["Sensitivity"],
      Specificity = confusion_matrix$byClass["Specificity"],
      F1_Score = confusion_matrix$byClass["F1"],
      AUC = auc_value,
      Predictions = data.frame(Actual = actual, Predicted = predictions)
    )
  }
  
  # Return results
  result <- list(
    model = model,
    predictions = predictions,
    metrics = metrics,
    parameters = list(
      kernel = kernel,
      cost = cost,
      gamma = gamma,
      scale = scale,
      regression = regression
    ),
    data_info = list(
      train_size = nrow(train_data),
      test_size = nrow(test_data),
      target_variable = target_var
    )
  )
  
  return(result)
}
,
svm_rf_ensemble = function(train_data, test_data, target_var, input_vars, 
                            regression = TRUE, ntree=500) {
  
  library(e1071)
  library(randomForest)
  library(caret)
  
  # Train SVM model
  if (regression) {
    svm_model <- svm(
      x = train_data[, input_vars],
      y = train_data[[target_var]],
      type = "eps-regression",
      kernel = "radial"
    )
  } else {
    svm_model <- svm(
      x = train_data[, input_vars],
      y = train_data[[target_var]],
      type = "C-classification",
      kernel = "radial",
      probability = TRUE
    )
  }
  
  # Train Random Forest model
  rf_model <- randomForest(
    x = train_data[, input_vars],
    y = train_data[[target_var]],
    ntree = ntree,
    importance = TRUE
  )
  
  # Make predictions from both models
  svm_pred <- predict(svm_model, test_data[, input_vars])
  rf_pred <- predict(rf_model, test_data[, input_vars])
  
  # Simple average ensemble
  if (regression) {
    ensemble_pred <- (svm_pred + rf_pred) / 2
  } else {
    # For classification, use voting
    ensemble_pred <- ifelse((as.numeric(svm_pred) + as.numeric(rf_pred)) > 2.5, 
                            levels(train_data[[target_var]])[2],
                            levels(train_data[[target_var]])[1])
    ensemble_pred <- factor(ensemble_pred, levels = levels(train_data[[target_var]]))
  }
  
  # Calculate metrics
  actual <- test_data[[target_var]]
  
  if (regression) {
    residuals <- actual - ensemble_pred
    mae <- mean(abs(residuals))
    rmse <- sqrt(mean(residuals^2))
    r_squared <- cor(actual, ensemble_pred)^2
    
    metrics <- list(
      MAE = mae,
      RMSE = rmse,
      R_squared = r_squared,
      Ensemble_Predictions = data.frame(
        Actual = actual,
        SVM = svm_pred,
        RF = rf_pred,
        Ensemble = ensemble_pred
      )
    )
  } else {
    confusion_matrix <- confusionMatrix(ensemble_pred, actual)
    metrics <- list(
      Confusion_Matrix = confusion_matrix,
      Accuracy = confusion_matrix$overall["Accuracy"],
      Ensemble_Predictions = data.frame(
        Actual = actual,
        SVM = svm_pred,
        RF = rf_pred,
        Ensemble = ensemble_pred
      )
    )
  }
  
  return(list(
    svm_model = svm_model,
    rf_model = rf_model,
    predictions = ensemble_pred,
    metrics = metrics,
    individual_predictions = list(svm = svm_pred, rf = rf_pred)
  ))
}
,
suport_vector_regression_model = function(train_data, test_data, target_var, input_vars, 
                                           kernel = "radial", cost = 1, gamma = 0.1,
                                           epsilon = 0.1, scale = TRUE) {
  
  library(e1071)
  library(caret)
  
  # Create formula
  formula <- as.formula(paste(target_var, "~", paste(input_vars, collapse = " + ")))
  
  # Train SVR model
  model <- svm(
    formula,
    data = train_data,
    type = "eps-regression",  # This makes it SVR instead of SVM classification
    kernel = kernel,
    cost = cost,
    gamma = gamma,
    epsilon = epsilon,        # Epsilon parameter for SVR
    scale = scale
  )
  
  # Make predictions
  predictions <- predict(model, test_data)
  
  # Get actual values
  actual <- test_data[[target_var]]
  
  # Calculate comprehensive performance metrics
  valid_idx <- !is.na(actual) & !is.na(predictions)
  actual_clean <- actual[valid_idx]
  predicted_clean <- predictions[valid_idx]
  
  residuals <- actual_clean - predicted_clean
  mae <- mean(abs(residuals))
  rmse <- sqrt(mean(residuals^2))
  r_squared <- cor(actual_clean, predicted_clean)^2
  mape <- mean(abs(residuals / actual_clean)) * 100
  mean_error <- mean(residuals)
  median_ae <- median(abs(residuals))
  
  metrics <- list(
    MAE = mae,
    RMSE = rmse,
    R_squared = r_squared,
    MAPE = mape,
    Mean_Error = mean_error,
    Median_AE = median_ae,
    Bias_Direction = ifelse(mean_error > 0, "Under-predicting", "Over-predicting"),
    Predictions = data.frame(
      Actual = actual, 
      Predicted = predictions,
      Residuals = actual - predictions
    )
  )
  
  # Return results
  result <- list(
    model = model,
    predictions = predictions,
    metrics = metrics,
    parameters = list(
      kernel = kernel,
      cost = cost,
      gamma = gamma,
      epsilon = epsilon,
      scale = scale
    ),
    data_info = list(
      train_size = nrow(train_data),
      test_size = nrow(test_data),
      target_variable = target_var
    )
  )
  
  return(result)
}
,
gbm_model = function(train_data, test_data, target_var, input_vars, 
                      n_trees = 100, interaction_depth = 3, shrinkage = 0.1,
                      regression = TRUE) {
  
  library(gbm)
  library(caret)
  
  # Create formula
  formula <- as.formula(paste(target_var, "~", paste(input_vars, collapse = " + ")))
  
  # Train GBM model
  if (regression) {
    model <- gbm(
      formula,
      data = train_data,
      distribution = "gaussian",
      n.trees = n_trees,
      interaction.depth = interaction_depth,
      shrinkage = shrinkage,
      cv.folds = 5,
      verbose = FALSE
    )
    
    # Find optimal number of trees
    best_trees <- gbm.perf(model, method = "cv", plot.it = FALSE)
    
  } else {
    model <- gbm(
      formula,
      data = train_data,
      distribution = "bernoulli",
      n.trees = n_trees,
      interaction.depth = interaction_depth,
      shrinkage = shrinkage,
      cv.folds = 5,
      verbose = FALSE
    )
    
    best_trees <- gbm.perf(model, method = "cv", plot.it = FALSE)
  }
  
  # Make predictions
  predictions <- predict(model, test_data, n.trees = best_trees)
  
  # Get actual values
  actual <- test_data[[target_var]]
  
  # Calculate metrics
  if (regression) {
    valid_idx <- !is.na(actual) & !is.na(predictions)
    actual_clean <- actual[valid_idx]
    predicted_clean <- predictions[valid_idx]
    
    residuals <- actual_clean - predicted_clean
    mae <- mean(abs(residuals))
    rmse <- sqrt(mean(residuals^2))
    r_squared <- cor(actual_clean, predicted_clean)^2
    
    metrics <- list(
      MAE = mae,
      RMSE = rmse,
      R_squared = r_squared,
      Best_Trees = best_trees,
      Predictions = data.frame(Actual = actual, Predicted = predictions)
    )
  } else {
    # For classification, convert probabilities to classes
    class_predictions <- ifelse(predictions > 0.5, 
                                levels(actual)[2], 
                                levels(actual)[1])
    class_predictions <- factor(class_predictions, levels = levels(actual))
    
    confusion_matrix <- confusionMatrix(class_predictions, actual)
    metrics <- list(
      Confusion_Matrix = confusion_matrix,
      Accuracy = confusion_matrix$overall["Accuracy"],
      Predictions = data.frame(Actual = actual, 
                               Predicted = class_predictions,
                               Probability = predictions)
    )
  }
  
  # Get variable importance
  var_importance <- summary(model, plotit = FALSE)
  
  return(list(
    model = model,
    predictions = predictions,
    metrics = metrics,
    variable_importance = var_importance,
    parameters = list(
      n_trees = best_trees,
      interaction_depth = interaction_depth,
      shrinkage = shrinkage,
      regression = regression
    )
  ))
}
,
bart_model = function(train_data, test_data, target_var, input_vars,
                      n_tree = 50, n_iter = 1000, n_burn = 250) {
  
  library(BART)
  
  # Prepare data
  train_x <- train_data[, input_vars]
  train_y <- train_data[[target_var]]
  test_x <- test_data[, input_vars]
  
  # Convert to BART format
  x_train <- as.matrix(train_x)
  y_train <- as.numeric(train_y)
  x_test <- as.matrix(test_x)
  
  # Train BART model
  model <- wbart(x_train, y_train, x_test,
                 ntree = n_tree,
                 nskip = n_burn,
                 ndpost = n_iter)
  
  # Get predictions (posterior mean)
  predictions <- colMeans(model$yhat.test)
  
  # Calculate metrics
  actual <- test_data[[target_var]]
  valid_idx <- !is.na(actual) & !is.na(predictions)
  actual_clean <- actual[valid_idx]
  predicted_clean <- predictions[valid_idx]
  
  residuals <- actual_clean - predicted_clean
  mae <- mean(abs(residuals))
  rmse <- sqrt(mean(residuals^2))
  r_squared <- cor(actual_clean, predicted_clean)^2
  
  metrics <- list(
    MAE = mae,
    RMSE = rmse,
    R_squared = r_squared,
    Predictions = data.frame(Actual = actual, Predicted = predictions)
  )
  
  # Get variable importance
  var_importance <- colMeans(model$varcount)
  names(var_importance) <- input_vars
  
  return(list(
    model = model,
    predictions = predictions,
    metrics = metrics,
    variable_importance = var_importance,
    parameters = list(
      n_tree = n_tree,
      n_iter = n_iter,
      n_burn = n_burn
    )
  ))
}
,
elastic_net_model = function(train_data, test_data, target_var, input_vars,
                              alpha = 0.5, tune_lambda = TRUE) {
  
  library(glmnet)
  library(caret)
  
  # Prepare data
  train_x <- as.matrix(train_data[, input_vars])
  train_y <- train_data[[target_var]]
  test_x <- as.matrix(test_data[, input_vars])
  
  # Train Elastic Net with cross-validation
  if (tune_lambda) {
    cv_model <- cv.glmnet(train_x, train_y, alpha = alpha)
    best_lambda <- cv_model$lambda.min
  } else {
    cv_model <- glmnet(train_x, train_y, alpha = alpha)
    best_lambda <- cv_model$lambda[length(cv_model$lambda)]  # Most regularized
  }
  
  # Final model
  model <- glmnet(train_x, train_y, alpha = alpha, lambda = best_lambda)
  
  # Make predictions
  predictions <- predict(model, test_x)
  
  # Calculate metrics
  actual <- test_data[[target_var]]
  valid_idx <- !is.na(actual) & !is.na(predictions)
  actual_clean <- actual[valid_idx]
  predicted_clean <- predictions[valid_idx]
  
  residuals <- actual_clean - predicted_clean
  mae <- mean(abs(residuals))
  rmse <- sqrt(mean(residuals^2))
  r_squared <- cor(actual_clean, predicted_clean)^2
  
  metrics <- list(
    MAE = mae,
    RMSE = rmse,
    R_squared = r_squared,
    Predictions = data.frame(Actual = actual, Predicted = predictions)
  )
  
  # Get coefficients
  coefficients <- coef(model)
  
  return(list(
    model = model,
    predictions = predictions,
    metrics = metrics,
    coefficients = coefficients,
    parameters = list(
      alpha = alpha,
      lambda = best_lambda
    )
  ))
}
,
pls_model = function(train_data, test_data, target_var, input_vars,
                     ncomp = NULL, scale = TRUE) {
  
  library(pls)
  library(caret)
  
  # Create formula
  formula <- as.formula(paste(target_var, "~", paste(input_vars, collapse = " + ")))
  
  # Train PLS model
  model <- plsr(
    formula,
    data = train_data,
    scale = scale,
    validation = "CV"
  )
  
  # Find optimal components if not specified
  if (is.null(ncomp)) {
    ncomp <- which.min(RMSEP(model)$val[1,,]) - 1
  }
  
  # Make predictions
  predictions <- predict(model, test_data, ncomp = ncomp)[,,1]
  
  # Calculate metrics
  actual <- test_data[[target_var]]
  valid_idx <- !is.na(actual) & !is.na(predictions)
  actual_clean <- actual[valid_idx]
  predicted_clean <- predictions[valid_idx]
  
  residuals <- actual_clean - predicted_clean
  mae <- mean(abs(residuals))
  rmse <- sqrt(mean(residuals^2))
  r_squared <- cor(actual_clean, predicted_clean)^2
  
  metrics <- list(
    MAE = mae,
    RMSE = rmse,
    R_squared = r_squared,
    Best_Components = ncomp,
    Predictions = data.frame(Actual = actual, Predicted = predictions)
  )
  
  # Get variable importance (loadings)
  loadings <- model$loadings[, 1:ncomp]
  
  return(list(
    model = model,
    predictions = predictions,
    metrics = metrics,
    loadings = loadings,
    parameters = list(
      ncomp = ncomp,
      scale = scale
    )
  ))
}
  )
)