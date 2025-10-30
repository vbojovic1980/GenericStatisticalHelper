# GenericStatisticalHelper
Just automatizing some statistical methods

# GenericHelper R6 Class

A comprehensive R6 class providing statistical analysis, visualization, and machine learning methods for data science workflows.

## Methods Overview

### 1. Correlation & Regression Analysis

#### `correlations_and_regressions()`
Comprehensive correlation analysis with visualization and regression modeling.

**Parameters:**
- `data`: Data frame with numeric variables
- `outFolder`: Output directory for results
- `filePrefix`: Prefix for output files
- `cor_threshold`: Correlation threshold
- `p_threshold`: P-value threshold
- `min_observations`: Minimum observations required
- `plot_width`, `plot_height`, `plot_res`: Plot dimensions and resolution

**Results:**
- Correlation matrices (CSV)
- P-value matrices (CSV)
- Regression results (CSV)
- Multiple visualization plots (correlation matrix, heatmap, network plot, regression plots)

#### `correlation_variable_selection()`
Feature selection based on correlation analysis.

**Parameters:**
- `df`: Data frame containing variables
- `variables`: Character vector of variable names
- `correlation_method`: "pearson", "spearman", or "kendall"
- `threshold`: Correlation threshold
- `alpha`: Significance level
- `use_method`: Method for handling missing values

**Results:**
- `removed_vars`: Variables removed due to high correlation
- `remaining_vars`: Variables retained
- `similar_vars`: Mapping of removed variables to similar ones
- Summary statistics of reduction

#### `generate_correlation_plots_scatter()`
Generate scatter plots for significant correlations.

**Parameters:**
- `data`: Input data frame
- `variables`: Variables to analyze
- `correlation_method`: Correlation method
- `pvalue`: P-value threshold
- `threshold`: Correlation threshold
- `output_folder`: Output directory
- `width`, `height`, `dpi`: Plot specifications
- `smooth_method`: Smoothing method for plots

**Results:**
- Individual scatter plots for significant correlations
- Results data frame with correlation statistics

### 2. Visualization Methods

#### `create_group_boxplots_overview()`
Create overview of boxplots grouped by a categorical variable.

**Parameters:**
- `df`: Input data frame
- `output_folder`: Output directory
- `group_var`: Grouping variable
- `target_vars`: Target variables to plot
- `dpi`, `width`, `height`: Plot specifications

**Results:**
- Combined boxplot overview PNG file
- Grid layout optimized for number of variables

#### `create_variable_boxplots_overview()`
Create overview of ungrouped boxplots for multiple variables.

**Parameters:**
- `df`: Input data frame
- `output_file`: Output file path
- `target_vars`: Target variables to plot
- `dpi`, `width`, `height`: Plot specifications

**Results:**
- Combined boxplot overview PNG file

#### `create_individual_boxplots()`
Create individual boxplot files for each variable.

**Parameters:**
- `df`: Input data frame
- `output_folder`: Output directory
- `group_var`: Optional grouping variable
- `target_vars`: Target variables
- `dpi`, `width`, `height`: Plot specifications

**Results:**
- Individual PNG files for each variable
- List of created file paths

#### `create_individual_histograms()`
Create individual histogram files for each variable.

**Parameters:**
- `df`: Input data frame
- `output_folder`: Output directory
- `group_var`: Optional grouping variable
- `target_vars`: Target variables
- `dpi`, `width`, `height`: Plot specifications
- `bins`: Number of histogram bins

**Results:**
- Individual histogram PNG files
- List of created file paths

#### `create_histograms_overview()`
Create overview of histograms for multiple variables.

**Parameters:**
- `df`: Input data frame
- `output_file`: Output file path
- `target_vars`: Target variables
- `dpi`, `width`, `height`: Plot specifications
- `bins`: Number of histogram bins

**Results:**
- Combined histogram overview PNG file

#### `create_scatter_matrix()`
Create scatter plot matrix with correlation coefficients.

**Parameters:**
- `data`: Input data frame
- `var_names`: Variables to include
- `output_path`: Output file path
- `correlation_method`: Correlation method
- `resolution`, `width`, `height`: Plot specifications
- Various text size and display parameters

**Results:**
- Scatter matrix PNG file with correlations
- Invisible plot object

#### `create_correlation_plot()`
Create correlation matrix visualization.

**Parameters:**
- `data`: Input data frame
- `var_names`: Variables to include
- `output_file`: Output file path
- `width_inches`, `height_inches`, `dpi`: Plot specifications

**Results:**
- Correlation plot PNG file
- Correlation matrix

#### `create_correlation_plot_sig()`
Create correlation matrix with significance indicators.

**Parameters:**
- `data`: Input data frame
- `var_names`: Variables to include
- `output_file`: Output file path
- `correlation_method`: Correlation method
- `width_inches`, `height_inches`, `dpi`: Plot specifications
- `alpha`: Significance level
- Various display parameters for non-significant correlations

**Results:**
- Enhanced correlation plot with significance
- Comprehensive statistics and matrices

### 3. Statistical Tests

#### `regression_models()`
Perform regression analysis for highly correlated variable pairs.

**Parameters:**
- `data`: Input data frame
- `var_names`: Variables to analyze
- `threshold`: Correlation threshold
- `cor_method`: Correlation method

**Results:**
- Data frame with regression results
- Correlation coefficients and p-values
- Regression coefficients and statistics

#### `create_regression_plot()`
Create individual regression plot for variable pairs.

**Parameters:**
- `data`: Input data frame
- `var_x`, `var_y`: Variables for regression
- `output_file`: Output file path
- `width_inches`, `height_inches`, `dpi`: Plot specifications

**Results:**
- Regression plot PNG file
- Regression statistics and model

#### `shapiro_test_comprehensive()`
Comprehensive normality testing using Shapiro-Wilk test.

**Parameters:**
- `data`: Input data frame
- `variables`: Variables to test (default: all numeric)
- `alpha`: Significance level

**Results:**
- Data frame with normality test results
- Additional statistics (mean, SD, skewness, kurtosis)
- Normality classification

#### `kruskal_wallis_comprehensive()`
Comprehensive Kruskal-Wallis non-parametric ANOVA.

**Parameters:**
- `data`: Input data frame
- `iv_vars`: Independent (grouping) variables
- `dv_vars`: Dependent variables
- `alpha`: Significance level
- `posthoc`: Whether to perform post-hoc tests
- `effect_size`: Whether to calculate effect size

**Results:**
- Data frame with test results
- Optional post-hoc comparisons and effect sizes
- Significance indicators

#### `dunn_tests()`
Perform Dunn's post-hoc tests for Kruskal-Wallis.

**Parameters:**
- `data`: Input data frame
- `dv_names`: Dependent variables
- `iv_names`: Independent variables
- `method`: P-value adjustment method
- `alpha`: Significance level

**Results:**
- Data frame with Dunn's test results
- Pairwise comparisons with adjusted p-values
- Kruskal-Wallis statistics

#### `permanova()`
PERMANOVA for multivariate analysis.

**Parameters:**
- `data`: Input data frame
- `dv_vars`: Dependent variables (multivariate)
- `iv_vars`: Independent variables
- `permutations`: Number of permutations

**Results:**
- PERMANOVA results data frame
- Multivariate F-statistics and p-values

### 4. Clustering & Machine Learning

#### `hierarchical_clustering_plot()`
Hierarchical clustering with visualization.

**Parameters:**
- `data_frame`: Input data frame
- `var_names`: Variables for clustering
- `data_label_var`: Variable for labels
- `distance_method`: Distance metric
- `output_file`: Output file path
- `width`, `height`, `dpi`: Plot specifications
- `linkage_method`: Clustering linkage method
- `num_clusters`: Number of clusters (auto-detected if NULL)
- Various display parameters

**Results:**
- Dendrogram PNG file
- Cluster assignments and model objects

#### `binarize_variables_by_mean()`
Binarize numeric variables based on mean threshold.

**Parameters:**
- `data_frame`: Input data frame
- `var_names`: Variables to binarize
- `suffix`: Suffix for new binary variables
- `na_action`: How to handle NA values

**Results:**
- Modified data frame with binary variables
- Means summary as attribute

#### `generate_c50_tree()`
Generate C5.0 decision tree with visualization.

**Parameters:**
- `data_frame`: Input data frame
- `input_var_names`: Predictor variables
- `target_var_name`: Target variable
- `output_png_file`: Output file path
- `width_inches`, `height_inches`, `dpi`: Plot specifications

**Results:**
- Decision tree PNG file
- C5.0 and rpart models
- Accuracy and variable importance

#### `random_forest()`
Random Forest implementation with visualization.

**Parameters:**
- `data_frame`: Input data frame
- `input_var_names`: Predictor variables
- `target_var_name`: Target variable
- `output_png_file`: Output file path
- `width_inches`, `height_inches`, `dpi`: Plot specifications
- `ntree`: Number of trees

**Results:**
- Random Forest model
- Error rate and variable importance plots
- OOB error statistics

#### `KNN_label()`
K-Nearest Neighbors for classification.

**Parameters:**
- `train_data`: Training data
- `test_data`: Test data
- `target_var`: Target variable
- `predictor_vars`: Predictor variables
- `k`: Number of neighbors
- `l`: Minimum vote threshold

**Results:**
- Predictions and probabilities
- Accuracy and confusion matrix
- Comprehensive classification metrics

#### `KNN_reg()`
K-Nearest Neighbors for regression.

**Parameters:**
- `train_data`: Training data
- `test_data`: Test data
- `target_var`: Target variable
- `predictor_vars`: Predictor variables
- `k`: Number of neighbors

**Results:**
- Predictions and residuals
- Regression metrics (MSE, RMSE, MAE, R-squared)
- Comprehensive regression results

### 5. Utility Methods

#### `minmax_normalize()`
Min-max normalization of a vector.

**Parameters:**
- `x`: Numeric vector to normalize

**Results:**
- Normalized vector between 0 and 1

#### `normalized_variance()`
Calculate variance of min-max normalized data.

**Parameters:**
- `x`: Numeric vector

**Results:**
- Variance of normalized data

#### `filter_variables_usingXY()`
Filter variables based on correlation with target variables.

**Parameters:**
- `data`: Input data frame
- `y_vars`: Target variables
- `x_vars`: Predictor variables to filter
- `cor_threshold`: Correlation threshold
- `pvalue_threshold`: P-value threshold

**Results:**
- Lists of removed and remaining variables
- Filtering criteria summary

#### `get_fingerprint_maccs_binary()**
Generate MACCS fingerprints from SMILES strings.

**Parameters:**
- `smiles`: SMILES string

**Results:**
- Binary fingerprint string and vector
- Fingerprint length information