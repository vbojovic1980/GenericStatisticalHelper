# GenericHelper R6 Class Documentation

A comprehensive R6 class providing statistical analysis, visualization, and machine learning methods for data science workflows.

## Statistical Analysis Methods

### correlations_and_regressions
Performs comprehensive correlation analysis and regression modeling between numeric variables with automated visualization and reporting.

**Arguments:**
- data: Input data frame
- outFolder: Output directory path
- filePrefix: Prefix for output files
- cor_threshold: Minimum correlation coefficient threshold
- p_threshold: Statistical significance threshold
- min_observations: Minimum observations required for analysis
- plot_width, plot_height, plot_res: Visualization dimensions and resolution

**Returns:** List containing correlation matrices, regression results, file paths, and analysis summary

### correlation_variable_selection
Performs correlation-based feature selection to remove highly correlated variables while preserving variance.

**Arguments:**
- df: Input data frame
- variables: Character vector of variable names to analyze
- correlation_method: Correlation calculation method ("pearson", "spearman", "kendall")
- threshold: Correlation threshold for variable removal
- alpha: Statistical significance level
- use_method: Missing value handling method

**Returns:** List with removed variables, remaining variables, similarity mappings, and selection summary

### shapiro_test_comprehensive
Performs Shapiro-Wilk normality tests on multiple variables with comprehensive descriptive statistics.

**Arguments:**
- data: Input data frame
- variables: Variables to test (defaults to all numeric variables)
- alpha: Significance level for normality assessment

**Returns:** Data frame with test statistics, p-values, normality classification, and descriptive statistics

### kruskal_wallis_comprehensive
Performs Kruskal-Wallis tests for multiple independent and dependent variable combinations with optional post-hoc analysis.

**Arguments:**
- data: Input data frame
- iv_vars: Independent (grouping) variables
- dv_vars: Dependent variables
- alpha: Significance level
- posthoc: Boolean for post-hoc test execution
- effect_size: Boolean for effect size calculation

**Returns:** Data frame with test statistics, p-values, group counts, and significance indicators

### dunn_tests
Performs Dunn's post-hoc tests following Kruskal-Wallis analysis for pairwise group comparisons.

**Arguments:**
- data: Input data frame
- dv_names: Dependent variable names
- iv_names: Independent variable names
- method: P-value adjustment method
- alpha: Significance level

**Returns:** Data frame with pairwise comparison results, z-values, and adjusted p-values

### permanova
Performs Permutational Multivariate Analysis of Variance (PERMANOVA) for multivariate data analysis.

**Arguments:**
- data: Input data frame
- dv_vars: Multivariate response variables
- iv_vars: Independent variables
- permutations: Number of permutations for significance testing

**Returns:** PERMANOVA results data frame with variance components and significance tests

## Data Visualization Methods

### create_group_boxplots_overview
Creates comprehensive boxplot overviews grouped by categorical variables with automated layout optimization.

**Arguments:**
- df: Input data frame
- output_folder: Output directory path
- group_var: Grouping variable name
- target_vars: Target variables for visualization
- dpi, width, height: Plot resolution and dimensions

**Returns:** Combined boxplot visualization with group legends

### create_variable_boxplots_overview
Generates overview of ungrouped boxplots for multiple variables in grid layout.

**Arguments:**
- df: Input data frame
- output_file: Output file path
- target_vars: Variables to visualize
- dpi, width, height: Plot specifications

**Returns:** Multi-panel boxplot visualization

### create_individual_boxplots
Creates individual boxplot files for each specified variable with grouping options.

**Arguments:**
- df: Input data frame
- output_folder: Output directory
- group_var: Optional grouping variable
- target_vars: Variables to plot
- dpi, width, height: Plot parameters

**Returns:** List of created file paths

### create_individual_histograms
Generates individual histogram files for specified variables with grouping capabilities.

**Arguments:**
- df: Input data frame
- output_folder: Output directory
- group_var: Optional grouping variable
- target_vars: Variables for histogram creation
- dpi, width, height, bins: Plot specifications

**Returns:** List of generated histogram file paths

### create_histograms_overview
Creates multi-panel histogram overview for rapid variable distribution assessment.

**Arguments:**
- df: Input data frame
- output_file: Output file path
- target_vars: Variables to include
- dpi, width, height, bins: Visualization parameters

**Returns:** Combined histogram grid visualization

### create_scatter_matrix
Generates scatter plot matrix with correlation coefficients and distribution diagnostics.

**Arguments:**
- data: Input data frame
- var_names: Variables to include in matrix
- output_path: Output file path
- correlation_method: Correlation calculation method
- Various styling and formatting parameters

**Returns:** Comprehensive scatter matrix plot

### create_correlation_plot
Creates correlation matrix visualization using corrplot package.

**Arguments:**
- data: Input data frame
- var_names: Variables for correlation analysis
- output_file: Output file path
- width_inches, height_inches, dpi: Plot dimensions

**Returns:** Correlation matrix plot and matrix object

### create_correlation_plot_sig
Generates significance-highlighted correlation matrix with statistical testing.

**Arguments:**
- data: Input data frame
- var_names: Variables to analyze
- output_file: Output file path
- correlation_method: Correlation calculation method
- Various significance and styling parameters

**Returns:** Enhanced correlation plot with significance indicators and summary statistics

### generate_correlation_plots_scatter
Creates individual scatter plots for significant correlations with regression lines.

**Arguments:**
- data: Input data frame
- variables: Variables to analyze
- correlation_method: Correlation calculation method
- pvalue, threshold: Significance criteria
- output_folder: Output directory
- Various plot customization parameters

**Returns:** Data frame of significant correlation results and generated plot files

## Machine Learning Methods

### generate_c50_tree
Implements C5.0 decision tree algorithm with comprehensive visualization and model assessment.

**Arguments:**
- data_frame: Input data frame
- input_var_names: Predictor variables
- target_var_name: Target variable
- output_png_file: Visualization output path
- width_inches, height_inches, dpi: Plot dimensions

**Returns:** List containing C5.0 and rpart models, accuracy metrics, and variable importance

### random_forest
Trains random forest models with error rate visualization and variable importance analysis.

**Arguments:**
- data_frame: Input data frame
- input_var_names: Predictor variables
- target_var_name: Target variable
- output_png_file: Output visualization path
- width_inches, height_inches, dpi: Plot specifications
- ntree: Number of trees in forest

**Returns:** Random forest model object with performance metrics

### random_forest_model
Comprehensive random forest implementation with detailed performance metrics and validation.

**Arguments:**
- train_data: Training dataset
- test_data: Testing dataset
- target_var: Target variable name
- input_vars: Predictor variables
- ntree: Number of trees
- mtry: Variables sampled at each split
- regression: Boolean for regression vs classification

**Returns:** Comprehensive model results including predictions, metrics, and variable importance

### xgboost_model
Extreme Gradient Boosting implementation with cross-validation and advanced performance tracking.

**Arguments:**
- train_data: Training dataset
- test_data: Testing dataset
- target_var: Target variable
- input_vars: Predictor variables
- nrounds: Boosting iterations
- params: XGBoost parameters
- regression: Regression or classification mode
- cv_folds: Cross-validation folds

**Returns:** XGBoost model with performance metrics, variable importance, and CV results

### dnn_model
Deep Neural Network implementation using Keras/TensorFlow backend with comprehensive configuration.

**Arguments:**
- train_data: Training dataset
- test_data: Testing dataset
- target_var: Target variable
- input_vars: Predictor variables
- layers: Hidden layer architecture
- activation: Activation functions
- dropout: Regularization parameters
- Various training hyperparameters

**Returns:** DNN model with training history, predictions, metrics, and variable importance

### knn_model
K-Nearest Neighbors implementation with automatic parameter tuning and scaling options.

**Arguments:**
- train_data: Training dataset
- test_data: Testing dataset
- target_var: Target variable
- input_vars: Predictor variables
- k: Number of neighbors
- scale_data: Feature scaling flag
- regression: Regression or classification mode
- tune_k: Automatic parameter tuning

**Returns:** KNN model results with performance metrics and tuning information

### extra_trees_ranger
Extremely Randomized Trees implementation using ranger package with advanced tuning capabilities.

**Arguments:**
- train_data: Training dataset
- test_data: Testing dataset
- target_var: Target variable
- input_vars: Predictor variables
- ntree: Number of trees
- mtry: Variable sampling parameter
- min_node_size: Tree stopping criterion
- Various algorithm-specific parameters

**Returns:** Extra Trees model with comprehensive performance assessment

### svm_model
Support Vector Machine implementation with multiple kernel options and comprehensive metrics.

**Arguments:**
- train_data: Training dataset
- test_data: Testing dataset
- target_var: Target variable
- input_vars: Predictor variables
- kernel: SVM kernel type
- cost, gamma: Model parameters
- scale: Feature scaling flag
- regression: Regression or classification mode

**Returns:** SVM model with performance metrics and predictions

### suport_vector_regression_model
Support Vector Regression implementation with epsilon-insensitive loss function.

**Arguments:**
- train_data: Training dataset
- test_data: Testing dataset
- target_var: Target variable
- input_vars: Predictor variables
- kernel: SVR kernel type
- cost, gamma, epsilon: Model parameters
- scale: Feature scaling flag

**Returns:** SVR model with regression metrics and residual analysis

### svm_rf_ensemble
Hybrid ensemble model combining Support Vector Machine and Random Forest predictions.

**Arguments:**
- train_data: Training dataset
- test_data: Testing dataset
- target_var: Target variable
- input_vars: Predictor variables
- regression: Problem type flag
- ntree: Number of trees in RF component

**Returns:** Ensemble model with individual and combined predictions

## Data Preprocessing Methods

### minmax_normalize
Performs min-max normalization on numeric vectors to [0,1] range.

**Arguments:**
- x: Numeric vector for normalization

**Returns:** Normalized vector with preserved structure

### normalized_variance
Calculates variance on min-max normalized data for scale-invariant dispersion measurement.

**Arguments:**
- x: Numeric vector for analysis

**Returns:** Normalized variance value

### binarize_variables_by_mean
Converts continuous variables to binary based on mean thresholding.

**Arguments:**
- data_frame: Input data frame
- var_names: Variables to binarize
- suffix: Naming suffix for new variables
- na_action: Missing value handling strategy

**Returns:** Data frame with added binary variables and summary statistics

### filter_variables_usingXY
Filters predictor variables based on correlation with target variables and significance.

**Arguments:**
- data: Input data frame
- y_vars: Target variables
- x_vars: Predictor candidates
- cor_threshold: Correlation threshold
- pvalue_threshold: Statistical significance threshold

**Returns:** Filtered variable list with removal reasons and summary

## Specialized Analytical Methods

### hierarchical_clustering_plot
Performs hierarchical clustering with dendrogram visualization and cluster analysis.

**Arguments:**
- data_frame: Input data frame
- var_names: Variables for clustering
- data_label_var: Observation labeling variable
- distance_method: Distance calculation method
- Various clustering and visualization parameters

**Returns:** Clustering results, dendrogram, and cluster assignments

### get_fingerprint_maccs_binary
Converts SMILES strings to MACCS molecular fingerprints in binary representation.

**Arguments:**
- smiles: SMILES string input

**Returns:** Binary fingerprint representation with structural information

### KNN_label
K-Nearest Neighbors classification implementation with probability estimates.

**Arguments:**
- train_data: Training dataset
- test_data: Testing dataset
- target_var: Classification target
- predictor_vars: Feature variables
- k: Number of neighbors
- l: Decision threshold parameter

**Returns:** Classification results with probabilities and performance metrics

### KNN_reg
K-Nearest Neighbors regression implementation with comprehensive error metrics.

**Arguments:**
- train_data: Training dataset
- test_data: Testing dataset
- target_var: Regression target
- predictor_vars: Feature variables
- k: Number of neighbors

**Returns:** Regression results with predictions and error analysis

### regression_models
Performs automated regression analysis on highly correlated variable pairs.

**Arguments:**
- data: Input data frame
- var_names: Variables for analysis
- threshold: Correlation threshold
- cor_method: Correlation calculation method

**Returns:** Data frame of regression results with coefficients and significance

### create_regression_plot
Creates detailed scatter plots with regression lines and statistical annotations.

**Arguments:**
- data: Input data frame
- var_x: Independent variable
- var_y: Dependent variable
- output_file: Output path
- Plot dimension parameters

**Returns:** Regression visualization with model statistics and equation