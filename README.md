# GenericStatisticalHelper
Just automatizing some statistical methods

# GenericHelper Class Documentation

## Overview
The `GenericHelper` class provides comprehensive statistical analysis and visualization methods for data science workflows in R.
## Example usage
```R
source("https://raw.githubusercontent.com/vbojovic1980/GenericStatisticalHelper/main/GenericHelper.R")
gh=GenericHelper$new()
gh$someMethod(someParams)
```
## Methods

### `correlations_and_regressions()`
Performs comprehensive correlation analysis and regression modeling on numeric variables.

**Features:**
- Calculates pairwise correlations with p-values
- Generates correlation matrix visualizations (corrplot and ggplot2 heatmap)
- Performs linear regression on significant correlations
- Creates regression plots for strong correlations (r ≥ 0.8)
- Exports results to CSV files (correlations, p-values, regression results)
- Generates network plots of significant correlations

**Parameters:**
- `data`: Data frame containing the data
- `outFolder`: Output directory for results
- `filePrefix`: Prefix for output files
- `cor_threshold`: Correlation threshold (default: 0.7)
- `p_threshold`: P-value threshold (default: 0.05)
- `min_observations`: Minimum observations required (default: 3)
- `plot_width`, `plot_height`, `plot_res`: Plot dimensions and resolution

---

### `minmax_normalize()`
Performs min-max normalization on a numeric vector.

**Features:**
- Scales values to [0,1] range
- Handles constant variables (returns zeros)
- Robust to missing values

**Parameters:**
- `x`: Numeric vector to normalize

**Returns:** Normalized vector

---

### `normalized_variance()`
Calculates variance of min-max normalized data.

**Features:**
- Useful for comparing variability across different scales
- Handles missing values

**Parameters:**
- `x`: Numeric vector

**Returns:** Variance of normalized data

---

### `correlation_variable_selection()`
Performs correlation-based variable selection to remove highly correlated variables.

**Features:**
- Identifies variable pairs with correlation above threshold
- Removes variables with lower normalized variance
- Considers statistical significance of correlations
- Provides detailed reporting of removed variables and reasons

**Parameters:**
- `df`: Data frame containing variables
- `variables`: Character vector of variable names to analyze
- `threshold`: Correlation threshold (default: 0.8)
- `alpha`: Significance level (default: 0.05)

**Returns:** List with removed variables, remaining variables, and similarity information

---

### `create_group_boxplots_overview()`
Creates a comprehensive overview of boxplots grouped by a categorical variable.

**Features:**
- Generates grid of boxplots for multiple target variables
- Automatic grid layout optimization
- Shared legend and group level information
- Enhanced visual styling with jitter points
- Customizable dimensions and resolution

**Parameters:**
- `df`: Input data frame
- `output_folder`: Directory to save the plot
- `group_var`: Grouping variable name
- `target_vars`: Vector of target variable names
- `dpi`, `width`, `height`: Plot resolution and dimensions

---

### `create_variable_boxplots_overview()`
Creates an overview of ungrouped boxplots for multiple variables.

**Features:**
- Single-panel boxplots arranged in grid
- Automatic layout optimization
- Consistent styling across variables
- Customizable output dimensions

**Parameters:**
- `df`: Input data frame
- `output_folder`: Output directory
- `target_vars`: Variables to plot
- `dpi`, `width`, `height`: Plot specifications

---

### `create_individual_boxplots()`
Creates individual boxplot files for each variable.

**Features:**
- Optionally grouped by categorical variable
- High-quality individual plots
- Enhanced dot visibility
- Flexible output customization

**Parameters:**
- `df`: Input data frame
- `output_folder`: Output directory
- `group_var`: Optional grouping variable
- `target_vars`: Variables to plot
- `dpi`, `width`, `height`: Plot specifications

---

### `create_individual_histograms()`
Creates individual histogram files for each variable.

**Features:**
- Optionally grouped by categorical variable
- Customizable bin count
- High-quality output
- Consistent styling

**Parameters:**
- `df`: Input data frame
- `output_folder`: Output directory
- `group_var`: Optional grouping variable
- `target_vars`: Variables to plot
- `dpi`, `width`, `height`, `bins`: Plot specifications

---

### `create_histograms_overview()`
Creates a grid overview of histograms for multiple variables.

**Features:**
- Compact grid layout
- Automatic dimension calculation
- Consistent styling
- Efficient visualization of multiple distributions

**Parameters:**
- `df`: Input data frame
- `output_folder`: Output directory
- `target_vars`: Variables to plot
- `dpi`, `width`, `height`, `bins`: Plot specifications

---

### `create_scatter_matrix()`
Creates a scatter plot matrix (pairs plot) for multiple variables.

**Features:**
- Correlation coefficients in upper triangle
- Scatter plots in lower triangle
- Density plots on diagonal
- Fully customizable text sizes
- High-resolution output

**Parameters:**
- `data`: Input data
- `var_names`: Variables to include
- `output_path`: Output file path
- `resolution`, `width`, `height`: Output specifications
- Various text size controls

---

### `create_correlation_plot()`
Creates a correlation matrix visualization using corrplot.

**Features:**
- Hierarchical clustering of variables
- Color-coded correlation values
- Coefficient display
- Customizable dimensions

**Parameters:**
- `data`: Input data
- `var_names`: Variables to correlate
- `output_file`: Output path
- `width_inches`, `height_inches`, `dpi`: Output specifications

---

### `create_correlation_plot_sig()`
Creates a correlation plot highlighting statistically significant correlations.

**Features:**
- Marks non-significant correlations with crosses
- Customizable significance level
- Summary statistics of significant correlations
- Enhanced visual distinction

**Parameters:**
- `data`: Input data
- `var_names`: Variables to correlate
- `output_file`: Output path
- `width_inches`, `height_inches`, `dpi`: Output specifications
- `alpha`: Significance level

---

### `analyze_strong_correlations()`
Identifies and analyzes strongly correlated variable pairs.

**Features:**
- Calculates correlation coefficients and p-values
- Performs linear regression on strong correlations
- Returns comprehensive results table
- Sorted by correlation strength

**Parameters:**
- `data`: Input data
- `var_names`: Variables to analyze
- `threshold`: Correlation threshold (default: 0.7)

**Returns:** Data frame with correlation and regression results

---

### `create_regression_plot()`
Creates a detailed scatter plot with regression line.

**Features:**
- Linear regression with equation
- R-squared and p-value display
- High-quality visualization
- Statistical summary

**Parameters:**
- `data`: Input data
- `var_x`, `var_y`: X and Y variables
- `output_file`: Output path
- `width_inches`, `height_inches`, `dpi`: Output specifications

---

### `shapiro_test_comprehensive()`
Performs Shapiro-Wilk normality tests with comprehensive statistics.

**Features:**
- Tests multiple variables for normality
- Calculates additional statistics (mean, SD, skewness, kurtosis)
- Handles insufficient data gracefully
- Sorted by p-value

**Parameters:**
- `data`: Input data
- `variables`: Variables to test (default: all numeric)
- `alpha`: Significance level (default: 0.05)

**Returns:** Data frame with test results and statistics

---

### `kruskal_wallis_comprehensive()`
Performs Kruskal-Wallis tests with optional post-hoc analysis.

**Features:**
- Tests multiple IV-DV combinations
- Optional effect size calculation (epsilon squared)
- Optional post-hoc tests with PMCMRplus
- Comprehensive results reporting

**Parameters:**
- `data`: Input data
- `iv_vars`: Independent variables (factors)
- `dv_vars`: Dependent variables (numeric)
- `alpha`: Significance level
- `posthoc`: Whether to perform post-hoc tests
- `effect_size`: Whether to calculate effect size

**Returns:** Data frame with test results

---

### `dunn_tests()`
Performs Dunn's post-hoc tests following Kruskal-Wallis.

**Features:**
- Automatic test for significant Kruskal-Wallis results
- Multiple comparison correction
- Detailed pairwise comparisons
- Comprehensive results structure

**Parameters:**
- `data`: Input data
- `dv_names`: Dependent variables
- `iv_names`: Independent variables
- `method`: P-value adjustment method
- `alpha`: Significance level

**Returns:** Data frame with Dunn's test results

---

### `permanova()`
Performs PERMANOVA (Permutational Multivariate Analysis of Variance).

**Features:**
- Multivariate analysis using Bray-Curtis distance
- Handles multiple independent variables
- Comprehensive results table
- Robust to missing data

**Parameters:**
- `data`: Input data
- `dv_vars`: Dependent variables (multivariate response)
- `iv_vars`: Independent variables
- `permutations`: Number of permutations

**Returns:** PERMANOVA results data frame

---

### `hierarchical_clustering_plot()`
Creates hierarchical clustering visualization with dendrogram.

**Features:**
- Multiple distance and linkage methods
- Automatic cluster number suggestion
- Customizable coloring and labeling
- Horizontal or vertical orientation
- High-quality output

**Parameters:**
- `data_frame`: Input data
- `var_names`: Variables for clustering
- `data_label_var`: Variable for labeling
- Various clustering and output options

**Returns:** Clustering results and assignments

---

### `binarize_variables_by_mean()`
Binarizes numeric variables based on their mean values.

**Features:**
- Values above mean become 1, below become 0
- Multiple NA handling strategies
- Summary statistics output
- Flexible naming convention

**Parameters:**
- `data_frame`: Input data
- `var_names`: Variables to binarize
- `suffix`: Suffix for new variable names
- `na_action`: How to handle NAs

**Returns:** Data frame with binarized variables

---

### `generate_c50_tree()`
Generates C5.0 decision tree with visualization.

**Features:**
- Handles special characters in variable names
- Creates rpart visualization
- Variable importance analysis
- Model performance metrics

**Parameters:**
- `data_frame`: Input data
- `input_var_names`: Predictor variables
- `target_var_name`: Target variable
- Output specifications

**Returns:** Model objects and analysis results

---

### `get_fingerprint_maccs_binary()`
Generates MACCS fingerprint as binary string from SMILES.

**Features:**
- Converts SMILES to molecular structure
- Generates 166-bit MACCS fingerprint
- Returns binary string and vector representations

**Parameters:**
- `smiles`: SMILES string

**Returns:** Fingerprint representation

---

### `random_forest()`
Builds random forest model with visualization.

**Features:**
- Robust handling of special characters in variable names
- Error rate and variable importance plots
- Out-of-bag error estimation
- High-quality visualization output

**Parameters:**
- `data_frame`: Input data
- `input_var_names`: Predictor variables
- `target_var_name`: Target variable
- Output and model specifications

**Returns:** Random forest model object

---

### `filter_variables_usingXY()`
Filters variables based on correlation with target variables.

**Features:**
- Removes zero-variance variables
- Filters based on correlation threshold and significance
- Comprehensive reporting of filtering process
- Handles both numeric and non-numeric variables

**Parameters:**
- `data`: Input data
- `y_vars`: Target variables
- `x_vars`: Predictor variables to filter
- `cor_threshold`: Minimum correlation
- `pvalue_threshold`: Maximum p-value

**Returns:** Filtering results and remaining variables