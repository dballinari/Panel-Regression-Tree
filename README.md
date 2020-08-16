# Panel regression decision tree
## Introduction
The code implements a decision tree with a fixed-effect panel regression model in each leaf. This 
approach extends the classical regression tree to panel structured data and allows to fit more 
complex models in the leaves.

The tree is grown in the same way as a classical regression tree, i.e. by minimizing the variance. The main
difference with respect to a normal regression tree, is that we have to ensure that each "individual" has at
least one observation in each node. The minimum number of observations for each individual can be controlled
through the input `min_ids`.

## Usage
First, a `GoogleTrendsScraper` object needs to be initialized with (optional) parameters, such as
the `sleep` time used to when long time windows are scraped, the `headless` option which defines
whether the browser opened by Selenium is visible. In a second step, the trends for a specific 
keyword, time range and region can be obtained by running the `GoogleTrendsScrpaer.get_trends`
method. 

```R
soruce("panel_regression_tree.R)

fitted <- build_tree(data = iris, formula = "Sepal.Length ~ Sepal.Width + Petal.Length", 
split_variables = c("Sepal.Width", "Petal.Length", "Petal.Width"), max_depth = 4, id = "Species", min_obs = 5, min_ids = 3)
```

```R
data <- tibble(a = 1:5, b = 6:10)
data %>% export_data_to_latex(path="data_for_latex.tex")
```
