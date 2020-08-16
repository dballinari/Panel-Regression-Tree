# Panel regression tree
## Introduction
The code implements a decision tree with a fixed-effect panel regression model in each leaf. This 
approach extends the classical regression tree to panel structured data and allows to fit more 
complex models in the leaves.

The tree is grown in the same way as a classical regression tree, i.e. by minimizing the variance. The main
difference with respect to a normal regression tree, is that we have to ensure that each "individual" has at
least one observation in each node. The minimum number of observations for each individual can be controlled
through the input `min_ids`.

## Usage
The main function is `build_tree` which constructs and fits the fixed effect panel model. The user has to define a regression function `formula`, 
a vector of variables `split_variables` used for splitting the data, the name of the panel identifier `id`, the maximal depth of the tree `max_depth`, 
the minimal number of observations in each leaf `min_obs` and the minimal number of observations for each "individual" required in each leaf `min_ids`.

```R
soruce("panel_regression_tree.R")

# Fit the tree model
fitted <- build_tree(data = iris, formula = "Sepal.Length ~ Sepal.Width + Petal.Length", split_variables = c("Sepal.Width", "Petal.Length", "Petal.Width"), max_depth = 4, id = "Species", min_obs = 5, min_ids = 3)

# Print the tree
print_tree(tree = fitted)

# Prune the tree and print it
pruned <- prune_tree(tree = fitted)
print_tree(tree = fitted)

# Make predictions
predict_tree(tree = pruned, newdata = iris)
```

