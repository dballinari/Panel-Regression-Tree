# Title: R implementation of a fixed effects regression tree model
# Date: 27.07.2020
# Author: Daniele Ballinari
# R version 4.0.1 (2020-06-06)

if (!require("fixest")) {
  install.packages("fixest")
}


# Evaluate the loss of the model
model_loss <- function(data, formula, id) {

  model <- feols(fml=formula, data = data)
  
  sse <- sum(residuals(model)^2)
  
  return(sse)
}

# Function that evaluates the split
evaluate_split <- function(groups, formula, id) {
  sse_left <- model_loss(groups[["left"]], formula, id)
  sse_right <- model_loss(groups[["right"]], formula, id)
  return(list("left"=sse_left, "right"=sse_right))
}

# Split a dataset based on a specific column and a value
split_dataset <- function(data, variable, value, id, min_obs, min_ids) {
  # Determine the split
  left_index <- data[,variable]<value 
  
  # Divide dataset 
  left <- data[left_index,]
  right <- data[!left_index, ]
  
  # Check if split is valid:
  if (length(unique(left[,id])) < length(unique(data[,id])) | 
      nrow(left) < min_obs |
      min(table(left[,id])) < min_ids) left <- NULL
  if (length(unique(right[,id])) < length(unique(data[,id])) | 
      nrow(right) < min_obs |
      min(table(right[,id])) < min_ids) right <- NULL
  
  return(list("left"=left, "right"=right))
}

# Get the best split for a dataset 
get_split <- function(data, formula, split_variables, id, min_obs, min_ids, mesh = 8) {
  
  best_loss <- best_loss_left <- best_loss_right <- Inf
  best_index <- NULL
  best_value <- NULL
  best_variable <- NULL
  best_groups <- NULL
  
  alphas <- (1:(mesh-1))/mesh
  
  for (variable in split_variables) {
    values <- quantile(tdigest::tdigest(data[,variable]), alphas)
    for (value in values) {
      data_groups <- split_dataset(data, variable, value, id, min_obs, min_ids)
      
      if (is.null(data_groups[["left"]]) | is.null(data_groups[["right"]])) next
      
      losses <- evaluate_split(data_groups, formula, id)
      loss <- losses[["left"]] + losses[["right"]]
      
      if (loss < best_loss) {
        best_loss <- loss
        best_loss_left <- losses[["left"]]
        best_loss_right <- losses[["right"]]
        best_variable <- variable
        best_value <- value
        best_groups <- data_groups
      }
    }
  }
  return(list("variable"=best_variable, 
              "value"=best_value,
              "groups"=best_groups,
              "loss_left"=best_loss_left,
              "loss_right"=best_loss_right))
}

# Return the predictive model for the final leaf:
to_terminal <- function(data, formula, id) {

  model <- feols(fml=formula, data)
  
  return(model)
}

# Main recursive algorithm that sequentially splits the data in regions
split <- function(node, formula, split_variables, id, num_ids,
                  max_depth, min_obs, min_ids, depth, node_num=0) {
  if (is.null(node[["groups"]])) {
    return(NULL)
  }
  left <- node[["groups"]][["left"]]
  right <- node[["groups"]][["right"]]
  node$groups <- NULL
  # Add models
  node$model_left <- to_terminal(left, formula, id)
  node$model_right <- to_terminal(right, formula, id)
  node$id <- node_num
  # Check if max depth is reached 
  if (depth >= max_depth) {
    return(list("node"=node, "node_num"=node_num))
  }
  # Process left child:
  if (nrow(left) >= (min_obs*2) & 
      min(table(left[,id])) >= (2*min_ids) & 
      length(unique(left[,id])) == num_ids ) {
    new_split <- split(get_split(left, formula, split_variables, id, min_obs, min_ids), 
                       formula, split_variables, id, num_ids, max_depth, min_obs, min_ids, 
                       depth+1, node_num+1)
    if (is.null(new_split)) {
      return(list("node"=node, "node_num"=node_num))
    }
    node$left <- new_split$node
    node_num <- new_split$node_num
  }
  # Process right child:
  if (nrow(right) >= (min_obs*2) & 
      min(table(right[,id])) >= 2 & 
      length(unique(right[,id])) == num_ids) {
    new_split <- split(get_split(right, formula, split_variables, id, min_obs, min_ids), 
                       formula, split_variables, id, num_ids, max_depth, min_obs, min_ids,
                       depth+1, node_num+1)
    if (is.null(new_split)) {
      return(list("node"=node, "node_num"=node_num))
    }
    node$right <- new_split$node
    node_num <- new_split$node_num
  }
  return(list("node"=node, "node_num"=node_num))
}

# Main function that starts the tree-growing procedure
build_tree <- function(data, formula, split_variables, id, max_depth, min_obs, min_ids = 1) {
  # Check that data and formula are in the correct format
  data <- as.data.frame(data)
  formula <- as.formula(paste(formula, id, sep="|"))
  # Make sure that the variable to be predicted is in the first column
  base_loss <- model_loss(data, formula, id)
  # Number of unique ids
  num_ids <- length(unique(data[,id]))
  # If the maximum depth is zero or the dataset is to small, return tree without splits
  if (max_depth == 0 | nrow(data)<(2*min_obs) ) {
    root <- list()
  } else {
    root <- split(get_split(data, formula, split_variables, id, min_obs, min_ids), 
                  formula, split_variables, id, num_ids, max_depth, min_obs, min_ids,
                  1, 1) 
  }
  root$base_loss <- base_loss
  root$base_model <- to_terminal(data, formula, id)
  root$nobs <- nrow(data)
  root$nids <- num_ids
  root$variables <- colnames(data)
  return(root)
}

# Function that returns the total sum o squared residuals of the tree
get_tree_loss <- function(tree, loss=0) {
  # If the input is the base tree, extract the first node to start the recursion
  if (!is.null(tree[["base_loss"]])) {
    # If there is no node, i.e. the tree as a depth of 0, we return the base loss
    if (is.null(tree[["node"]])) {
      return(tree[["base_loss"]])
    }
    tree <- tree[["node"]]
  }
  # Handle the left child
  if (is.null(tree[["left"]])) {
    loss <- loss + tree[["loss_left"]]
  } else {
    loss <- get_tree_loss(tree[["left"]], loss)
  }
  # Handle the right child
  if (is.null(tree[["right"]])) {
    loss <- loss + tree[["loss_right"]]
  } else {
    loss <- get_tree_loss(tree[["right"]], loss)
  }
  
  return(loss)
}

# Function that returns the reduction in loss achieved by each split
get_variable_importance <- function(tree, loss=NA, variable_importance = list()) {
  # If the input is the base tree, extract the first node to start the recursion
  if (!is.null(tree[["base_loss"]])) {
    # If there is no node, i.e. the tree as a depth of 0, we return the base loss
    if (is.null(tree[["node"]])) {
      return(NA)
    }
    loss <- tree[["base_loss"]]
    tree <- tree[["node"]]
  }
  # Get the splitting variable name
  variable <- tree$variable
  # Compute loss after split:
  loss_new <- tree$loss_left + tree$loss_right
  # Compute the improvement in loss achieved by the splitting variable
  improvement <- loss - loss_new
  # If the splitting variable was already previously used, add the improvement
  if (!is.null(variable_importance[[variable]])) {
    improvement <- variable_importance[[variable]] + improvement
  }
  # Save the improvement
  variable_importance[[variable]] <- improvement
  
  # Handle the left child
  if (!is.null(tree[["left"]])) {
    variable_importance <- get_variable_importance(tree = tree[["left"]], 
                                                   loss = tree$loss_left, 
                                                   variable_importance = variable_importance)
  }
  # Handle the right child
  if (!is.null(tree[["right"]])) {
    variable_importance <- get_variable_importance(tree = tree[["right"]], 
                                                   loss = tree$loss_right, 
                                                   variable_importance = variable_importance)
  }
  
  return(variable_importance)
}

# Function that returns the number of internal nodes in the tree
get_number_internal_nodes <- function(tree, num_internal=1) {
  # If the input is the base tree, extract the first node to start the recursion
  if (!is.null(tree[["base_loss"]])) {
    # If there is no node, i.e. the tree as a depth of 0
    if (is.null(tree[["node"]])) {
      return(0)
    }
    tree <- tree[["node"]]
  }
  # Handle the left child
  if (!is.null(tree[["left"]])) {
    num_internal <- get_number_internal_nodes(tree[["left"]], num_internal+1)
  }
  # Handle the right child
  if (!is.null(tree[["right"]])) {
    num_internal <- get_number_internal_nodes(tree[["right"]], num_internal+1)
  }
  return(num_internal)
}

# Function that returns a sub-tree by folding all nodes from a specific array of node-ids
get_sub_tree <- function(tree, id_remove) {
  
  if (is.null(id_remove)) {
    return(tree)
  }
  
  if (!is.null(tree[["node"]])) {
    tree[["node"]] <- get_sub_tree(tree[["node"]], id_remove)
    return(tree)
  }
  
  
  if (tree$id %in% id_remove) {
    return(NULL)
  } 
  
  if (!is.null(tree[["left"]])) {
    tree[["left"]] <- get_sub_tree(tree[["left"]], id_remove)
  }
  
  if (!is.null(tree[["right"]])) {
    tree[["right"]] <- get_sub_tree(tree[["right"]], id_remove)
  }
  
  return(tree)
}

# Computes the complexity cost of the tree. Currently AIC or BIC
get_complexity_cost <- function(tree, ic = "aic") {
  sse <- get_tree_loss(tree)
  num_internal <- get_number_internal_nodes(tree)
  num_terminal <- num_internal + 1
  nobs <- tree$nobs
  nparam <- num_terminal*(length(coef(tree$base_model)) + tree$nids)
  k <- switch (ic,
               "aic" = 2,
               "bic" = log(nobs)
  )
  # Complexity cost:
  complexity_cost <- k*nparam + nobs*log(sse)
  return(complexity_cost)
}

# Prune the tree to find the optimal tre based on a cost function
prune_tree <- function(tree, ic = "aic") {
  # Determine number of internal nodes:
  num_internal <- tree$node_num
  # Determine all combinations of nodes to drop:
  nodes_fold <- c()
  for (i in 1:num_internal) {
    nodes_fold <- c(nodes_fold, combn(1:num_internal, i, simplify = FALSE))
  }
  # Add the possibility that no node has to be removed:
  nodes_fold <- c(nodes_fold, list(NULL))
  cost_criterion <- sapply(nodes_fold, 
                           FUN = function(ids_remove) get_complexity_cost(get_sub_tree(tree, id_remove = ids_remove), ic = ic))
  pruned_tree <- get_sub_tree(tree, nodes_fold[[which.min(cost_criterion)]])
  return(pruned_tree)
}

# Print the resulting tree
print_tree <- function(tree, shift="", var_names = NULL) {
  # If the imput is the base tree, extract the first node to start the recursion
  if (!is.null(tree[["base_loss"]])) {
    # If there is no node, i.e. the tree as a depth of 0
    if (is.null(tree[["node"]])) {
      cat(shift, "Root:\n", coefficients(tree$base_model))
      return()
    }
    var_names <- tree[["variables"]]
    tree <- tree[["node"]]
  }
  # Handle the left child
  if (!is.null(tree[["left"]])) {
    cat(shift, tree$variable, "<", tree$value, ":\n")
    print_tree(tree[["left"]], paste0(shift, "\t"), var_names)
  } else {
    cat(shift, tree$variable, "<", tree$value, ":", coefficients(tree$model_left), "\n")
  }
  # Handle the right child
  
  if (!is.null(tree[["right"]])) {
    cat(shift, tree$variable, ">=", tree$value, ":\n")
    print_tree(tree[["right"]], paste0(shift, "\t"), var_names)
  } else {
    cat(shift, tree$variable, ">=", tree$value, ":", coefficients(tree$model_right), "\n")
  }
}

# Get the logical expressions that lead to each of the leaf-models
get_path_to_leaves <- function(tree, path_to_leaves=list(), path=NULL) {
  # If the input is the base tree, extract the first node to start the recursion
  if (!is.null(tree[["base_loss"]])) {
    # If there is no node
    if (is.null(tree[["node"]])) {
      return("")
    }
    tree <- tree[["node"]]
  }
  # Handle the left child:
  new_path_left <- paste(tree[["variable"]], "<", tree[["value"]])
  new_path_left <- ifelse(is.null(path), new_path_left, paste(path, new_path_left, sep=" & "))
  if (!is.null(tree[["left"]])) {
    path_to_leaves <- get_path_to_leaves(tree[["left"]], path_to_leaves, new_path_left)
  } else {
    leaf_name <- paste(tree[["id"]], "left", sep="_")
    leaf <- list("path"=new_path_left, "model"=tree[["model_left"]])
    path_to_leaves[[leaf_name]] <- leaf
  }
  # Handle the right child:
  new_path_right <- paste(tree[["variable"]], ">=", tree[["value"]])
  new_path_right <- ifelse(is.null(path), new_path_right, paste(path, new_path_right, sep=" & "))
  if (!is.null(tree[["right"]])) {
    path_to_leaves <- get_path_to_leaves(tree[["right"]], path_to_leaves, new_path_right)
  } else {
    leaf_name <- paste(tree[["id"]], "right", sep="_")
    leaf <- list("path"=new_path_right, "model"=tree[["model_right"]])
    path_to_leaves[[leaf_name]] <- leaf
  }
  return(path_to_leaves)
}

# Make predictions with a tree model
predict_tree <- function(tree, newdata) {
  # Ensure that the new data is a data frame object
  newdata <- as.data.frame(newdata)
  # Add an identifier to the data: this helps in putting the predictions in the correct order
  newdata[,"__identifier__"] <- 1:nrow(newdata)
  # Obtain the simplified tree (logical paths to the leaves and the respective models)
  simple_tree <- get_path_to_leaves(tree)
  # Initialize the vector of predictions
  predictions <- vector(mode="numeric", length = nrow(newdata))
  # For each leaf, make the predictions
  for (leaf_i in simple_tree) {
    # Filter out all observations that are in the current leaf
    newdata_i <- subset(newdata, eval(parse(text=leaf_i$path)))
    # Make the predictions
    predictions[newdata_i[,"__identifier__"]] <- predict(leaf_i$model, newdata_i)
  }
  
  return(predictions)
}

