#' Apply function that can safely accept a vector
#' 
#' @param arr An array of arbitrary dimensions. Can be a vector.
#' @param mar A vector giving the subscripts which the function will be applied 
#' over. Ignored when `arr` is a vector. 
#' @param FUN The function to be applied.
#' @param ... Optional arguments to FUN.
#' @return 
#' When `arr` is a vector, returns `FUN(arr, ...)`. When `arr` is an array with 
#' 2 or more dimensions, returns `apply(arr, mar, FUN, ...)`. 
safe_apply <- function(arr, mar, FUN, ...) {
  if (is.null(dim(arr))) {
    return(FUN(arr, ...))
  }
  else {
    return(apply(arr, mar, FUN, ...))
  }
}


#' Replace values with another value (`NA` by default)
#' 
#' @details 
#' The expected usage is in the context of an apply function. 
#' @param na_val Value to replace. 
#' @param rep Replacement value. The default is `NA`. 
#' @return A function that takes a single parameter, `x`. `x` should be 
#' specified as a vector that is compatible with `na_val` and `rep`. Evaluating 
#' the return function will return `x` with values of `na_val` replaced by 
#' `rep`. 
fix_with_rep <- function(na_val, rep = NA) {
  function(x) {
    x[x %in% na_val] <- rep
    x[x < 0] <- rep
    x
  }
}


#' Extract the cancers from the column names of a pedigree data frame
#' 
#' @param fam A pedigree data frame. 
#' @return Character vector of short (abbreviated) names of cancers extracted 
#' directly from the `isAFf*` columns of a pedigree. Does not ensure that a 
#' corresponding `Age*` column exists. 
.getCancersFromFam <- function(fam) {
  cn <- dimnames(fam)[[2]]
  cancers <- regmatches(cn, regexpr("(?<=^isAff).*", cn, perl = T))
  cancers <- cancers[!cancers %in% "Any"]
  return(cancers)
}


#' Subset array by a dimension of choice and optionally reduce that dimension
#' 
#' @param arrayInput A multi-dimensional array. 
#' @param axisName A character string giving one or more named axis dimensions 
#' of `arrayInput` of which to index. 
#' @param choice A list of indices to subset for each dimension. 
#' @param ... Additional arguments passed in to `asub`. 
#' @return 
subset_array <- function(arrayInput, axisName, choice, ...) {
  dimNames <- dimnames(arrayInput)
  axisId <- which(names(dimnames(arrayInput)) %in% axisName)
  abind::asub(arrayInput, choice, axisId, ...)
}


#' Get possible genotypes
#' 
#' Generates all possible genotypes based on a list of gene variants, 
#' restricted to a maximum number of mutations. 
#' @param gene A character vector of gene names
#' @param max_mut The maximum number of simultaneous mutations allowed. 
#' @param homo_genes A list of genes with heterozygous and homozygous mutations. 
#' Each component is named after a gene and contains a character vector of the 
#' full gene names for that gene. The default is an empty list `list()`, which 
#' indicates that no genes have both heterozygous and homozygous mutations. 
#' @return A list with two components: 
#' - `df`: A data frame of possible genotypes for the genes in `gene`, 
#' restricted to at most `max_mut` simultaneous mutations. Each row represents 
#' a genotype and each column represents a gene variant. The data frame is 
#' comprised of `0` and `1` values that encode the genotypes. 
#' - `list`: A character vector of genotype names corresponding to the rows of 
#' `df`. 
#' @seealso \code{\link{.ppPossibleGenotypes}}
.getPossibleGenotype <- function(gene, max_mut = 2, 
                                 homo_genes = list()) {
  
  # Indices for gene with heterozygous and homozygous mutations
  if (length(homo_genes) == 0) {
    # Empty list if there are no genes with heterozygous and homozygous 
    # mutations
    homo_idx <- list()
  } else {
    # Otherwise, identify the index pairs for heterozygous and homozygous 
    # mutations
    homo_idx <- lapply(homo_genes, function(g) {
      which(gene %in% g)
    })
  }
  
  # Get data frame of possible genotypes
  genodf <- .ppPossibleGenotypes(length(gene), max_mut, 
                                 homo_idx, collapse = FALSE)
  colnames(genodf) <- gene
  
  # Generate corresponding names for possible genotypes
  gps <- apply(genodf, 1, 
               function(v) paste0(colnames(genodf)[as.logical(v)], 
                                  collapse = "."))
  gps[1] <- "noncarrier"
  
  return(list(df = as.data.frame(genodf), list = gps))
}


#' Drop array dimensions that have only one level
#' 
#' @param arr An array. 
#' @return 
#' An array or vector with reduced dimensions from dropping dimensions with 
#' only one level. 
adrop_one <- function(arr) {
  abind::adrop(arr, drop = which(dim(arr) == 1))
}


#' Ensure that pedigree data frame column only contains allowed values
#' 
#' Used to check categorical and logical columns in \code{\link{checkFam}}. 
#' @param ped A pedigree data frame. 
#' @param col_spec The name of a column if `ped` with a specific attribute. 
#' @param allowed A vector of the values that are allowed in `ped[[col_spec]]`. 
#' @param default_NA_rep The value to default to if `NA`s are encountered. The 
#' default is `NULL`. 
#' @param NA_Ok A logical value indicating whether `NA`s can be present in 
#' `ped[[col_spec]]`. The default is `TRUE`. 
#' @param vec A vector. The default is `NULL`. If `vec` is not `NULL`, the 
#' function will check `vec` instead of `ped[[col_spec]]`. 
#' @param `force_NA` A logical value indicating whether non-permissible should 
#' be replaced with `NA`s. The default is `TRUE`. 
#' @param `force_remove` A logical value indicating whether non-permissible 
#' values should be removed. The default is `FALSE`. 
#' @return The vector `ped[[col_spec]]` (or `vec` if `vec` is not `NULL`) with 
#' non-permissible values either removed and replaced. 
forcingNA_ifnot_contains <- function(ped, col_spec, allowed,
                                     default_NA_rep = NULL,
                                     NA_Ok = TRUE, vec = NULL,
                                     force_NA = TRUE, force_remove = FALSE) {
  if (is.null(vec)) {
    vec <- ped[[col_spec]]
  }
  contains_check <- .ifContains(vec, allowed, NA_Ok)
  if (!is.null(contains_check)) {
    if (force_NA & !force_remove) {
      if (sum(vec %in% contains_check) > 0) {
        msg <- sprintf(
          "Column %s contains %s that is not recognized. Only %s are valid inputs. Forcing default/NA ..",
          col_spec,
          paste0(unique(contains_check), collapse = ","),
          paste0(allowed, collapse = ",")
        )
        rlang::warn(msg, level = "InvalidInputForcingNA")
        vec[vec %in% contains_check] <- NA
      }
    }
    if (force_remove & !force_NA) {
      if (sum(!vec %in% contains_check) > 0) {
        msg <- sprintf(
          "Column %s contains %s that is not recognized. Only %s are valid inputs. Forcing removal ..",
          col_spec,
          paste0(unique(contains_check), collapse = ","),
          paste0(allowed, collapse = ",")
        )
        rlang::warn(msg, level = "InvalidInputForcingNA")
        vec <- vec[!vec %in% contains_check]
      }
    }
  }
  if (!is.null(default_NA_rep)) {
    if (sum(is.na(vec)) > 0) {
      msg <- sprintf(
        "Column %s does not support NA. Forcing to %s ..",
        col_spec,
        default_NA_rep
      )
      rlang::warn(msg, level = "InvalidNA")
      vec[is.na(vec)] <- default_NA_rep
    }
  }
  return(vec)
}


#' Base function used to map values using a dictionary
#' 
#' @param dict A list with two named components that map to each other. The 
#' default is `PanelPRO:::CANCER_NAME_MAP`, which maps short and long cancer 
#' names. 
#' @return A function where the input should be a keyword argument 
#' corresponding to one of the named components of `dict`. This will be treated 
#' as the "from" component. The return value will be the input mapped to the 
#' values of the other named component of `dict` (the "to" component).
#' @seealso \code{\link{.mapCancerNames}}, \code{\link{.mapGenderNames}} 
.mapDict <- function(dict = CANCER_NAME_MAP) {
  function(...) {
    args <- list(...)
    if (any(!names(args) %in% names(dict))) {
      rlang::abort("Some arguments cannot be found in dictionary",
                   level = "MatchError"
      )
    }
    if (length(args) > 1) {
      rlang::abort("Argument list has a length larger than one",
                   level = "MatchError"
      )
    }
    input_key <- names(args)[1]
    key_to_match <- names(dict)[!names(dict) %in% input_key]
    matched_indices <- match(args[[1]], dict[[input_key]])
    
    if (any(is.na(matched_indices))) {
      rlang::warn("Returned indices contains NA, please double check!",
                  level = "MatchWarning"
      )
    }
    return(dict[[key_to_match]][matched_indices])
  }
}


#' Convert cancer short names to long names or vice versa
#' 
#' @param ... A keyword argument that represents the "from" component of the 
#' dictionary, one of: 
#' * `short`: A character vector of short cancer names that should be mapped to 
#' long cancer names. 
#' * `long`: A character vector of long cancer names that should be mapped to  
#' short cancer names. 
#' @return A vector of the input values mapped to the "to" component of the 
#' dictionary (either `long` or `short`). 
#' @seealso \code{\link{.mapDict}}
.mapCancerNames <- .mapDict(CANCER_NAME_MAP)


#' Convert between numeric-valued sex and character-valued sex
#' 
#' @param ... A keyword argument that represents the "from" component of the 
#' dictionary, one of: 
#' * `number`: A numeric vector with values that are a subset of 
#' `c(0, 1, NA, -999)` and will be mapped to 
#' `c("Female", "Male", "All_Sexes", "All_Sexes")`. 
#' * `character`: A character vector with values that are a subset of 
#' `c("Female", "Male", "All_Sexes", "All_Sexes")` and will be mapped to 
#' `c(0, 1, NA, -999)`. 
#' @return A vector of the input values mapped to the "to" component of the 
#' dictionary (either `character` or `number`). 
#' @seealso \code{\link{.mapDict}}
.mapGenderNames <- .mapDict(list(
  number = c(0, 1, NA, -999),
  character = c("Female", "Male", rep("All_Sexes", 2))
))


#' Identifies values that are not allowed in a vector
#' 
#' @param vec A vector. 
#' @param allowed_vals A vector of the values that are allowed in `vec`. The 
#' default is `c(0, 1)`. 
#' @param NA_ok A logical value indicating whether `NA`s can be present in 
#' `vec`. The default is `TRUE`. 
#' @return Vector of values in `vec` that are not in `allowed_vals` (or `NA`, 
#' if `NA_ok = TRUE`). If no values were found, returns `NULL`. 
.ifContains <- function(vec, allowed_vals = c(0, 1),
                        NA_ok = TRUE) {
  if (NA_ok) {
    vec <- vec[!is.na(vec)]
  }
  ifin <- vec %in% allowed_vals
  if (!all(ifin)) {
    who_not_in <- vec[!ifin]
    if (length(who_not_in) != 0) {
      return(who_not_in)
    }
  }
  else {
    return(NULL)
  }
}


#' Evaluate function while ignoring `NA` and another predefined value
#' 
#' @param FUN The function to be applied. 
#' @param vec Vector over which to apply FUN. 
#' @param na_code Value to ignore in addition to `NA`. The default is `-999`. 
#' @return The output of `FUN(vec)`, while ignoring values in `c(NA, na_code)`. 
.safeGet <- function(FUN, vec, na_code = -999) {
  if (class(vec) == "list") {
    tryCatch(vec <- unlist(vec),
             error = function(c) "Input is not a unlistable list!"
    )
  }
  # If all components are NA
  if (all(is.na(vec))) {
    return(NA)
  }
  
  # If input is empty
  if (is.null(vec)) {
    return(NA)
  }
  
  # If there is NA
  vec <- vec[!is.na(vec)]
  
  # If there is matching NA code
  if (!is.null(na_code)) vec <- vec[vec != na_code]
  
  # If then there is nothing left
  if (length(vec) == 0) {
    return(NA)
  } else {
    return(FUN(vec))
  }
}


#' Force `NA` to be `FALSE`
#' 
#' @param vec A logical vector. 
#' @return `vec` with values of `NA` set to `FALSE`. 
.forceFalse <- function(vec) {
  vec[is.na(vec)] <- FALSE
  return(vec)
}


#' Get the random state
#' 
#' Uses `get0` to return `NULL` when `.Random.seed` doesn't exist, with 
#' `inherits = FALSE` to only search the global environment and not its parent. 
#' @return `.Random.seed` object. 
get_rand_state <- function() {
  get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
}


#' Set a random state
#' 
#' @details 
#' Assigning `state = NULL` might lead to unwanted consequences. 
#' @param state Current state. 
#' @return None. 
set_rand_state <- function(state) {
  if (!is.null(state)) {
    assign(".Random.seed", state, envir = .GlobalEnv, inherits = FALSE)
  }
}
