
# Function to calculate Jaccard similarity between gene sets
calculate_jaccard <- function(set1, set2) {
  intersection_size <- length(intersect(set1, set2))
  union_size <- length(union(set1, set2))
  return(intersection_size / union_size)
}

# Function to filter redundant pathways based on gene overlap and NES scores
filter_redundant_pathways <- function(results_df, pathways, similarity_threshold = 0.5) {
  # Convert results to data frame if not already
  results_df <- as.data.frame(results_df)
  
  # Sort by absolute NES score
  results_df <- results_df[order(abs(results_df$NES), decreasing = TRUE), ]
  
  # Initialize vector to store indices of pathways to keep
  keep_pathways <- c(1)  # Always keep the pathway with highest absolute NES
  
  # Compare each pathway with already selected ones
  for(i in 2:nrow(results_df)) {
    current_pathway <- results_df$pathway[i]
    current_genes <- pathways[[current_pathway]]
    
    # Check overlap with all kept pathways
    is_redundant <- FALSE
    for(j in keep_pathways) {
      compare_pathway <- results_df$pathway[j]
      compare_genes <- pathways[[compare_pathway]]
      
      # Calculate Jaccard similarity
      similarity <- calculate_jaccard(current_genes, compare_genes)
      
      # If similarity is above threshold and pathways have same direction of enrichment
      if(similarity > similarity_threshold && 
         sign(results_df$NES[i]) == sign(results_df$NES[j])) {
        is_redundant <- TRUE
        break
      }
    }
    
    # If not redundant, add to keep list
    if(!is_redundant) {
      keep_pathways <- c(keep_pathways, i)
    }
  }
  
  # Return filtered results
  return(results_df[keep_pathways, ])
}
