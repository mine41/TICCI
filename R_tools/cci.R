# Extract ligand and receptor
extract_gene <- function(input_string) {
  # Define the extraction function
  extract_part <- function(part) {
    part <- gsub("[()]", "", part)  # Remove parentheses
    part_split <- strsplit(part, "\\+")[[1]]
    return(part_split)
  }
  
  # Use regular expressions to extract A and B
  matches <- regmatches(input_string, gregexpr("[^-]+", input_string))[[1]]
  
  # Extract A
  A <- matches[1]
  current_ligand <- extract_part(A)
  
  # Extract B
  B <- matches[2]
  current_receptor <- extract_part(B)
  
  result <- list(current_ligand = current_ligand, current_receptor = current_receptor)
  return(result)
}


# Find CCI for each communication result
find_each_cci <- function(my_data, my_meta, current_source, current_target, current_ligand, current_receptor, prob) {
  # Filter sender cells
  sender_cells <- rownames(my_data)[my_meta$Classification == current_source]
  sender_cells <- sender_cells[rowSums(my_data[sender_cells, current_ligand, drop = FALSE] > 1) == length(current_ligand)]
  # print("Sender")
  # print(sender_cells)
  
  # Filter receiver cells
  receiver_cells <- rownames(my_data)[my_meta$Classification == current_target]
  receiver_cells <- receiver_cells[rowSums(my_data[receiver_cells, current_receptor,drop = FALSE] > 1) == length(current_receptor)]
  # print("Receiver")
  # print(receiver_cells)
  
  # Generate all combinations
  cci_list <- expand.grid(sender_cells = sender_cells, receiver_cells = receiver_cells)
  cci_list$ligand <- apply(my_data[cci_list$sender_cells, current_ligand, drop = FALSE], 1, function(x) paste(names(x), collapse = "_"))
  cci_list$receptor <- apply(my_data[cci_list$receiver_cells, current_receptor, drop = FALSE], 1, function(x) paste(names(x), collapse = "_"))
  if (nrow(cci_list) > 0) {
    cci_list$prob <- prob
    cci_list$ligand <- paste(current_ligand, collapse = "_")
    cci_list$receptor <- paste(current_receptor, collapse = "_")
  } else {
    # Handle the case of an empty data frame, you can choose to create a new data frame or perform other logic
    print("None Communication cell")
  }
  
  return(cci_list)
}

# Find all CCIs
get_all_cci <- function(my_data, my_meta, df.net) {
  col_names <- c("sender_cells", "receiver_cells", "ligand", "receptor", "prob")
  CCIList_all <- data.frame(matrix(ncol = length(col_names), nrow = 0))
  colnames(CCIList_all) <- col_names
  gc()
  for (i in 1:nrow(df.net)) {
    # for (i in 1:10) {
    current_row <- df.net[i, ]
    extract_current_interaction_name_2 <- extract_gene(gsub(" ", "", df.net[i,'interaction_name_2']))
    current_source <- current_row['source'][[1]] # Source cluster
    current_target <- current_row['target'][[1]] # Target cluster
    current_ligand <- extract_current_interaction_name_2$current_ligand # Ligand genes
    current_receptor <- extract_current_interaction_name_2$current_receptor # Receptor genes
    current_prob <- current_row['prob'][[1]] # Communication strength
    print(i)
    print(df.net[i,'interaction_name_2'])
    # print(current_source, current_target, current_prob)
    
    # Call the function
    temp_cci <- find_each_cci(my_data, my_meta, current_source, current_target, current_ligand, current_receptor, current_prob)
    CCIList_all <- rbind(CCIList_all, temp_cci)
  }
  return(CCIList_all)
}
