###
#
# This is all code modified/extended from the stmCorrViz package (https://github.com/cran/stmCorrViz/tree/master/R).
# Thank you Antonio Coppola, Margaret Roberts, Brandon Stewart, and Dustin Tingley.
#
###

# requires questions, question_names
find_cluster_docs <- function(topic_nums, raw_docs, question_num, min_total_prev = 0.3) {
  responses_theta <- data.frame(raw_docs, thetas[[question_names[question_num]]][-1])
  names(responses_theta) <- c('responses', names(responses_theta)[2:(ncol(responses_theta))])

  min_prevalence <- min_total_prev/length(topic_nums)
  responses_theta <- data.frame(responses_theta[1], responses_theta[-1] %>% select(topic_nums))

  filter_string <- ''
  for(i in seq(ncol(responses_theta)-1)) {
    filter_string <- paste(filter_string, names(responses_theta)[i+1], '>', min_prevalence, '&')
  }
  filter_string <- filter_string %>% substr(2, nchar(filter_string)-2)

  responses_sum <- data.frame(responses_theta,
                              variance = responses_theta[-1] %>% apply(1, function(x)(diff(range(x)))) %>% as.vector,
                              sum = responses_theta[-1] %>% apply(1, function(x)(sum(x))) %>% as.vector) %>%
    filter_(filter_string) %>%
    filter(sum >= max(min(length(topic_nums)/10, 0.9), 0.2)) %>%
    filter(variance < 0.1) %>%
    arrange(variance) %>%
    select(1, sum, variance) %>%
    slice(1:50)

  responses_sum$responses %<>% as.character
  return(responses_sum)
}
cluster_travel <- function(datajs, raw_docs, question_num) {
  if(!is.null(datajs$children)) {
    docs <- find_cluster_docs(as.numeric(datajs$topic_no), raw_docs, question_num)
    datajs$thoughts <- docs$responses
    datajs$thought_proportions <- docs$sum
    datajs$thought_variances <- docs$variance
    for(i in seq(datajs$children)) {
      datajs$children[[i]] <- cluster_travel(datajs$children[[i]], raw_docs, question_num)
    }
  }
  return(datajs)
}
create_json <- function(stm, documents_raw, documents_matrix, column_name,
                        title='STM Tree', clustering_thresh=Inf, labels_number=7,
                        verbose=F, instant=F, topic_labels=NULL, cluster_labels=NULL, directory=NULL,
                        question_num=NULL)
{
  # with instant = True, the file is written to data.js
  # with instant = False, the file is named <question>_<#ofTopics>_data.js
  json <- stmJSON(mod = stm, documents_raw = documents_raw,
                  documents_matrix = documents_matrix,
                  topic_labels = topic_labels,
                  cluster_labels = cluster_labels,
                  title = title,
                  clustering_threshold = clustering_thresh,
                  labels_number = labels_number,
                  verbose = verbose) %>%
    .$json

  name = paste(column_name, ncol(stm$theta),'data.js', sep='_')
  if(instant){ name = 'data.js'}
  if(!is.null(directory)){
    if(substr(directory, nchar(directory), nchar(directory)) != '/') {
      directory <- paste0(directory, '/')
    }
    name <- paste0(directory, name)
  }
  write(json, name, sep='')
  if(!is.null(question_num)) {
      json <- cluster_travel(read_json(name, simplifyVector=T, simplifyDataFrame = F,
                                   simplifyMatrix = F), documents_raw, question_num)
      json <- toJSON(json)
  }
  json <- paste0('var stm_data = ', json)
  write(json, name, sep='')
}

stmJSON <-
  function(mod, documents_raw=NULL, documents_matrix=NULL,
           topic_labels=NULL, cluster_labels=NULL,
           title="STM Model",clustering_threshold=1.5,
           labels_number=3, verbose=T){

    # Generate baseline topic list
    out <- list()

    if(verbose==TRUE)
      cat("Performing hierarchical topic clustering ... \n")

    # Run hclust subroutine
    clust <- clusterAnalysis(mod, labels_number, topic_labels)
    if(clustering_threshold == Inf) {
      clustering_threshold <- max(clust$height) - 0.1
    }

    if(verbose==TRUE)
      cat("Generating JSON representation of the model ... \n")

    # Extra data objects
    thoughts <- stm::findThoughts(mod, documents_raw, n=50)
    full_labels <- stm::labelTopics(mod)
    topic_proportions <- colMeans(mod$theta)

    # Find aggregation points
    K <- mod$settings$dim$K
    topic_to_topic_splits <- c()
    for(i in seq(K-1))
      if(clust$merge[i,1] <0 && clust$merge[i,2] < 0)
        topic_to_topic_splits <- c(topic_to_topic_splits, i)
    aggregate <- setdiff(which(clust$height <= clustering_threshold),
                         topic_to_topic_splits)

    # Produce merge list
    merge_list <- list()
    for(i in seq(K-1))
      merge_list[[i]] <- clust$merge[i,]
    names(merge_list) <- 1:(K-1)

    # Collapse merge list
    merge_list <- collapseMergeList(merge_list, clust$merge, aggregate, K)

    # Build collapsed-clusters data structure
    top_layer <- merge_list[paste(K-1)][[1]]
    out$children <- buildClusters(list(), current = top_layer, merge_list,
                                  labels=clust$labels, full_labels,
                                  thoughts, topic_proportions, mod$theta)

    # Implementation with diagnostics
    #   out$children <- buildClusters(list(), current = top_layer, merge_list,
    #                                 labels=clust$labels, full_labels,
    #                                 thoughts, exclusivity_scores,
    #                                 semcoh_scores)

    # Get beta weights for model
    beta_weights <- getBetaWeights(mod, documents_matrix)

    # Assign cluster names
    if(!is.null(cluster_labels)) {
      sequence <- sapply(cluster_labels, function(x)(x$clustNum)) %>% as.vector
      temp <- c()
      for(i in seq(sequence)) {
        temp <- c(temp, list(cluster_labels[[sequence[i]]]))
      }
      cluster_labels <- temp
    }
    out <- assignClusterNames(out, labels_number, beta_weights, mod$vocab, cluster_labels)

    # Root Information
    out$name <- title
    out$this_root <- TRUE
    out$summary <- utils::capture.output(mod)
    out$proportions <- topic_proportions

    # Convert structure to JSON
    out_JSON <- jsonlite::toJSON(out, force=TRUE)
    return(list(json=out_JSON, n_merge=length(merge_list)))
  }

clusterAnalysis <- function (stmobj, labels_number = 3, topic_labels=NULL)
{
  labels <- stm::labelTopics(stmobj)
  theta <- stmobj$theta
  d <- stats::dist(stats::cor(theta))
  clust <- stats::hclust(d, method = "complete")
  K <- stmobj$settings$dim$K
  clust$labels <- rep(NA, K)
  for (i in seq(K)) {
    if(length(topic_labels[[i]]) > 0) {
      clust$labels[i] <- topic_labels[[i]]$name
    }
    else {
      l <- labels$frex[i, ]
      l <- paste(l[1:labels_number], collapse = ", ")
      clust$labels[i] <- l
    }
  }
  return(clust)
}

collapseMergeList <-
  function(merge_list, merge_matrix, aggregate, K){
    calls <- which(merge_matrix %in% aggregate)

    # Generate deletion sequence
    old_refs <- c()
    row_indices <- c()
    col_indices <- c()
    for(i in seq(length(calls))){
      call <- calls[i]
      if(call <= K-1){
        row_indices[i] <- call
        col_indices[i] <- 1
      } else {
        row_indices[i] <- call - K + 1
        col_indices[i] <- 2
      }
      old_refs[i] <- merge_list[row_indices[i]][[1]][col_indices[i]]
    }
    deletion_seq <- data.frame(old_refs=old_refs,
                               row_indices=row_indices,
                               col_indices=col_indices)
    deletion_seq <- deletion_seq[order(old_refs),]

    # Perform collapsing: Insert new sequences
    for(i in seq(nrow(deletion_seq))){
      row_index <- deletion_seq$row_indices[i]
      col_index <- deletion_seq$col_indices[i]
      old_ref <- deletion_seq$old_refs[i]
      #merge_list[row_index][[1]] <- merge_list[row_index][[1]][-col_index]
      merge_list[row_index][[1]] <- c(merge_list[row_index][[1]],
                                      merge_list[old_ref][[1]])
    }

    # Perform collapsing: Delete old references
    for(row_index in names(merge_list)){
      delete <- which(merge_list[paste(row_index)][[1]] %in% aggregate)
      if(any(delete))
        merge_list[row_index][[1]] <- merge_list[row_index][[1]][-delete]
    }
    merge_list <- merge_list[-aggregate]
    return(merge_list)
  }

buildClusters <-
  function(out, current, merge_list, labels, full_labels,
           thoughts, topic_proportions, theta){

    # Recursive definition
    for(i in seq(length(current))){
      out[[i]] <- list()
      out[[i]]$name <- current[i]
      if(current[i] > 0){
        out[[i]]$children <- buildClusters(list(),
                                           merge_list[paste(current[i])][[1]],
                                           merge_list, labels=labels, full_labels,
                                           thoughts, topic_proportions, theta)
        out[[i]]$name <- current[i]
      } else {
        out[[i]]$size <- 1800 # if removed, code silently fails
        out[[i]]$name <- labels[-current[i]]
        out[[i]]$topic_no <- -current[i]
        out[[i]]$thoughts <- c()
        out[[i]]$thought_proportions <- c()
        for(j in seq(50)) {
          if(0.2 <= theta[thoughts$index[[-current[i]]][j], -current[i]]) {
            out[[i]]$thoughts <- c(out[[i]]$thoughts, iconv(thoughts$docs[[-current[i]]][j], to='utf-8', sub=""))
            out[[i]]$thought_proportions <- c(out[[i]]$thought_proportions, theta[thoughts$index[[-current[i]]][j], -current[i]])
          }
          else {
            # minimum of 10 documents
            if(j >= 10) break
            for(k in j:10) {
              out[[i]]$thoughts <- c(out[[i]]$thoughts, iconv(thoughts$docs[[-current[i]]][k], to='utf-8', sub=""))
              out[[i]]$thought_proportions <- c(out[[i]]$thought_proportions, theta[thoughts$index[[-current[i]]][k], -current[i]])
            }
            break
          }
        }
        out[[i]]$prob <- paste(full_labels$prob[-current[i],], collapse = ", ")
        out[[i]]$frex <- paste(full_labels$frex[-current[i],], collapse = ", ")
        out[[i]]$lift <- paste(full_labels$lift[-current[i],], collapse = ", ")
        out[[i]]$score <- paste(full_labels$score[-current[i],], collapse = ", ")
        out[[i]]$proportion <- format(round(topic_proportions[-current[i]], 2))
      }
    }
    return(out)
  }

assignClusterNames <-
  function(out, lab_no, beta_weights, vocab, cluster_labels=NULL){

    members <- c()

    for(i in seq(length(out$children))){
      if (!('size' %in% names(out$children[[i]]))) # if child i is cluster
        out$children[[i]] <- assignClusterNames(out$children[[i]], lab_no,
                                                beta_weights, vocab, cluster_labels)
    }

    for(i in seq(length(out$children))){
      members <- c(members, out$children[[i]]$topic_no)
    }

    out$topic_no <- members
    margins <- marginalize(members, beta_weights$beta, beta_weights$weights)
    labels <- vocab[margins$indices[1:lab_no]]

    out$name <- paste(sample(labels, lab_no), collapse=", ")
    if(!is.null(cluster_labels)) {
      for (j in seq(cluster_labels)) {
        if(sum(members %in% cluster_labels[[j]]$topics) == length(members %in% cluster_labels[[j]]$topics) &
           length(members) == length(cluster_labels[[j]]$topics)) {
          if(!is.null(cluster_labels[[j]]$name[1]))
            out$name <- cluster_labels[[j]]$name[1]
          break
        }
      }
    }

    return(out)
  }

getBetaWeights <-
  function(model, documents=NULL) {
    logbeta <- model$beta$logbeta
    K <- model$settings$dim$K
    vocab <- model$vocab

    #Let's start by marginalizing
    margbeta <- exp(logbeta[[1]])
    if(length(logbeta) > 1) {
      weights <- model$settings$covariates$betaindex
      tab <- table(weights)
      weights <- tab/sum(tab)
      #marginalize
      margbeta <- margbeta*weights[1]
      for(i in 2:length(model$beta$logbeta)) {
        margbeta <- margbeta + exp(model$beta$logbeta[[i]])*weights[i]
      }
    }

    ##
    # figure out how to weight the topics.
    # NB: if they didn't provide topics use naive weights
    #     otherwise calibrate thetas by the total counts
    #     per document.
    if(is.null(documents)) {
      weights <- colSums(model$theta)
    } else {
      D.n <- unlist(lapply(documents, function(x) sum(x[2,])))
      weights <- colSums(D.n*model$theta)
    }

    return(list(beta=margbeta, weights=weights))
  }

marginalize <-
  function(members, beta, weights) {
    w <- weights[members]/sum(weights[members])
    pvec <- colSums(beta[members,,drop=FALSE]*w)
    words <- order(pvec, decreasing=TRUE)
    return(list(beta=pvec, indices=words))
  }
