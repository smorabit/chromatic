

InterpretChromatic <- function(
    seurat_obj, # doesn't do anything right now
    chromHMM_states,
    mat,
    target,
    xgb_params,
    state_col = "name",
    nrounds = 1000,
    early_stopping_rounds = 25,
    verbose = 0
){

    print("1. Setting up the data...")

    # TODO: Checks?
    states_keep <- unique(mcols(chromHMM_states)[,state_col])
    states_keep <- states_keep[states_keep %in% colnames(mat)]
    X <- as.matrix(mat[,states_keep])

    # update later if needed?
    y_true <- target

    #---------------------------------------------------------------
    # XGBoost modeling
    #---------------------------------------------------------------
    
    print("2. XGBoost modeling...")

    # setup the XGBoost data matrix
    dtrain <- xgb.DMatrix(data = X, label = y_true)

    n_cells <- nrow(mat)

    # train with early stopping using a small validation split
    # TODO: add parameters here for the train/test split?
    train_idx <- sample(seq_len(n_cells), size = floor(0.8 * n_cells))
    dtr <- xgb.DMatrix(data = X[train_idx, ], label = y_true[train_idx])
    dval <- xgb.DMatrix(data = X[-train_idx, ], label = y_true[-train_idx])

    watchlist <- list(train = dtr, val = dval)

    # run the xgb model on the train split
    # TODO: add parameters here?
    model <- xgb.train(
        params = xgb_params,
        data = dtr,
        nrounds = nrounds,
        watchlist = watchlist,
        early_stopping_rounds = early_stopping_rounds,
        verbose = verbose
    )

    #---------------------------------------------------------------
    # SHAP values 
    #---------------------------------------------------------------

    print("3. Calculating SHAP values...")

    # get the SHAP values per cell and per feature
    shap_contribs <- predict(model, X, predcontrib = TRUE) 
    colnames(shap_contribs) <- c(colnames(X), "BIAS")

    # # 2) aggregate importance: mean absolute SHAP per feature
    # mean_abs_shap <- apply(abs(shap_contribs[, colnames(X), drop = FALSE]), 2, mean)
    # mean_abs_shap <- sort(mean_abs_shap, decreasing = TRUE)
    # print(mean_abs_shap)

    #---------------------------------------------------------------
    # setup the output (probably will change later)
    #---------------------------------------------------------------

    output <- list()
    output$model <- model 
    output$SHAP <- shap_contribs
    output$X <- X

    return(output)

}