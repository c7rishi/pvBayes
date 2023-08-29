#' Basis function for computing pFDR
#' @importFrom data.table :=
#' @importFrom data.table %like%
#' @noRd
pFDR0 <- function(lambda_draws,
                  test_stat,
                  optim = FALSE,
                  alpha = NULL,
                  k = 1.0){

  # test_stat is an user-defined function to calculate the
  # test statistic (critical value).

  tmp_draws <- posterior::as_draws_df(lambda_draws) %>%
    data.table::setDT() %>%
    data.table::melt(
      id.vars = c(".chain", ".iteration", ".draw"),
      variable.name = "parameter",
      value.name = "value"
    )
  # reshape the lambda_draws into a long format data.table
  # identified by chain, iteration, and draw,
  # (iteration and draw are the same actually)
  # and keep the joint name (AE, drug) as parameter.

  tmp_draws[
    ,
    obs_test_stat := test_stat(value),
    keyby = parameter
  ]
  # add new col named obs_test_stat
  # computed by the test_stat function, referenced by 'value'


  test_stat_mat <- tmp_draws[
    ,
    .(parameter, obs_test_stat)
  ] %>% # select columns
    unique() %>% # delete duplicate rows by AE-drug names
    {.[
      ,
      nm_AE_Drug := parameter %>% # extract AE and Drug into a vector of chr
        gsub("x\\[", "", x = .) %>%
        gsub("\\]", "", x = .) %>%
        strsplit("\\,")
    ][,
      `:=`( # the first element is AE and the second one is Drug
        AE = vapply(nm_AE_Drug, function(x) x[1], "name"),
        # %>% factor(levels = row.names()),
        Drug = vapply(nm_AE_Drug, function(x) x[2], "name")
      ) # add two new columns AE and Drug
    ][,
      .(AE, Drug, obs_test_stat)
    ]} %>%
    data.table::dcast.data.table( # reshape to wide format
      AE ~ Drug,
      value.var = "obs_test_stat"
    ) %>% # the rownames and colnames were disordered, fixed by adding the following two pipe
    {.[
     order(match(AE,rownames(lambda_draws)))
    ] } %>%
    data.table::setcolorder(.,c("AE", colnames(lambda_draws))) %>%
    {
      rn <- .$AE
      tmp_dt <- .
      tmp_dt[, AE := NULL] %>%
        data.matrix() %>% # convert to matrix
        magrittr::set_rownames(rn) %>%
        as.matrix()
    }

  dscov_mat <- tmp_draws[
    ,
    .(parameter, obs_test_stat)
  ] %>% # select columns
    unique() %>%
    .[, dscov := as.numeric(obs_test_stat > k)] %>% # greater than k is a dscov
    .[, obs_test_stat := NULL] %>% # remove obs_test_stat
    {.[
      ,
      nm_AE_Drug := parameter %>% # extract AE and Drug into a vector of Chr
        gsub("x\\[", "", x = .) %>%
        gsub("\\]", "", x = .) %>%
        strsplit("\\,")
    ][,
      `:=`( # the first element is AE and the second one is Drug
        AE = vapply(nm_AE_Drug, function(x) x[1], "name"),
        Drug = vapply(nm_AE_Drug, function(x) x[2], "name")
      ) # add two new columns AE and Drug
    ][,
      .(AE, Drug, dscov) # table with dscov 0 or 1
    ]} %>%
    data.table::dcast.data.table( # reshape to wide format
      AE ~ Drug,
      value.var = "dscov",
      fill = 0 # fill NA with 0
    ) %>% # the rownames and colnames were disordered, fixed by adding the following two pipe
    {.[
      order(match(AE,rownames(lambda_draws)))
    ] } %>%
    data.table::setcolorder(.,c("AE", colnames(lambda_draws)))%>%
    {
      rn <- .$AE
      tmp_dt <- .
      tmp_dt[, AE := NULL] %>%
        data.matrix() %>% # convert to matrix
        magrittr::set_rownames(rn) %>%
        as.matrix()
    }
  # generate the dscov matrix with names

  range_test_stat <- range(tmp_draws$obs_test_stat)

  tmp_prob_null_mat_full <- matrix(
    NA,
    nrow(lambda_draws),
    ncol(lambda_draws),
    dimnames = dimnames(lambda_draws)
  ) # empty table for store the probabilities, 0 means non-rejected cell

  tmp_prob_null_reject_sub_dt <- tmp_draws[
    obs_test_stat >= k, # for those rejected cell
    .(prob_null = mean(value <= 1.0)), # calculate the Bayes TIE
    keyby = parameter
  ] # long table, but not complete

  tmp_prob_null_mat_sub <- tmp_prob_null_reject_sub_dt %>%
    {
      if (nrow(.) >= 1) {
        tmp_in <- .
        tmp_in[
          ,
          nm_AE_Drug := parameter %>%
            gsub("x\\[", "", x = .) %>%
            gsub("\\]", "", x = .) %>%
            strsplit("\\,") # vector of AE and Drug names
        ][,
          `:=`(
            AE = vapply(nm_AE_Drug, function(x) x[1], "name"),
            Drug = vapply(nm_AE_Drug, function(x) x[2], "name")
          )
        ][,
          .(AE, Drug, prob_null)
        ] %>%
          data.table::dcast.data.table(
            AE ~ Drug,
            value.var = "prob_null",
            fill = NA
          ) %>%
          {
            rn <- .$AE
            tmp_dt <- .
            tmp_dt[, AE := NULL] %>%
              data.matrix() %>%
              magrittr::set_rownames(rn) %>%
              as.matrix()
          }
      } else {
        matrix(NA, nrow = 0, ncol = 0)
      }
    }
  # a matrix with Bayes TIE, those non-rejected cells
  # are filled with 0

  tmp_prob_null_mat_full[
    rownames(tmp_prob_null_mat_sub),
    colnames(tmp_prob_null_mat_sub)
  ] <- tmp_prob_null_mat_sub
  # table for output

  pFDR <- ifelse(
    nrow(tmp_prob_null_reject_sub_dt) >= 1,
    mean(tmp_prob_null_reject_sub_dt$prob_null,
         na.rm = TRUE),
    NA
  )

  # lambda_s1 <- lambda_draws  %>%
  #   {posterior::Pr(. <= 1)}
  #
  # BayesTIE <- lambda_s1 %>%
  #   {ifelse( test_stat > k, ., 1)}
  # # in old version, test_stat is an input matrix, and cannot
  # # be adjusted by user. In current version, the command is
  # # not able to be used for calculation. But I believe the
  # # method is correct. I will just keep the code in comment and
  # # use it in the future.
  #
  # pFDR <- lambda_s1 %>%
  #   .[test_stat > k] %>%
  #   mean() %>%
  #   {ifelse( is.na(.), 0, .)}


  return(
    list(
      k = k,
      optim = optim,
      pFDR = pFDR,
      test_stat = test_stat_mat,
      BayesTIE = tmp_prob_null_mat_full,
      range_test_stat = range_test_stat,
      sig_pfdr = dscov_mat
    )
  )
  # The output contains a matrix of Bayes type-I-error,
  # a matrix of dscov (0 or 1).
  # To simplify this function, these two results must be
  # included. I think using the posterior package can simplify
  # this function and make it more readable. I will do it in
  # the next version.
}

utils::globalVariables(c("."))

