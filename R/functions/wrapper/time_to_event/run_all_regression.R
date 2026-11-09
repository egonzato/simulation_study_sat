run_all_regression=function(data,methods,landmark=730){
  out=bind_rows(
    lapply(
      seq_len(length(methods)),
      function(i){
        run_regression_methods(
          data=data,
          methods=methods[i],
          landmark=landmark
        )
      }
    )
  )
}

