run_all_regression=function(data,methods,id=NULL,landmark=2){
  out=bind_rows(
    lapply(
      seq_len(length(methods)),
      function(i){
        run_regression_methods(
          data=data,
          id=id,
          methods=methods[i],
          landmark=landmark
        )
      }
    )
  )
}

