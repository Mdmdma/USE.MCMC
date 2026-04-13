# Check that column names exist in a dataset

Check that column names exist in a dataset

## Usage

``` r
check_columns_exist(
  data,
  columns,
  data_name = "dataset",
  col_name = "dimensions"
)
```

## Arguments

- data:

  The data.frame or sf object to check columns in

- columns:

  Character vector of column names to look for

- data_name:

  Name of the data argument (for error messages)

- col_name:

  Name of the columns argument (for error messages)
