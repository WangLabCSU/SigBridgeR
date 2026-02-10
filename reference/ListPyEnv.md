# List Available Python Environments

Discovers and lists available Python environments of various types on
the system. This generic function provides a unified interface to find
Conda environments and virtual environments (venv) through S3 method
dispatch.

Default method that lists all Python environments by combining results
from Conda and virtual environment discovery methods.

Discovers Conda environments using multiple detection strategies for
maximum reliability. First attempts to use system Conda commands, then
falls back to reticulate's built-in Conda interface if Conda command is
unavailable or fails. Returns empty data frame if Conda is not available
or no environments are found.

Discovers virtual environments by searching common venv locations
including user directories (`~/.virtualenvs`, `~/.venvs`) and project
folders (`./venv`, `./.venv`). Supports custom search paths through the
`venv_locations` parameter. Returns empty data frame if no virtual
environments are found in the specified locations.

## Usage

``` r
ListPyEnv(
  env_type = c("all", "conda", "venv", "virtualenv"),
  timeout = 30000L,
  venv_locations = c("~/.virtualenvs", "~/.venvs", "./venv", "./.venv"),
  verbose = TRUE,
  ...
)
```

## Arguments

- env_type:

  Character string specifying the type of environments to list. One of:
  `"all"`, `"conda"`, `"venv"`. Defaults to `"all"`.

- timeout:

  The maximum timeout time when using system commands, only effective
  when `env_type=conda`.

- venv_locations:

  Character vector of directory paths to search for virtual
  environments. Default includes standard locations and common project
  directories.

- verbose:

  Logical indicating whether to print verbose output.

- ...:

  For future use.

## Value

A data frame with the following columns:

- `name` - Character vector of environment names

- `python` - Character vector of paths to Python executables

- `type` - Character vector indicating environment type (`"conda"` or
  `"venv"`)

Returns an empty data frame with these columns if no environments are
found.

## Details

The function uses S3 method dispatch to handle different environment
types:

- **`"all"`**: Combines results from all environment types using
  [`rbind()`](https://rdrr.io/r/base/cbind.html)

- **`"conda"`**: Searches for Conda environments using multiple methods:

  - Primary:
    [`reticulate::conda_list()`](https://rstudio.github.io/reticulate/reference/conda-tools.html)
    for reliable environment detection

  - Fallback: System `conda info --envs` command for broader
    compatibility

- **`"venv"`**: Searches common virtual environment locations including
  user directories and project folders

Each method includes comprehensive error handling and will return empty
results with informative warnings if no environments are found or if
errors occur during discovery.

## Examples

``` r
if (FALSE) { # \dontrun{
# List all Python environments
ListPyEnv("all")

# List only Conda environments
ListPyEnv("conda")

# List only virtual environments with custom search paths
ListPyEnv("venv", venv_locations = c("~/my_envs", "./project_env"))
} # }
```
