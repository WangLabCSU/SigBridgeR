# Create or Use Python Environment with Required Packages

Sets up a Python environment with specified packages for DEGAS screening
methods. This function can create new environments or reuse existing
ones, supporting both Conda and venv environment types. It ensures all
required dependencies are properly installed and verified.

Default method for unsupported environment types. Throws an informative
error with supported environment types.

## Usage

``` r
SetupPyEnv(env_type = c("conda", "venv"), ...)
```

## Arguments

- env_type:

  Character string specifying the type of Python environment to create
  or use. One of: \`"conda"\`, \`"venv"\`.

- ...:

  Additional parameters passed to specific environment methods.

## Value

A data frame containing verification results for the environment setup,
including installation status of all required packages. Invisibly
returns the verification results.

## Details

This function provides a comprehensive solution for Python environment
management in R projects, particularly for machine learning workflows
requiring TensorFlow. Key features include:

\- \*\*Environment Creation\*\*: Automatically creates new environments
or reuses existing ones with the same name - \*\*Package Management\*\*:
Installs specified Python packages with version pinning support -
\*\*Verification\*\*: Validates environment setup and package
installations - \*\*Flexible Methods\*\*: Supports different backend
methods for environment creation (reticulate vs system calls)

The function uses S3 method dispatch to handle different environment
types, allowing for extensible support of additional environment
managers in the future.

## See also

\[reticulate::conda_create()\], \[reticulate::virtualenv_create()\] for
underlying environment creation functions.

## Examples

``` r
if (FALSE) { # \dontrun{
# Setup a Conda environment with default parameters
SetupPyEnv("conda")

# Setup a venv environment
SetupPyEnv("venv")
} # }
```
