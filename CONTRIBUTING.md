# Contributing to SigBridgeR

Thank you for your interest in contributing to SigBridgeR! We welcome contributions in all forms — bug reports, feature suggestions, documentation improvements, and code contributions.

Please take a moment to review this guide to make the process smooth and efficient for everyone.

## Table of Contents

- [Contributing to SigBridgeR](#contributing-to-sigbridger)
  - [Table of Contents](#table-of-contents)
  - [Code of Conduct](#code-of-conduct)
  - [How to Contribute](#how-to-contribute)
    - [Report Bugs](#report-bugs)
    - [Suggest Features](#suggest-features)
    - [Improve Documentation](#improve-documentation)
    - [Contribute Code](#contribute-code)
  - [Development Setup](#development-setup)
  - [Pull Request Process](#pull-request-process)
  - [Getting Help](#getting-help)

## Code of Conduct

By participating in this project, you agree to abide by the [Contributor Covenant](https://www.contributor-covenant.org/). Please be respectful and constructive in all interactions.

## How to Contribute

### Report Bugs

If you find a bug, please open an [issue](https://github.com/WangLabCSU/SigBridgeR/issues) with the following information:

- A clear, descriptive title.
- Steps to reproduce the issue (include code snippets and data if possible).
- Expected vs. actual behavior.
- Your R session info (`sessionInfo()`).
- Package version (`packageVersion("SigBridgeR")`).

### Suggest Features

Have an idea for a new feature or an improvement? Open an [issue](https://github.com/WangLabCSU/SigBridgeR/issues) and label it as a feature request. Describe:

- The problem you're trying to solve.
- How you envision the feature working.
- Any relevant context or use cases.

### Improve Documentation

Documentation improvements are always welcome! This includes:

- Fixing typos or unclear wording in function docs or vignettes.
- Adding examples or clarifying edge cases.
- Translating content.

Simply submit a pull request with your changes.

### Contribute Code

We encourage contributions to extend SigBridgeR's functionality. The main areas where you can contribute code are:

| Area | Description |
|------|-------------|
| **New screening methods** | Register custom phenotype-associated cell screening algorithms via the extensible registry system. |
| **New annotation methods** | Register custom cell-type annotation algorithms. |
| **Bug fixes & optimizations** | Fix issues or improve performance of existing code. |

For detailed instructions on how to implement and register custom algorithms — including format requirements, validation, and registration — please refer to the [Extending SigBridgeR: A Guide for Custom Extensions](https://wanglabcsu.github.io/SigBridgeR/articles/Extending_ScreenMethod.html) vignette (`vignettes/Extending.Rmd`).

## Development Setup

1. **Fork and clone** the repository:

   ```bash
   git clone https://github.com/<your-username>/SigBridgeR.git
   cd SigBridgeR
   ```

2. **Install development dependencies** (recommended):

   code formatting: [air](https://github.com/posit-dev/air)   
   lints for R (choose one): 
      - [flir](https://github.com/etiennebacher/flir)
      - [lintr](https://github.com/r-lib/lintr)
      - [jarl](https://github.com/etiennebacher/jarl)
  
   install R code tools:

   ```r
   pak::pak(c(
      "tictoc",
      "codetools",
      "lintr",
      "rstudioapi"
   ))
   ```


3. **Load the package** in your R session:

   ```r
   devtools::document()
   ```

4. **Checking** everything is working:

   ```r
   devtools::check()
   ```

## Pull Request Process

1. **Create a feature branch** from `main`:

   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**, keeping them focused and well-scoped. Follow the existing code style.

3. **Write or update tests** as appropriate.

4. **Commit your changes** with a clear and descriptive commit message.

5. **Push to your fork** and open a pull request against the `main` branch.

   In your pull request description, please:
   - Summarize the changes and motivation.
   - Reference any related issues.
   - Mention any breaking changes or special considerations.

6. A maintainer will review your PR and may request changes. Please respond promptly to feedback.

## Getting Help

- **Issues:** Use [GitHub Issues](https://github.com/WangLabCSU/SigBridgeR/issues) for bug reports and feature requests.
- **Discussions:** Use [GitHub Discussions](https://github.com/WangLabCSU/SigBridgeR/discussions) for questions and general help.
- **Vignette:** See the [Extending SigBridgeR](https://wanglabcsu.github.io/SigBridgeR/) vignette for detailed implementation guidance.

---

Thank you for helping make SigBridgeR better!
