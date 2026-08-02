You are an R package maintainer.

Generate a NEWS.md changelog entry for an R package based on the following git commits.

Requirements:

- Use Markdown.
- The top-level heading must be exactly: `# {{PACKAGE}} {{VERSION}}`
- The version number must not include the `v` prefix.
- Use second-level headings for sections, for example: `## NEW FEATURES`.
- Organize changes into the following sections when applicable:
  - NEW FEATURES
  - BUG FIXES
  - MINOR IMPROVEMENTS
  - DOCUMENTATION
  - DEPRECATED
  - DEFUNCT
  - BREAKING CHANGES
  - PERFORMANCE
  - TESTING
  - INTERNAL CHANGES
- Omit any section that has no relevant content.
- Do not mechanically copy commit messages one by one.
- Merge related commits into a single changelog bullet (e.g. "* Fixed bug") where appropriate.
- Ignore trivial commits that are not relevant to users, unless they belong in INTERNAL CHANGES.
- Do not invent information that is not supported by the commits.
- Do not include commit hashes, pull request numbers, or internal metadata unless they are essential.
- Use bullet points under each section.
- Prefer R package terminology where relevant, such as functions, arguments, datasets, vignettes, examples, tests, dependencies, CRAN checks, and documentation.
- You can include commit messages in the changelog entry. If a commit is ambiguous, describe it according to the diff summary.
- Besides the new NEWS.md entry, you can use following rules to NEWS.md
  - You can include markdown syntax (tables, links, bold font, italics, list, etc.) outside the changelog entry when necessary.
  - You can include code examples showing what has been changed or added, but do not include the entire file, code or function usage. The code chunk should be less than 30 lines, each line should be less than 80 characters. The code should be formatted with air style.
  - You can use emojis and kaomoji. Generally, one or two emojis per paragraph are fine.
- All news writtern in English with markdon syntax

Package: {{PACKAGE}}
Version: {{VERSION}}

Commits:

{{COMMITS}}

Changed files:

{{CHANGED_FILES}}

Diff summary:

{{DIFF_STAT}}