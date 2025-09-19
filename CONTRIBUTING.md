# Contributing to OPM Simulators

Thank you for your interest in contributing to the Open Porous Media (OPM) Simulators project! This guide will help you get started with contributing code, documentation, or bug reports.

## Table of Contents

- [Getting Started](#getting-started)
- [Development Process](#development-process)
- [Code Style and Formatting](#code-style-and-formatting)
- [Automated Formatting Tools](#automated-formatting-tools)
- [Testing](#testing)
- [Submitting Changes](#submitting-changes)
- [Reporting Issues](#reporting-issues)
- [Getting Help](#getting-help)

## Getting Started

Before contributing, please:

1. Read the [OPM Project Documentation](https://opm-project.org)
2. Ensure you have the necessary [prerequisites](https://opm-project.org/?page_id=239) installed

## Development Process

### Fork and Clone

1. Fork the repository on GitHub to your account
2. Clone the official OPM repository:
   ```bash
   git clone https://github.com/OPM/opm-simulators.git
   cd opm-simulators
   ```

3. Add your fork as a remote:
   ```bash
   git remote add fork https://github.com/YOUR-USERNAME/opm-simulators.git
   ```

   Note: If you already have a clone of the official repo, you can simply add your fork
   as a remote. The default 'origin' will track the official OPM repository.

### Building

Before making changes, ensure you can build the project successfully. For detailed build instructions, see the [OPM build documentation](https://opm-project.org/?page_id=231).

Quick build steps:
```bash
mkdir build
cd build
cmake ..
make -j$(nproc)
```

### Create a Feature Branch

Always create a new branch for your changes:
```bash
git checkout -b feature/your-feature-name
```

## Code Style and Formatting

Maintaining consistent code style is crucial for readability and maintainability. Please follow these guidelines:

### C++ Code Style

We use `clang-format` for C++ code formatting. Key settings include:

- **Line length**: Maximum 120 characters
- **Indentation**: 4 spaces (no tabs)
- **Braces**: Linux style (opening brace on next line for functions and classes)
- **Based on**: WebKit style with modifications

Apply clang-format to your C++ code:

**For new files:**
```bash
clang-format -i path/to/your/new_file.cpp
```

**For existing files:**
Only format the lines you modify to avoid creating noisy diffs. Most editors can format selections:
- **VS Code**: Select code → Right-click → "Format Selection"
- **vim**: Visual select → `=` to format
- **emacs**: Select region → `M-x clang-format-region`
- **Command line**: Use `git clang-format` to format only staged changes:
  ```bash
  git add your_modified_file.cpp
  git clang-format
  ```
  Note: `git clang-format` is a separate tool (Python script) that uses `clang-format` internally
  to format only the changed lines. Both tools need to be installed.

This ensures your functional changes aren't obscured by formatting changes elsewhere in the file.

### Include File Organization

Organize `#include` directives in the following order, with blank lines between groups:

1. **System headers** (C standard library, C++ STL):
   ```cpp
   #include <algorithm>
   #include <cmath>
   #include <iostream>
   ```

2. **External library headers** (Dune, etc.):
   ```cpp
   #include <dune/common/fmatrix.hh>
   #include <dune/common/fvector.hh>
   ```

3. **OPM headers from other modules**:
   ```cpp
   #include <opm/common/Exceptions.hpp>
   #include <opm/material/Constants.hpp>
   ```

4. **Local OPM headers** (from this module):
   ```cpp
   #include <opm/models/common/multiphasebaseproperties.hh>
   #include <opm/models/flash/flashproperties.hh>
   ```

Within each group, maintain **alphabetical order**.

### Python, Shell, and YAML Files

- Use **4 spaces** for indentation
- Follow [PEP 8](https://peps.python.org/pep-0008/) for Python code
- Use **lowercase with underscores** for shell script variables
- Maintain consistent quoting style in YAML

### CMake Files

- Use **2 spaces** for indentation
- Use lowercase for commands (e.g., `add_library` not `ADD_LIBRARY`)
- Group related commands together with blank lines for readability

## Automated Formatting Tools

We recommend using pre-commit hooks to automatically format your code before commits.

### What is Pre-commit?

Pre-commit is a framework for managing git hooks that automatically run code formatters and linters before each commit.
It only modifies files you're about to commit (staged files), not your entire codebase. The hooks fix common issues
like trailing whitespace and missing newlines automatically. When pre-commit makes changes, it will fail the commit
and report which files were modified - you'll need to review the changes, re-stage the files (`git add`), and commit
again. This ensures you see exactly what was changed before it's committed.

### Pre-commit Setup

1. **Install pre-commit**:
   ```bash
   pip install pre-commit
   # Or: conda install -c conda-forge pre-commit
   # Or: apt install pre-commit (Ubuntu/Debian)
   # Or: brew install pre-commit (macOS)
   ```

2. **Install hooks in your local repository**:
   ```bash
   pre-commit install
   ```

3. **What it does**:
   - Removes trailing whitespace (preserves Markdown line breaks)
   - Ensures files end with exactly one newline
   - Runs automatically on every commit
   - Supports C/C++, Python, Shell, YAML, CMake, and Markdown files

4. **Manual usage**:
   ```bash
   # Run on all files
   pre-commit run --all-files

   # Skip hooks for a specific commit
   git commit --no-verify

   # Update hook versions
   pre-commit autoupdate
   ```

## Testing

Before submitting your changes:

1. **Build the project**:
   ```bash
   mkdir build && cd build
   cmake ..
   make
   ```

2. **Run the test suite**:
   ```bash
   make test
   # Or for more verbose output:
   ctest --output-on-failure
   ```

3. **Run specific tests**:
   ```bash
   ctest -R test_name
   ```

4. **Ensure your changes don't break existing functionality**

## Submitting Changes

### Commit Messages

Write clear, descriptive commit messages:

- First line: Brief summary (50 characters or less)
- Blank line
- Detailed description (wrap at 72 characters)
- Reference relevant issues: `Fixes #123` or `Relates to #456`

Example:
```
Add support for polymer flooding in Flow

This commit introduces polymer flooding capabilities to the Flow
simulator, including:
- Polymer concentration tracking
- Viscosity modifications
- Adsorption modeling

Fixes #789
```

### Pull Request Process

1. **Update your branch** with latest upstream changes:
   ```bash
   git fetch origin
   git rebase origin/master
   ```

2. **Push your changes** to your fork:
   ```bash
   git push fork feature/your-feature-name
   ```

3. **Create a Pull Request** on GitHub:
   - Provide a clear title and description
   - Reference any related issues
   - Include test results if applicable
   - Be responsive to reviewer feedback

4. **Code Review**:
   - Address reviewer comments promptly
   - Add new commits for changes (don't force-push during review)
   - Once approved, your PR will be merged

## Reporting Issues

When reporting issues, please provide:

1. **Clear description** of the problem
2. **Steps to reproduce** the issue
3. **Expected behavior** vs actual behavior
4. **System information**:
   - OS and version
   - Compiler and version
   - OPM version or commit hash
5. **Build logs** if relevant (see below)
6. **Input data deck** demonstrating the issue, if possible

### Capturing Build Logs

To capture a build log for issue reporting:

```bash
LOGFILE=$(date +%Y%m%d-%H%M-)build.log
cmake -E cmake_echo_color --cyan --bold "Log file: $LOGFILE"
script -q $LOGFILE -c 'cmake .. -DCMAKE_BUILD_TYPE=Debug' &&
script -q $LOGFILE -a -c 'make -j 4' ||
cat CMakeCache.txt CMakeFiles/CMake*.log >> $LOGFILE
```

Upload the log file to [gist.github.com](https://gist.github.com) and include the link in your issue.

## Getting Help

### Resources

- **Documentation**: [opm-project.org](https://opm-project.org)
- **Issue Tracker**: [GitHub Issues](https://github.com/OPM/opm-simulators/issues)
- **Mailing List**: [OPM Mailing List](https://opm-project.org/?page_id=358)
- **Build Instructions**: [Build Guide](http://opm-project.org/?page_id=36)

### Community Guidelines

- Be respectful and constructive in discussions
- Search existing issues before creating new ones
- Provide context and be specific when asking for help
- Acknowledge contributions from others

## License

By contributing to OPM Simulators, you agree that your contributions will be licensed under the GNU General Public License v3.0 or later (GPLv3+).

---

Thank you for contributing to OPM Simulators! Your efforts help advance open-source reservoir simulation technology.
