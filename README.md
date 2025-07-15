[![npm version](https://badge.fury.io/js/%40exabyte-io%2Fade.js.svg)](https://badge.fury.io/js/%40exabyte-io%2Fade.js)
[![License: Apache](https://img.shields.io/badge/License-Apache-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

# ade.js

ade.js houses entity definitions for use in the Mat3ra platform.


### Installation

For usage within a javascript project:

```bash
npm install @exabyte-io/ade.js
```

For development:

```bash
git clone https://github.com/Exabyte-io/ade.js.git
```


### Contribution

This repository is an [open-source](LICENSE.md) work-in-progress and we welcome contributions.

We regularly deploy the latest code containing all accepted contributions online as part of the
[Mat3ra.com](https://mat3ra.com) platform, so contributors will see their code in action there.

See [ESSE](https://github.com/Exabyte-io/esse) for additional context regarding the data schemas used here.

Useful commands for development:

```bash
# run linter without persistence
npm run lint

# run linter and save edits
npm run lint:fix

# compile the library
npm run transpile

# run tests
npm run test

# run tests with coverage
npm run test:coverage

# run tests with coverage and check thresholds
npm run test:coverage:check

# generate HTML coverage report
npm run test:coverage:html
```

## Code Coverage

This project includes comprehensive code coverage reporting with multiple viewing options:

### Local Coverage
- Run `npm run test:coverage:html` to generate an HTML coverage report locally
- Open `coverage/index.html` in your browser to view the report

### GitHub Integration
The project uses GitHub Actions to automatically generate and display coverage reports:

1. **PR Coverage Comments**: Every pull request automatically gets a coverage report comment showing:
   - Overall coverage percentages
   - Coverage changes compared to the base branch
   - Detailed file-by-file coverage breakdown

2. **Coverage Artifacts**: Coverage reports are uploaded as GitHub artifacts for each PR and commit
   - Download from the Actions tab in GitHub
   - Available for 30 days for main branch, 7 days for PRs

3. **GitHub Pages** (Optional): Coverage reports are published to GitHub Pages for easy browser viewing
   - Available at: `https://exabyte-io.github.io/ade.js/`
   - Updated on every push to main branch

### Coverage Thresholds
The project enforces minimum coverage thresholds:
- **Statements**: 85%
- **Branches**: 80%
- **Functions**: 80%
- **Lines**: 85%

### External Coverage Services
- **Codecov**: Coverage data is automatically uploaded to Codecov for historical tracking and trend analysis

ADe
===

The`ADe` package sits just below the `WoDe` package in the Mat3ra workflow
ecosystem, where `ADe` houses entity definitions for:

- `Application` - uniquely determined by `name, [version], [build]`
- `Executable` - defined for a given application and accessible from application by name
- `Flavor` - defined for a given executable and accessible from executable by name
- `Template` - a jinja template for an application input file

The relevant data parameterizing these entities is housed in
the [Application Flavors](https://github.com/Exabyte-io/exabyte-application-flavors)
repository. This includes the supported applications, executables, flavors,
and defined templates.

Templates themselves are organized by application in a top-level `assets`
directory in `application-flavors` and the API for loading and working with templates can be found in
each application's `assets.js` module.
At build time, all templates are loaded and compiled into a single monolithic
JS file using `build_templates.js` so that it can be used in the client as well as in NodeJS.
This is how templates are consumed from `applicaton-flavors` in `ADe`.
