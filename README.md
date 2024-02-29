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
```

ADe
===

The`ADe` package sits just below the `WoDe` package in the Mat3ra workflow
ecosystem, where `ADe` houses entity definitions for:

- `Application` - uniquely determined by `name, [version], [build]`
- `Executable` - defined for a given application and accessible from application by name
- `Flavor` - defined for a given executable and accessible from executable by name
- `Template` - a jinja template for an application input file

The relevant data parameterizing these entities is housed in
the [Application Flavors](https://github.com/Exabyte-io/application-flavors)
repository. This includes the supported applications, executables, flavors,
and defined templates.

Templates themselves are organized by application in a top-level `assets`
directory in `application-flavors` and the API for loading and working with templates can be found in
each application's `assets.js` module.
At build time, all templates are loaded and compiled into a single monolithic
TS file using `build_application_trees.ts` so that it can be used in the client as well as in NodeJS.
This is how templates are consumed from `applicaton-flavors` in `ADe`.
