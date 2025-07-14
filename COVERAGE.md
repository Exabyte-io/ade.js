# Code Coverage Guide

This document explains how to generate and view code coverage reports for the ade.js project.

## Quick Start

### Local Development
```bash
# Generate and view coverage report
npm run test:coverage:html
npm run test:coverage:view
```

### Check Coverage Thresholds
```bash
# Run tests and verify coverage meets minimum thresholds
npm run test:coverage:check
```

## Coverage Options

### 1. Local HTML Report
- **Command**: `npm run test:coverage:html`
- **Output**: `coverage/index.html`
- **View**: `npm run test:coverage:view` (opens in browser automatically)

### 2. GitHub Actions Coverage

#### PR Coverage Comments
- **Trigger**: Every pull request
- **What it does**: 
  - Runs tests with coverage
  - Posts a detailed coverage report as a PR comment
  - Shows coverage changes compared to base branch
  - Uploads coverage artifacts for download

#### Coverage Artifacts
- **Location**: GitHub Actions → Artifacts tab
- **Retention**: 30 days for main branch, 7 days for PRs
- **Download**: Click on artifact to download and view locally

#### GitHub Pages (Optional)
- **URL**: `https://exabyte-io.github.io/ade.js/`
- **Trigger**: Every push to main branch
- **Features**: 
  - Always up-to-date coverage report
  - No authentication required
  - Easy sharing with team

### 3. External Services

#### Codecov Integration
- **Service**: [Codecov](https://codecov.io)
- **Features**:
  - Historical coverage tracking
  - Coverage trends and graphs
  - PR coverage comparison
  - Coverage badges
- **Setup**: Automatic via GitHub Actions

## Coverage Thresholds

The project enforces these minimum coverage levels:

| Metric | Threshold |
|--------|-----------|
| Statements | 85% |
| Branches | 80% |
| Functions | 80% |
| Lines | 85% |

## Workflow Files

- `.github/workflows/coverage.yml` - Basic coverage with artifacts
- `.github/workflows/coverage-pages.yml` - GitHub Pages publishing
- `.github/workflows/pr-coverage.yml` - PR coverage comments

## Configuration

### nyc Configuration (`.nycrc`)
```json
{
  "reporter": ["text", "html", "lcov"],
  "check-coverage": true,
  "branches": 80,
  "lines": 85,
  "functions": 80,
  "statements": 85
}
```

### Coverage Files
- `coverage/lcov.info` - LCOV format for external services
- `coverage/index.html` - HTML report for browser viewing
- `coverage/` - Complete coverage report directory

## Troubleshooting

### Coverage Report Not Found
```bash
# Generate coverage first
npm run test:coverage:html

# Then view it
npm run test:coverage:view
```

### GitHub Pages Not Working
1. Check repository settings → Pages
2. Ensure GitHub Actions has Pages write permissions
3. Verify workflow is running on main branch pushes

### Coverage Thresholds Failing
1. Run `npm run test:coverage:check` to see specific failures
2. Add tests for uncovered code paths
3. Consider adjusting thresholds if appropriate

## Best Practices

1. **Always check coverage locally** before pushing
2. **Review PR coverage comments** to understand impact
3. **Use coverage artifacts** for detailed analysis
4. **Monitor coverage trends** over time
5. **Aim for meaningful coverage** rather than just hitting thresholds 
