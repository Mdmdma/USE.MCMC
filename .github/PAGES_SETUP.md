# GitHub Pages Setup for USE.MCMC

This document explains the GitHub Pages setup for automatically publishing package documentation.

## Overview

The package documentation is automatically built and published using:
- **pkgdown**: R package for building package websites
- **GitHub Actions**: Automated workflow that builds and deploys the site
- **GitHub Pages**: Hosting platform for the documentation

## Configuration

### Automatic Deployment Triggers

The documentation website is automatically rebuilt and deployed when:
1. **Push to main/master branch**: Any commit pushed to the main or master branch will trigger a rebuild
2. **Pull request to main/master**: Pull requests will build the site (but won't deploy) for preview
3. **Release published**: Creating a new release will trigger a rebuild
4. **Manual trigger**: You can manually trigger a rebuild from the Actions tab

### Workflow Location

The GitHub Actions workflow is located at: `.github/workflows/pkgdown.yaml`

### Website Configuration

The pkgdown configuration is in: `_pkgdown.yml`

The website URL is: https://mdmdma.github.io/USE.MCMC/

## Enabling GitHub Pages

To enable GitHub Pages for the first time (if not already enabled):

1. Go to your repository on GitHub
2. Click on **Settings**
3. Scroll down to **Pages** in the left sidebar
4. Under **Source**, select:
   - Branch: `gh-pages`
   - Folder: `/ (root)`
5. Click **Save**

## How It Works

1. When changes are pushed to main/master, the GitHub Actions workflow starts
2. The workflow:
   - Sets up R and required dependencies
   - Installs the package
   - Builds the documentation website using pkgdown
   - Deploys the generated site to the `gh-pages` branch
3. GitHub Pages serves the content from the `gh-pages` branch

## Viewing Workflow Runs

To check the status of documentation builds:
1. Go to the **Actions** tab in your GitHub repository
2. Click on the **pkgdown.yaml** workflow
3. View recent runs and their status

## Troubleshooting

If the documentation doesn't update:
1. Check the Actions tab for failed workflow runs
2. Ensure GitHub Pages is enabled in repository settings
3. Verify the `gh-pages` branch exists
4. Check that the workflow has write permissions

## Manual Rebuild

To manually rebuild the documentation:
1. Go to the **Actions** tab
2. Select **pkgdown.yaml** workflow
3. Click **Run workflow**
4. Select the branch (usually main/master)
5. Click **Run workflow**
