# Contributor's Guide

## Contributing

### Issues

You can report issues with `harsat` using Github. Simply use the **Issues** tab
and choose **New issue**.

If you report an issue, we'd like you to answer the following questions.

1. What version of `harsat` are you using? (It's in the `DESCRIPTION` file)
2. What version of R are you using?
3. What operating system are you using?
4. What are you observing that appears to be incorrect, and what are you
   expecting to see? 

On this last point, the more information you can give us, the better. If you 
can copy and paste R console logs that helps tremendously, but if they're more
than 30 lines or so, it's better to attach them as a file, rather than pasting
them into the issue directly.


### Pull requests

We welcome pull requests. The easiest way to do this is to fork the repository
on Github, make the changes you want in your fork, push them to Github, and create
a pull request using its web interface. This will then be sent to the
[CODEOWNERS](link to: CODEOWNERS) for review and the pull request
approved/rejected as appropriate.


## Documentation

All documentation is held within the Github repository. We use the following
flow.

1. `roxygen2` is used for source-code documentation. 
2. We have several vignettes in the `vignettes` directory. Some run `harsat` code (the ones
   matching `*.Rmd.orig`) and are therefore *precompiled*, because they can take 15-20 minutes
   to run. That turns `*.Rmd.orig` into corresponding `*.Rmd` files. Then, the 
   normal documentation building for vignettes on installation will deploy these
   files, along with the basic `*.Rmd` vignettes that do not run `harsat` code.
3. A web site is built using `pkgdown` -- this happens automatically when pull
   requests are merged to main. This documents the code *as is* -- it does not 
   run `roxygen2` -- if you need to do that, you should do it manually.
4. When we push to `main`, the web documentation is used to update our
   GitHub Pages site: https://osparcomm.github.io/HARSAT/ -- this should happen
   automatically if you use a `git flow` style release process.

Note that the `pkgdown` action -- stored in Github under `.github/workflows/pkgdown.yaml`
also creates zip files for some of the more common data and configuration setups. 
These are then copied for deployment through Github Pages, so they can be downloaded
directly from the web pages. The actions which drive this zipping require a little
care, although they will continue to work fine if you simply change files in the 
data and reference file directories.

### Documentation build process

The documentation is not entirely built automatically. This is intentional, as,
for example, running some of the tools (e.g., `roxygen2`) can change key package
files in a breaking way. Also, some of the vignettes are pre-rendered, so they
do not take a long time to install for users. 

To update the documentation, use the following steps in R. This does change the 
namespacing, so API changes can cause the vignette builds to break. If so, 
create an issue.

**Step 1. Run roxygen2**

```r
library(devtools)
roxygen2::roxygenize()
```

**Step 2. Precompile the vignettes**

And then from bash:

```bash
R --quiet --vanilla < vignettes/precompile.R
```

**Step 3. Build and check the code**

Back to R:

```r
library(devtools)
devtools::check()
```

### Release process

The release process is more or less as follows. This is adapted from 
[this Stack Overflow answer](https://stackoverflow.com/a/45064206).

Let's assume the current version is `x.y.z`

```
git checkout develop
git flow release start x.y.z+1
emacs DESCRIPTION ## Update version x.y.z-9000 -> x.y.z+1
R --quiet --vanilla < vignettes/precompile.R
R CMD build .
R CMD check --no-manual harsat_x.y.z+1.tar.gz
git add DESCRIPTION
git commit -m "Updated version to x.y.z+1
git flow release finish x.y.z+1
git checkout master
git push
git push --tags
git push upstream
git push --tags upstream
git checkout develop
emacs DESCRIPTION ## Bump version to x.y.(z+1)-9000
git commit -am "Bump develop version [ci skip]"
git push
git push upstream
```
