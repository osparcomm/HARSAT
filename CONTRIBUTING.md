# Contributor's Guide

This document will help you set up access to the main repository for
the `harsat` code, and help you contribute to the code.

## Summary

* Main repository on Github: https://github.com/osparcomm/HARSAT
* Repository is currently *private* -- for access permissions, contact 
  [Chris Moulton](https://github.com/moultonc) or the OSPAR team, as they
  are the ones who administer the repository on Github. You will need this access to 
  both read the code and to contribute.
* Web documentation on Github Pages: https://osparcomm.github.io/HARSAT/


## Packaging

Our aim is to make the `harsat` code work as an R package. It is not
going to be distributed on CRAN for the near future at least. Instead,
it can be installed directly from Github. 

To install the latest development version, use the `remotes` package:

```
library(remotes)
remotes::install_github("osparcomm/harsat", auth_token = 'XXXX')
```

This should install all the `harsat` code with all its dependencies.

> Note: during development the repository is marked as private on GitHub, so you
> will need a GitHub Personal Access Token (or PAT) to access it. Put it in the
> call above as the 'XXXX' string.


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
a pull request using its web interface. At that point it will show up on the
main repository and we can collaborate with you to integrate it into the 
code.

Please note that we are particularly keen to improve the code quality. 


## Documentation

All documentation is held within the Github repository. We use the following
flow.

1. `roxygen2` is used for source-code documentation. We particularly welcome
   pull requests to improve this documentation. 
2. We have several vignettes in the `vignettes` directory. Some actually run `harsat` code (the ones 
   matching `*.Rmd.orig`) and are therefore *precompiled*, because they can take 15-20 minutes
   to run. That turns `*.Rmd.orig` into corresponding `*.Rmd` files. Then, the 
   normal documentation building for vignettes on installation will deploy these
   files, along with the basic `*.Rmd` vignettes that do not run `harsat` code.
3. A web site is built using `pkgdown` -- this happens automatically when pull
   requests are merged to main. This documents the code *as is* -- it does not 
   run `roxygen2` -- if you need to do that, you should do it manually.
4. When we merge pull requests, the web documentation is used to update our
   GitHub Pages site: https://osparcomm.github.io/HARSAT/


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

