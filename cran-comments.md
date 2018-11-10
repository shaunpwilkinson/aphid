# aphid version 1.3.1

This is a patch release addressing two minor bugs in `align.list`

### Test environments

 * local ubuntu 16.04.4 x86_64-pc-linux-gnu; R version 3.5.1 (2018-07-02)
 * travis-ci ubuntu 14.04.5 x86_64-pc-linux-gnu; R version 3.5.1 (2017-01-27)
 * winbuilder R Under development (2018-11-08 r75566)

### R CMD check results

There were no ERRORs. 

There was one WARNING on local test:

checking compilation flags used ... WARNING
Compilation used the following non-portable flag(s):
  ‘-Wdate-time’ ‘-Werror=format-security’ ‘-Wformat’

According to this thread, this is a local issue and won't affect CRAN
<https://github.com/ropensci/stplanr/issues/260>


There was one NOTE:

Note_to_CRAN_maintainers
Maintainer: 'Shaun Wilkinson <shaunpwilkinson@gmail.com>'
  

### Downstream dependencies

There were no issues.
