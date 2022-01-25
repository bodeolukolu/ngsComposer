#!/bin/bash

red=$'\e[1;31m'
green=$'\e[1;32m'
yellow=$'\e[1;33m'
blue=$'\e[1;34m'
magenta=$'\e[1;35m'
cyan=$'\e[1;36m'
white=$'\e[0m'


main () {
  echo -e "${white}\n############################################## ${orange}\n- check for R installation ${white}\n##############################################${white}"
  if R --version; then
    :
  else
    module add R
    if R --version; then
      :
    else
      echo -e "${white}- install R before proceeding ${white}"
      echo -e "${white}- dependencies for R in linux: <sudo apt install libcurl4-openssl-dev> and <sudo apt install libssl-dev>"
    fi
  fi
}
main &>> ./log.out



main () {
echo -e "${blue}\n############################################## \n- installing R-package: ggplot2  ${blue}\n##############################################${white}"
  mkdir -p ./helpers/R_packages
  cd ./helpers/R_packages
  R -e 'install.packages("ggplot2", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("data.table", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("gtable", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("rlang", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("plyr", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("Rcpp", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("scales", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("R6", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("lifecycle", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("munsell", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("colorspace", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("glue", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("tibble", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("ellipsis", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("magrittr", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("pillar", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("fansi", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("utf8", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("vctrs", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("crayon", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("pkgconfig", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("withr", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("farver", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("digest", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  R -e 'install.packages("labeling", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
}
dirtool=./helpers/R_packages/ggplot2
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of R-package: ggplot2 ${white}"
  time main &>> ./log.out
fi
