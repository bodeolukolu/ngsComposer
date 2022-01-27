
red=$'\e[1;31m'
green=$'\e[1;32m'
yellow=$'\e[1;33m'
blue=$'\e[1;34m'
magenta=$'\e[1;35m'
cyan=$'\e[1;36m'
white=$'\e[0m'

######################################################################################################################################################
export ngsComposer_dir=${ngsComposer_dir%/*}
export projdir=${projdir%/*}

cluster="$(grep cluster config.sh)"
cluster=${cluster//*=}

######################################################################################################################################################
# Software defined parameters

Rout=$(R --version | head -n 3)
if [ -z "$Rout" ];then
  echo -e "${white}- R not available ${white}"
	module add R
	R --version | head -n 3
fi
Rout=$(R --version | head -n 3)
if [ -z "$Rout" ];then
  echo -e "${white}- install R before proceeding ${white}"
  echo -e "${white}- dependencies for R in linux: <sudo apt install libcurl4-openssl-dev> and <sudo apt install libssl-dev>"
fi

pythonout=$(python3 --version | head -n 3)
if [ -z "$pythonout" ];then
	module add python3
	python3 --version | head -n 3
fi
pythonout=$(python3 --version | head -n 3)
if [ -z "$pythonout" ];then
  echo -e "${white}- install Python3 before proceeding ${white}"
fi

######################################################################################################################################################
# tools
export anemone=${ngsComposer_dir}/tools/anemone.py
export crinoid=${ngsComposer_dir}/tools/crinoid.py
export krill=${ngsComposer_dir}/tools/krill.py
export porifera=${ngsComposer_dir}/tools/porifera.py
export rotifer=${ngsComposer_dir}/tools/rotifer.py
export scallop=${ngsComposer_dir}/tools/scallop.py

pythonout=$(python --version | head -n 3)
if [[ "$pythonout" =~ python3 ]]; then
  echo -e "${white}- Using $pythonout ${white}"
else
  mkdir ~/bin
  PATH=~/bin:$PATH
  ln -s /usr/bin/python3 ~/bin/python
fi


if command -v pigz &>/dev/null; then
  export gzip=pigz
  export gunzip=unpigz
  export zcat="unpigz -c"
else
  export gzip=gzip
  export gunzip=gunzip
  export zcat=zcat
fi
