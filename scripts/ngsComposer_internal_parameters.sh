
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
export cluster=${cluster//*=}

######################################################################################################################################################
# Software defined parameters

if [ "$cluster" == true ];then
	module unload R
  module add R
  Rversion=$((R --version) 2>&1)
  if [[ "$Rversion" =~ "R version" ]]; then
    echo -e "${white}\n- Using $Rversion\n ${white}"
  fi
fi
if [ "$cluster" == false ];then
  Rversion=$((R --version) 2>&1)
  if [[ "$Rversion" =~ "R version" ]]; then
    echo -e "${white}\n- Using $Rversion\n ${white}"
  else
    echo -e "${white}- install R before proceeding ${white}"
    echo -e "${white}- dependencies for R in linux: <sudo apt install libcurl4-openssl-dev> and <sudo apt install libssl-dev>"
  fi
fi



if [ "$cluster" == true ];then
	module unload python
  module add python/3.9.5
  pythonversion=$((python --version) 2>&1)
  if [[ "$pythonversion" =~ "Python 3" ]]; then
    echo -e "${white}\n- Using $pythonversion\n ${white}"
  else
    mkdir ~/bin
    PATH=~/bin:$PATH
		rm ~/bin/python
    ln -s /usr/bin/python3 ~/bin/python
  fi
fi
if [ "$cluster" == false ];then
  mkdir ~/bin
  PATH=~/bin:$PATH
	rm ~/bin/python
  ln -s /usr/bin/python3 ~/bin/python
  pythonversion=$((python --version) 2>&1)
  if [[ "$pythonversion" =~ "Python 3" ]]; then
    echo -e "${white}\n- Using $pythonversion\n ${white}"
  else
    echo -e "${white}- install python3 before proceeding ${white}"
  fi
fi


######################################################################################################################################################
# tools
export anemone=${ngsComposer_dir}/tools/anemone.py
export crinoid=${ngsComposer_dir}/tools/crinoid.py
export krill=${ngsComposer_dir}/tools/krill.py
export porifera=${ngsComposer_dir}/tools/porifera.py
export rotifer=${ngsComposer_dir}/tools/rotifer.py
export scallop=${ngsComposer_dir}/tools/scallop.py


if command -v pigz &>/dev/null; then
  export gzip=pigz
  export gunzip=unpigz
  export zcat="unpigz -c"
else
  export gzip=gzip
  export gunzip=gunzip
  export zcat=zcat
fi
