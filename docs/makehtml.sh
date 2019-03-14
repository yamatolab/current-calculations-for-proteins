#! /bin/bash



dir_path=$(cd $(dirname $0); pwd)
cd $dir_path

build_cmd=$dir_path/../buildout/bin/sphinx-build
src_path=$dir_path/source
html_path=$dir_path/../../bbwiki/html

$build_cmd -b html $src_path $html_path

if [ $# -eq 1 ]; then
    commit_message=$1
    cd $html_path
    hg addremove ./
    hg commit -m "$commit_message"
fi
    
