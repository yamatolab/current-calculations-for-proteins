#! /bin/bash



dir_path=$(cd $(dirname $0); pwd)
cd $dir_path

build_cmd=$dir_path/../buildout/bin/sphinx-build
src_path=$dir_path/source
dst_path=$dir_path/_build/doctest

$build_cmd -b doctest $src_path $dst_path

