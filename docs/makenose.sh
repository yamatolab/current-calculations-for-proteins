#! /bin/bash



dir_path=$(cd $(dirname $0); pwd)
cd $dir_path

build_cmd=$dir_path/../buildout/bin/nosetests
args=" --with-doctest --doctest-extension=.rst "
src_path=$dir_path/source

$build_cmd $args $src_path

