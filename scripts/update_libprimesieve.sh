#!/bin/bash

# Run from primecount root directory
test -e ../lib && cd ..

# Exit immediately if a command exits with a non-zero status
set -e

cd lib
rm -rf master.tar.gz*
rm -rf primesieve*
wget https://github.com/kimwalisch/primesieve/archive/master.tar.gz
tar xvf master.tar.gz
rm master.tar.gz*
mv primesieve-master primesieve

git status
echo ""

read -p "Do you want to commit and push the changes? [y/n] " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    git add .
    git commit -m "Update to latest libprimesieve"
    git push
fi
