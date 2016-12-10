#!/bin/sh

if [ $# -ne 2 ]
then
    echo "Usage example:"
    echo "$ ./update_version.sh 1.2.3"
    echo "Updates the primecount version to 1.2.3 in all files"

    exit 1
fi

# Run from primecount root directory
test -e ../src && cd ..

new_version=$1
old_version=$(grep "PRIMECOUNT_VERSION " include/primecount.hpp | cut -f2 -d'"')

new_major=$(echo $new_version | cut -f1 -d'.')
new_minor=$(echo $new_version | cut -f2 -d'.')

old_major=$(echo $old_version | cut -f1 -d'.')
old_minor=$(echo $old_version | cut -f2 -d'.')

new_year=$(date +'%Y')
old_year=$(grep "Copyright (c)" COPYING | cut -f5 -d' ' | cut -f1 -d',')

echo "New version: $new_version"
echo "Old version: $old_version"
echo ""
echo "New year: $new_year"
echo "Old year: $old_year"
echo ""

# Update version
for i in $(echo README.md \
                doc/primecount-MPI.md \
                include/primecount.hpp)
do
    echo "Update version in $i"
    sed "s/$old_major\.$old_minor/$new_version/g" $i > $i.tmp
    mv -f $i.tmp $i
done

# Update version
for i in $(echo include/primecount.hpp)
do
    sed "s/PRIMECOUNT_VERSION_MAJOR $old_major/PRIMECOUNT_VERSION_MAJOR $new_major/g" $i > $i.tmp
    mv -f $i.tmp $i
    sed "s/PRIMECOUNT_VERSION_MINOR $old_minor/PRIMECOUNT_VERSION_MINOR $new_minor/g" $i > $i.tmp
    mv -f $i.tmp $i
done

# Update year
for i in $(echo COPYING \
                src/app/help.cpp)
do
    echo "Update year in $i"
    sed "s/$old_year/$new_year/g" $i > $i.tmp
    mv -f $i.tmp $i
done

echo "Version has been updated!"
