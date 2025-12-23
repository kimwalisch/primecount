#!/bin/sh

if [ $# -ne 1 ]
then
    echo "Usage example:"
    echo "$ ./update_version.sh 3.5"
    echo "Update the primecount version to 3.5 in all files"
    exit 1
fi

# Run from primecount root directory
test -e ../src && cd ..

new_version=$1
old_version=$(grep "PRIMECOUNT_VERSION " include/primecount.hpp | cut -f2 -d'"')

old_major=$(echo $old_version | cut -f1 -d'.')
old_minor=$(echo $old_version | cut -f2 -d'.')

new_major=$(echo $new_version | cut -f1 -d'.')
new_minor=$(echo $new_version | cut -f2 -d'.')

old_year=$(grep "Copyright (c)" COPYING | cut -f5 -d' ' | cut -f1 -d',')
new_year=$(date +'%Y')

echo "Old version: $old_version"
echo "New version: $new_version"
echo ""
echo "Old year: $old_year"
echo "New year: $new_year"
echo ""

# 1. Update VERSION in CMakeLists.txt (Modern project() syntax)
echo "Update version in CMakeLists.txt"
sed "s/VERSION $old_version/VERSION $new_version/" CMakeLists.txt > CMakeLists.txt.tmp
mv -f CMakeLists.txt.tmp CMakeLists.txt

# Update version
for i in $(echo include/primecount.h \
                include/primecount.hpp)
do
    echo "Update version string in $i"
    sed "s/$old_major\.$old_minor/$new_version/g" $i > $i.tmp
    mv -f $i.tmp $i
done

# 3. Update version macros in headers (MAJOR/MINOR)
for i in $(echo include/primecount.h \
                include/primecount.hpp)
do
    echo "Update version macros in $i"
    sed "s/PRIMECOUNT_VERSION_MAJOR $old_major/PRIMECOUNT_VERSION_MAJOR $new_major/g" $i > $i.tmp
    mv -f $i.tmp $i
    sed "s/PRIMECOUNT_VERSION_MINOR $old_minor/PRIMECOUNT_VERSION_MINOR $new_minor/g" $i > $i.tmp
    mv -f $i.tmp $i
done

# 4. Update copyright year
for i in $(echo COPYING \
                src/app/help.cpp)
do
    echo "Update year in $i"
    sed "s/$old_year/$new_year/g" $i > $i.tmp
    mv -f $i.tmp $i
done

echo ""
echo "Version has been updated!"
