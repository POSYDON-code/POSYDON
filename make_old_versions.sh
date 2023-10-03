versions=`git tag -l`

for version in $versions; do
  git checkout tags/$version
  cd docs
  rm -rf build
  make html
  cd ..
  git checkout gh-pages-dev
  mkdir $version
  cp -r docs/build/html/* $version/
  git add $version
  git commit -m "generate documentation for version $version"  
done
