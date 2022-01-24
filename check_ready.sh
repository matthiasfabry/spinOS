version=$1

if grep -q "version-v$version" README.md; then
  echo "README ready for v$version"
else
  echo "README _not_ ready for v$version"
fi

if grep -q "v$version" svgs/spinOS.svg; then
  echo "SVG ready for v$version"
  ./make_logo.sh
else
  echo "SVG not ready for v$version"
fi

if grep -q "$version" modules/constants.py; then
  echo "modules ready for v$version"
else
  echo "modules not ready for v$version"
fi

