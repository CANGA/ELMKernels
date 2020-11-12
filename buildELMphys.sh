if [ -d "build" ]; then 
  rm -rf "build"
fi
ORIGIN_DIR=${pwd}
mkdir "build" ; cd "build"
cmake ../
  -DCMAKE_CXX_FLAGS="-W -Wall -Wextra" \
make VERBOSE=1
make install
cd $ORIGIN_DIR