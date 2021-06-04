
# Compile SWMM C
echo
echo Compiling SWMM 5.13 ...
echo

if [ $DEBUG_API = true ]
then
    rm -r debug
    mkdir debug

    # Compile tests
    cp "$API_TEST_DIR/tests.c" "$API_DIR/src/"
    cp "$API_TEST_DIR/tests.h" "$API_DIR/src/"

    # DEBUG source files
    DEBUG_SOURCES="$API_TEST_DIR/interface_tests.f08"
fi

rm "$API_DIR/src/interface.h"
rm "$API_DIR/src/interface.c"
rm "$API_DIR/src/Makefile"

if [ -f "$API_DIR/src/tests.c" ]
then
    rm "$API_DIR/src/tests.c"
fi

if [ -f "$API_DIR/src/tests.h" ]
then
    rm "$API_DIR/src/tests.h"
fi

cp "$API_DIR/interface.h" "$API_DIR/src/"
cp "$API_DIR/interface.c" "$API_DIR/src/"

# Insert new files in SWMM C Makefile

SCRIPTS="interface.o"
OBJECTS="interface.o   : headers.h interface.h\n"

if [ $DEBUG_API = true ]
then
    SCRIPTS="${SCRIPTS} tests.o"
    OBJECTS="${OBJECTS}tests.o       : headers.h tests.h\n"
    cp "$API_TEST_DIR/tests.h" "$API_DIR/src/"
    cp "$API_TEST_DIR/tests.c" "$API_DIR/src/"
fi

sed "s#{{SCRIPTS}}#$SCRIPTS#" "$API_DIR/../Makefile" > "$API_DIR/src/Makefile"
sed -i "s#{{OBJECTS}}#$OBJECTS#" "$API_DIR/src/Makefile"

cd "$API_DIR/src" && make
cd .. && cp src/libswmm5.so ../../../libswmm5.so
cd ../../..
