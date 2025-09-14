#!/bin/bash

SOURCE_PATHS=( 
    src/noise.cpp 
    src/math.cpp 
    src/image.cpp 
    src/readPNG/src/huffmanTree.cpp
    src/readPNG/src/readPNG.cpp 
    src/readPNG/src/safeRead.cpp 
    src/readPNG/src/datastream.cpp 
)

INCLUDE_DIRS=(

)

# Compile and run C++ program
echo "Compiling and running C++ program..."
# clang++ -std=c++20 "${SOURCE_PATHS[@]}" -O3 -o graphicsComp
clang++ -g -O0 -fsanitize=address,undefined -fno-omit-frame-pointer -std=c++23 graphics.cpp "${SOURCE_PATHS[@]}" -o graphicsComp
./graphicsComp

# Check if the C++ program ran successfully
if [ $? -eq 0 ]; then
    echo "C++ program executed successfully."

    # Run Python program
    echo "Running Python program..."
    source testenv/bin/activate
    python displayImg.py

    # Check if the Python program ran successfully
    if [ $? -eq 0 ]; then
        echo "Image program executed successfully."
    else
        echo "Error: Image program execution failed."
    fi
else
    echo "Error: C++ program execution failed."
fi

rm graphicsComp
