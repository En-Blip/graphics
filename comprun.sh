#!/bin/bash

# Compile and run C++ program
echo "Compiling and running C++ program..."
clang++ -std=c++17 graphics.cpp src/noise.cpp src/math.cpp src/image.cpp -o graphicsComp 
./graphicsComp

rm graphicsComp

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
