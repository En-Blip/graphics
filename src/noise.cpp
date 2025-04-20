#include "noise.h"
#include "math.h"

#define MAX_COEFF (pow(2, 0.5))
//#define DEBUG

Noise::Noise(int2 dimensions){
    dims = dimensions;
}

Noise::~Noise(){
}

unsigned char* Noise::genPerlin(int layers, int2 frequency, double (*weightCalc)(int nLayers, int layer), bool reversed){
    // get a pointer function to store our return
    int* weights = new int[dims.x * dims.y];

    // for each layer generate a perlin array
    // store it in a temporary array before 
    unsigned char* tempWeights = new unsigned char[dims.x * dims.y];
    int2 currentFrequency;
    for(int i = 0; i < layers; i++){
        // generate the perlin noise
        currentFrequency = {static_cast<int>(floor(frequency.x / (i+1))), static_cast<int>(floor(frequency.y / (i+1)))};
        if (currentFrequency.x == 0){ currentFrequency = {1, 1}; }
        tempWeights = perlinNoise(dims, currentFrequency);

        // add the perlin noise to the weights
        for(int j = 0; j < dims.x * dims.y; j++){
            if (i == layers-1) {
                weights[j] += int(tempWeights[j] * weightCalc(layers, i));
            } else {
                weights[j] += int(tempWeights[j] * weightCalc(layers, i+1));
            }
        }
    }

    delete[] tempWeights;

    // normalize the weights then return them
    // could always make it optional
    int2 bounds = {0, 255};
    for (int i = 0; i < dims.x * dims.y; i++){
        if (weights[i] > bounds.x){
            bounds.x = weights[i];
        }
        if (weights[i] < bounds.y){
            bounds.y = weights[i];
        }
    }

    unsigned char* noise = new unsigned char[dims.x * dims.y];

    // loop through all the pixels and normalize them
    for (int i = 0; i < dims.x * dims.y; i++){
        noise[i] = linearBlend(weights[i], {bounds.x, bounds.y}, {255*reversed, 255-255*reversed});
    }

    return noise;
}

    /**
     * Generates a worley noise image. The worley noise is a value that shows
     * the distance from the closest point in the grid to the current pixel.
     * The distance is then normalized to a value from 0 to 255.
     *
     * @attention for it to tile properly, make sure the highest frequency divides the tile size
     * @param layers The number of noise layers to generate.
     * @param tileSize The size of each the tile that will be tiled to fit the dimensions.
     * @param frequency The highest frequency of the noise.
     * @param weightCalc A function to calculate the weight of each layer.
     * @param reversed A boolean value to reverse the colors of the noise.
     * @return A pointer to an array of integers representing the noise.
     */
unsigned char* Noise::genWorley(int layers, int tileSize, int2 frequency, double (*weightCalc)(int nLayers, int layer), bool reversed){
    unsigned char* noise = new unsigned char[dims.x * dims.y];
    if(tileSize / frequency.x != (double) tileSize / frequency.x){
        std::cout << "frequency does not divide tilesize, tiling will cause artifacts" << std::endl;
    }

    // get a pointer function to store our return
    unsigned char* weights = new unsigned char[tileSize * tileSize];

    // for each layer generate a worley array
    // store it in a temporary array before 
    unsigned char* tempWeights = new unsigned char[tileSize * tileSize];
    int currentFrequency;
    for(int i = 0; i < layers; i++){
        // generate the worley noise
        currentFrequency = static_cast<int>(floor(frequency.x / (i+1)));
        if (currentFrequency == 0){ currentFrequency = 1; }
        tempWeights = worleyNoise({tileSize, tileSize}, currentFrequency);

        // add the worley noise to the weights
        for(int j = 0; j < tileSize * tileSize; j++){
            if (i == layers-1) {
                weights[j] += int(tempWeights[j] * weightCalc(layers, i));
            } else {
                weights[j] += int(tempWeights[j] * weightCalc(layers, i+1));
            }
        }
    }

    delete[] tempWeights;

    // normalize the weights then return them
    // could always make it optional
    int2 bounds = {0, 255};
    for (int i = 0; i < tileSize * tileSize; i++){
        if (weights[i] > bounds.x){
            bounds.x = weights[i];
        }
        if (weights[i] < bounds.y){
            bounds.y = weights[i];
        }
    }
    #ifdef DEBUG
    std::cout << "got to norm" << std::endl;
    #endif

    // loop through all the pixels and normalize them
    for (int i = 0; i < dims.x; i++){
        for (int j = 0; j < dims.y; j++){
            // set the noise to the normalized value of the tile (and tile it to the dims size)
            noise[j * dims.x + i] = linearBlend(weights[(j%tileSize) * tileSize + i%tileSize], {bounds.x, bounds.y}, {255*reversed, 255-255*reversed});
        }
    }


    delete[] weights;

    return noise;
}

unsigned char* Noise::perlinNoise(int2 dims, int2 cellcount){
    // create points based on the number of cells evenly spaced around the dimensions
    // and populate them with a random unit vector
    float2* pointVectors = new float2[cellcount.x * cellcount.y];
    int len_vectors = static_cast<int>(sizeof(vectors) / sizeof(float2));

    // get some basic values like cellsize
    int2 cellsize = {static_cast<int>(floor(dims.x/(cellcount.x))), static_cast<int>(floor(dims.y/(cellcount.y)))};

    // store a variable to find the minimum and maximum pixel values to normalize
    int2 bounds = {0, 255};


    for (int i = 0; i < cellcount.x * cellcount.y; i++){
        // choose a random unit vector
        
        pointVectors[i] = vectors[rand() % len_vectors];
    }
    
    // create the noise function based on the point vectors
    int* shade = new int[dims.x * dims.y];

    // populate the noise function
    for (int i = 0; i < dims.x * dims.y; i++){
        // calculate the pixel, and the cell that the current pixel is in
        int2 pixel = {i % dims.x, static_cast<int>(floor(i / dims.x))}; // check if int division works with this
        int2 cell;

        // for each of the 4 corners caclulate the dot product of the corner
        float* dotProducts = new float[4];
        int* relDistances = new int[4];
        int2* relPos = new int2[4];
        for (int y = 0; y < 2; y++){
            for (int x = 0; x < 2; x++){

                // calculate the cell that the pixel is in and loop around if necessary
                cell = {(static_cast<int>(ceil(pixel.x / cellsize.x)) + x), (static_cast<int>(ceil(pixel.y / cellsize.y)) + y)};

                // calculate the relative position of the pixel in the cell (for each corner)
                relPos[y*2+x] = {pixel.x - cellsize.x * cell.x, pixel.y - cellsize.y * cell.y};

                // loop the cell values now that we have them in the proper position
                cell = {cell.x % cellcount.x, cell.y % cellcount.y};

                // calculate the dot product, and the magnitude of each vector
                dotProducts[y*2+x] = dot(pointVectors[cell.y * cellcount.x + cell.x], {static_cast<float>(relPos[y*2+x].x), static_cast<float>(relPos[y*2+x].y)});
                relDistances[y*2+x] = sqrt(pow(relPos[y*2+x].x, 2) + pow(relPos[y*2+x].y, 2));
            } 
        }

        // cubicly interpolate the dot products
        float yContribution, xContribution;
        float* contributions = new float[4];
        for (int j = 0; j < 4; j++){
            yContribution = cubicBlend(abs(relPos[j].y), {0, static_cast<float>(cellsize.y)}, {0, 1});
            xContribution = cubicBlend(abs(relPos[j].x), {0, static_cast<float>(cellsize.x)}, {0, 1});
            contributions[3-j] = xContribution * yContribution;
        }

        int pixelValue = 0;
        for (int j = 0; j < 4; j++){
            pixelValue += dotProducts[j] * contributions[j];
        }
        shade[i] = pixelValue;

        if (shade[i] < bounds.y) {
            bounds.y = shade[i];
        }
        if (shade[i] > bounds.x) {
            bounds.x = shade[i];
        }

        delete[] dotProducts;
        delete[] relDistances;
        delete[] relPos;
        delete[] contributions;

    }
    unsigned char* noise = new unsigned char[dims.x * dims.y];

    // loop through all the pixels and normalize them
    for (int i = 0; i < dims.x * dims.y; i++){
        noise[i] = linearBlend(shade[i], {bounds.x, bounds.y}, {0, 255});
    }

    delete[] pointVectors; 

    return noise;
}

unsigned char* Noise::worleyNoise(int2 dims, int frequency){
    // requires dims (tilesize) to be square
    if (dims.x != dims.y){
        std::cerr << "noise.cpp: worleyNoise: dims must be square" << std::endl;
        return nullptr;
    }

    // create a variable to store our return
    unsigned char* noise = new unsigned char[dims.x * dims.y];

    // create a freq x freq grid to store random points 
    int2* grid = new int2[frequency * frequency];

    // variable for cellsize
    int2 cellsize = {dims.x / frequency, dims.y / frequency};

    // populate the grid
    for (int i = 0; i < frequency * frequency; i++){
        grid[i] = {rand() % (cellsize.x), rand() % (cellsize.y)};
    }

    // some variables to be declared 
    int2 distVector;
    int distance, closestDist = -1;
    int2 currentCell;
    int cellSearchDims = (frequency > 1)?2:1;

    // foreach pixel in a single tile
    for (int i = 0; i < dims.x; i++){
        for (int j = 0; j < dims.y; j++){
            // find the closest point in the grid
            // linearly interpolate using that distance and max (0, 255)
            for (int celli = -cellSearchDims; celli <= cellSearchDims; celli++){
                for (int cellj = -cellSearchDims; cellj <= cellSearchDims; cellj++){
                    // find current cell
                    currentCell = {i / cellsize.x + celli, j / cellsize.y + cellj};

                    // perform bound looping
                    if (currentCell.x < 0) { currentCell.x += frequency; }
                    if (currentCell.y < 0) { currentCell.y += frequency;  }
                    if (currentCell.x >= frequency) { currentCell.x -= frequency; }
                    if (currentCell.y >= frequency) { currentCell.y -= frequency; }

                    // offset the distance in this case to accuratly represent the distance
                    if((currentCell.x + (currentCell.y) * frequency) > (frequency * frequency)){
                        std::cerr << "noise.cpp: worleyNoise: index" << (currentCell.x + (currentCell.y) * frequency) << ':' << frequency << "out of bounds" << std::endl;
                        exit(1);
                    }
                    // calculate the distance vector (allowing for boundary looping
                    distVector = {cellsize.x * celli + grid[currentCell.x + (currentCell.y) * frequency].x - (i % cellsize.x),\
                                  cellsize.y * cellj + grid[currentCell.x + (currentCell.y) * frequency].y - (j % cellsize.y)};

                    //distVector = {grid[currentCell.x + (currentCell.y * frequency)].x - i + cellsize.x * currentCell.x,\
                        grid[currentCell.x + (currentCell.y * frequency)].y - j + cellsize.y * currentCell.y};
                    // calculate the closest distance
                    distance = pow(pow(distVector.x, 2) + pow(distVector.y, 2), 0.5);
                    if (distance < closestDist || closestDist < 0) {
                        closestDist = distance;
                    }
                }
            }

            // set the pixel value
            noise[j * dims.y + i] = flinearBlend(closestDist, {0, (float)MAX_COEFF*cellsize.x}, {0, 255});

            closestDist = -1;
        }
    }

    delete[] grid;

    return noise;
}

void Noise::convolve(int2 dims, unsigned char* image, char filter){
    int filterSize = 3;
    float* filterArray = new float[filterSize * filterSize];

    // set the filter
    if (filter == 'G'){
        float tempArray[] = {0.25, 0.5, 0.25, 0.5, 1, 0.5, 0.25, 0.5, 0.25};
    for (int i = 0; i < filterSize * filterSize; i++) {
            filterArray[i] = tempArray[i];
        }
    } else {
        for (int i = 0; i < filterSize * filterSize; i++) {
            //filterArray[i] = {1, 1, 1, 1, 1, 1, 1, 1, 1}[i];
        }
    }

    // sum variable
    int sum = 0;

    // filter amplitude based on size
    int filAmp = (filterSize - 1) / 2;
    // temporary variable to calculate which image pixel to use                                   
    int imagePix = 0;

    // convolve the image with the filter
    for (int i = 0; i < dims.x * dims.y; i++){

        // go through each pixel and add up the sums of each surrounding pixel
        for (int filterX = -filAmp; filterX <= filAmp; filterX++){
            for (int filterY = -filAmp; filterY <= filAmp; filterY++){
                imagePix = i + filterX + filterY * dims.x;

                if (imagePix < 0 || imagePix >= dims.x * dims.y || i % dims.x + filterX < 0 || i % dims.x + filterX >= dims.x || i / dims.x + filterY < 0 || i / dims.x + filterY >= dims.y){
                    sum += 0;
                } else {
                    sum += image[imagePix] * filterArray[filterX + filterY * filterSize];
                }
            }
        }

        sum = 0;
    }
}


unsigned char* Noise::gridTransform(unsigned char* domain, int multiplier){
    int2 dimensions = dims;
    int2* mappedPos = new int2[dimensions.x * dimensions.y];
    int totalXdist = 0;
    int* totalYdist = new int[dimensions.x];

    for(int j = 0; j < dimensions.y; j++){
        for(int i = 0; i < dimensions.x; i++){
            int lowerXBound = ((i-1) + dimensions.x) % dimensions.x;
            int lowerYBound = ((j-1) + dimensions.y) % dimensions.y;
            int upperXBound = (i+1) % dimensions.x;
            int upperYBound = (j+1) % dimensions.y;

            int xdist = totalXdist + pow(domain[i + j* dimensions.x], (float)1/multiplier);
            int ydist = totalYdist[i] + pow(domain[i + j* dimensions.x], (float)1/multiplier);

            totalYdist[i] = ydist;
            totalXdist = xdist;

            mappedPos[i + j*dimensions.x] = {xdist, ydist}; //{i, j};
        }
        totalXdist = 0;
    }

    unsigned char* gridImg = new unsigned char[dimensions.x * dimensions.y];

    for(int l = 0; l < dimensions.x; l++){
        for(int k = 0; k < dimensions.y; k++){
            int i = mappedPos[l + k*dimensions.x].x;
            int j = mappedPos[l + k*dimensions.x].y;
            gridImg[l + k*dimensions.x] = (cos((i)/(pi*10))*cos((j)/(pi*10)) >= 0) ? 255 : 0;
        }
    }

    return gridImg;
}


