#include <iostream>
#include <algorithm>
#include <chrono>
#include <vector>

#include "src/math.h"
#include "src/vectors.h"
#include "src/noise.h"
#include "src/error.h"
#include "src/readPNG/include/readPNG.hpp"
#include "src/image.h"


using namespace std;

// Dimensions of the image
// constexpr int2 dims = {3024, 1964};
int2 dims = {3024, 1964};

// DRAW STUFF
void drawSquare(int2 corner, int2 size, int3 color, int* red, int* green, int* blue){

  for (int i = corner.y; i < corner.y + size.y; i++){
    for (int j = corner.x; j < corner.x + size.x; j++){
      red[i * dims.x + j] = color.x;
      green[i * dims.x + j] = color.y;
      blue[i * dims.x + j] = color.z;
    }
  }
}
void drawSquareT(int2 center, int2 size, int* textureRed, int* textureGreen, int* textureBlue, int* red, int* green, int* blue){
  // REQUIRES: texture to be the size of the square
  int2 boundsX = {max(0, center.x - size.x/2), min(dims.x, center.x + size.x/2)};
  int2 boundsY = {max(0, center.y - size.y/2), min(dims.y, center.y + size.y/2)};
  for (int j = boundsY.x; j < boundsY.y; j++){
    for (int i = boundsX.x; i < boundsX.y; i++){
      red[j * dims.x + i] = textureRed[j * dims.x + i];
      green[j * dims.x + i] = textureGreen[j * dims.x + i];
      blue[j * dims.x + i] = textureBlue[j * dims.x + i];
    }
  }
}


// draw a circle based on a function to determine the color
void drawCircle(int2 center, int radius, color (*colorFunc)(int2 pos), color* outputRGB){
  int2 boundsX = {max(0, center.x - radius), min(dims.x, center.x + radius)};
  int2 boundsY = {max(0, center.y - radius), min(dims.y, center.y + radius)};
  for (int j = boundsY.x; j < boundsY.y; j++){
    for (int i = boundsX.x; i < boundsX.y; i++){
      if ((int)dist({i, j}, center) <= radius){
        outputRGB[j * dims.x + i] = colorFunc({i, j});
      }
    }
  }
}
// draw a circle based on a texture array
void drawCircleT(int2 center, int radius, color* texture, color* outputRGB){
  // textures are just aligned to the corner
  // REQUIRES: texture to be the size of the square with inscribed circle
  int2 boundsX = {max(0, center.x - radius), min(dims.x, center.x + radius)};
  int2 boundsY = {max(0, center.y - radius), min(dims.y, center.y + radius)};
  for (int j = boundsY.x; j < boundsY.y; j++){
    for (int i = boundsX.x; i < boundsX.y; i++){
      if ((int)dist({i, j}, center) <= radius){
        outputRGB[j * dims.x + i] = texture[j * dims.x + i];
      }
    }
  }
}

void drawSphere(int2 center, int radius, color (*sphereMapping)(int rad, int2 pos, color* noiseTexture, int2 texSize), color* noiseTexture, int2 texSize, color* outputRGB){
  int2 boundsX = {max(0, center.x - radius), min(dims.x, center.x + radius)};
  int2 boundsY = {max(0, center.y - radius), min(dims.y, center.y + radius)};
  for (int j = boundsY.x; j < boundsY.y; j++){
    for (int i = boundsX.x; i < boundsX.y; i++){
      if ((int)dist({i, j}, center) <= radius){
        outputRGB[j * dims.x + i] = sphereMapping(radius, {(i-center.x), (j-center.y)}, noiseTexture, texSize);
      }
    } 
  }
}

// layer weight calculations
double noiseWeightCalc(int nLayers, int layer){
  return (pow(2, layer-nLayers));
}
double noiseWeightCalc2(int nLayers, int layer){
  return 1/(pow(layer, 2));
}
double worleyWeightCalc(int nLayers, int layer){
  return (pow(2, layer-nLayers));
}

// color channel calculations
color colorFunc1(int2 pos){
  unsigned char shade = linearBlend(dist(pos, {dims.x/3, dims.y/4}), {0, (int)(dims.y*pow(2, 0.5))}, {0, 255});
  return {shade, shade, shade};
}
color riemmanMapping(int r, int2 pos, color* noiseTexture, int2 texSize){
  float z = sqrt(r*r - pos.x*pos.x - pos.y*pos.y);
  float diff = max((double)(r-pos.y), (double)1);
  float2 mappedPos = {static_cast<float>(pos.x * r) / diff, static_cast<float>(z * r) / diff};

  // normalize the mapped position
  mappedPos.x = (mappedPos.x >= 0) ? mappedPos.x : ((int)mappedPos.x % texSize.x) + texSize.x;
  mappedPos.y = (mappedPos.y >= 0) ? mappedPos.y : ((int)mappedPos.y % texSize.y) + texSize.y;

  int noise_map = ((int)(mappedPos.x) % texSize.x) + texSize.x * ((int)(mappedPos.y) % texSize.y);
  return noiseTexture[noise_map];
}
color geodesicMapping(int r, int2 pos, color* noiseTexture, int2 texSize){
  float z = sqrt(r*r-pos.x*pos.x - pos.y*pos.y);
  int2 r_div = {static_cast<int>(sqrt(pos.x*pos.x + z*z)), static_cast<int>(sqrt(pos.y*pos.y + z*z))};
  double2 angles = {asin((double)pos.x/r_div.x), atan((double)pos.y/r)};
  int2 mappedPos = {static_cast<int>((angles.x + pi/2.00) * texSize.x / pi), static_cast<int>((angles.y + pi/2.00) * texSize.y / pi)};

  return noiseTexture[((int)mappedPos.x % texSize.x) + texSize.x * ((int)mappedPos.y % texSize.y)];
}
color topdownMapping(int r, int2 pos, color* noiseTexture, int2 texSize){
  float z = sqrt(r*r-pos.x*pos.x - pos.y*pos.y);
  int2 r_div = {static_cast<int>(sqrt(pos.x*pos.x + z*z)), static_cast<int>(sqrt(pos.y*pos.y + z*z))};
  double2 angles = {asin((double)z/r), atan((double)pos.y/pos.x)};
  int2 mappedPos = {static_cast<int>((angles.x + pi/2.00) * texSize.x / pi), static_cast<int>((angles.y + pi/2.00) * texSize.y / pi)};

  return noiseTexture[((int)mappedPos.x % texSize.x) + texSize.x * ((int)mappedPos.y % texSize.y)];
}

color latlongMapping(int r, int2 pos, color* noiseTexture, int2 texSize){
  float z = sqrt(r*r-pos.x*pos.x - pos.y*pos.y);
  int2 r_div = {static_cast<int>(sqrt(pos.x*pos.x + z*z)), static_cast<int>(sqrt(pos.y*pos.y + z*z))};
  double2 angles = {asin((double)pos.x/r_div.x), atan((double)pos.y/r_div.y)};
  int2 mappedPos = {static_cast<int>((angles.x + pi/2.00) * texSize.x / pi), static_cast<int>((angles.y + pi/2.00) * texSize.y / pi)};

  return noiseTexture[((int)mappedPos.x % texSize.x) + texSize.x * ((int)mappedPos.y % texSize.y)];
}

color lightMapping(int r, int2 pos, color* noiseTexture, int2 texSize){
  float z = sqrt(pos.x*pos.x + pos.y*pos.y);
  int3 pos3 = {pos.x, pos.y, static_cast<int>(z)};

  float3 vec1 = {static_cast<float>(pos3.x - 0), static_cast<float>(pos3.y - 0), static_cast<float>(pos3.z - 500)};
  float3 vec2 = {static_cast<float>(pos3.x - 300), static_cast<float>(pos3.y + 600), static_cast<float>(pos3.z - 400)};
  unsigned char val = linearBlend(min(dot(normalize(vec1), normalize(vec2)), (float)255), {1, -1}, {0, 255});
  return {val, val, val};
}

int graphicsMain(void){

  // color* img = loadImgFromPng("output_image.png");
  std::vector<color> img = loadImgFromPng("peyto.png", dims);
  if(!img.size()){
    ERROR("image read failed");
    

    return 1;
  }
  
  // create arrays for each color channel
  color* rgb = new color[dims.x * dims.y];

  for (int i = 0; i < dims.x * dims.y; i++){
    if(img[i].b > 225) rgb[i] = {0,0,0};
    else{
      rgb[i].r = img[i].r;
      rgb[i].g = img[i].g;
      rgb[i].b = img[i].b;
    }
  }


#if 0
  // create the object to create noise
  Noise noise(dims);

  /* ---------CREATE TEXTURE LAYERS--------- */
  // perlin texture layer
  unsigned char* perlinNoise = noise.genPerlin(5, {10, 5}, noiseWeightCalc, false);
  unsigned char* perlinNoise2 = noise.genPerlin(3, {50, 25}, noiseWeightCalc, false);
  unsigned char* perlinNoise3 = noise.genPerlin(1, {500, 500}, noiseWeightCalc, false);
  //worley Noise layer
  // unsigned char* worleyNoise = noise.genWorley(1, 2048, {4, 4}, worleyWeightCalc, false);
  //noise.genCells(1600, 6, 0, 0);

  /* ---------ADD TEXTURES TO IMAGE (MASKING)--------- */
  // color* temp = new color[dims.x * dims.y];

  // for(int i = 0; i < 100; i++){
  //     for(int j = 0; j < 100; j++){
  //         unsigned char val = (127.5*abs(sin((i)/(pi*5))+cos((j)/(pi*5))) + 127.5);
  //         // unsigned char val = img[i*100 + j].r;
  //         temp[100 * i + j] = {val, val, val};
  //     }
  // }

  for (int i = 0; i < dims.x*dims.y; i++){
    if(img[i].b > 225){
      rgb[i].r = perlinNoise[i];
      rgb[i].g = perlinNoise2[i];
      rgb[i].b = perlinNoise3[i];
    }
    else{
      rgb[i].r = img[i].r;
      rgb[i].g = img[i].g;
      rgb[i].b = img[i].b;
    }
  }

  delete [] perlinNoise;
  delete [] perlinNoise2;
  delete [] perlinNoise3;

  /* ---------DRAW STUFF--------- */
  // drawSphere({dims.x/4*3, dims.y/2}, 300, riemmanMapping, temp, dims, rgb);

  //drawSphere({dims.x/2, dims.y/4}, 300, topdownMapping, temp, dims, rgb);

  //drawSphere({dims.x/2, dims.y/4*3}, 300, latlongMapping, temp, dims, rgb);

  //drawSphere({dims.x/4, dims.y/2}, 300, geodesicMapping, temp, dims, rgb);


#endif

  int* red = new int[dims.x * dims.y];
  int* green = new int[dims.x * dims.y];
  int* blue = new int[dims.x * dims.y];

  for (int i = 0; i < dims.x * dims.y; i++){
    red[i] = rgb[i].r;
    green[i] = rgb[i].g;
    blue[i] = rgb[i].b;
  }

  std::string outputFilename = "imgR.csv";
  writeImgToCSV(red, outputFilename, dims);
  outputFilename = "imgG.csv";
  writeImgToCSV(green, outputFilename, dims);
  outputFilename = "imgB.csv";
  writeImgToCSV(blue, outputFilename, dims);

  delete[] red;
  delete[] green;
  delete[] blue;
  delete[] rgb;

  return 0;
}

int timeFunc(int (*timed_func)(void)){
  auto start = std::chrono::high_resolution_clock::now();

  timed_func();

  auto end = std::chrono::high_resolution_clock::now();

  return (end - start).count(); 
}

int main(){
  int time = timeFunc(graphicsMain);

  std::cout << "Time to run: " << time/1000 << " milliseconds" << std::endl;
}

// TODO:
// make a discrete convolution function
// write a image function to load images from png's
// convert all the color arrays to color type

// DONE:
// create a layer system and a way to blend between layers (opaque, mean value, weighted interpolationb/w layers)


// noise references  
// https://www.youtube.com/watch?v=gsJHzBTPG0Y
