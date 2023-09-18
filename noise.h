#pragma once
#include <math.h>

class perlin {
	public:
		static float getPerlin(float x, float y);
		static float getPerlin(float x, float y, int octaves, float persistence);
		static float getNormal(float omax, float omin, float max, float min, float value);
		static float getNormalNoise(float x, float z);
		static float cosine_interpolate(float a, float b, float x);
		static float smooth_noise_2D(float x, float y);
		static float interpolated_noise(float x, float y);
		static float noise(int x, int y);
		
		static float octaves;
		static float persistence;
		static float lacunarity;
};

#define PI 3.1415927

float perlin::noise(int x, int y) {
    int n = x + y * 57;
    n = (n<<13) ^ n;
    return (1.0 - ( (n * ((n * n * 15731) + 789221) +  1376312589) & 0x7fffffff) / 1073741824.0);
}

float perlin::cosine_interpolate(float a, float b, float x) {
    float ft = x * PI;
    float f = (1 - cos(ft)) * 0.5;
    float result =  a*(1-f) + b*f;
    return result;
}

float perlin::smooth_noise_2D(float x, float y) {  
    float corners = ( noise(x-1, y-1)+noise(x+1, y-1)+noise(x-1, y+1)+noise(x+1, y+1) ) / 16;
    float sides   = ( noise(x-1, y)  +noise(x+1, y)  +noise(x, y-1)  +noise(x, y+1) ) /  8;
    float center  =  noise(x, y) / 4;

    return corners + sides + center;
}

float perlin::interpolated_noise(float x, float y) {
    int x_whole = (int) x;
    float x_frac = x - x_whole;

    int y_whole = (int) y;
    float y_frac = y - y_whole;

    float v1 = smooth_noise_2D(x_whole, y_whole); 
    float v2 = smooth_noise_2D(x_whole, y_whole+1); 
    float v3 = smooth_noise_2D(x_whole+1, y_whole); 
    float v4 = smooth_noise_2D(x_whole+1, y_whole+1); 

    float i1 = cosine_interpolate(v1,v3,x_frac);
    float i2 = cosine_interpolate(v2,v4,x_frac);

    return cosine_interpolate(i1, i2, y_frac);
}

float perlin::octaves = 8.0f;
float perlin::persistence = 0.5f;
float perlin::lacunarity = 1.0f;

float perlin::getPerlin(float x, float y, int octaves, float persistence) {
    float total = 0;

    for(int i=0; i<octaves-1; i++)
    {
        float frequency = pow(2,i);
        float amplitude = pow(persistence,i);
        total = total + interpolated_noise(x * frequency, y * frequency) * amplitude;
    }
    return total;
}

float perlin::getPerlin(float x, float y) {
    float total = 0;
	//x += int(1 << 20);
	//y += int(1 << 20);

    for(int i=0; i<octaves-1; i++)
    {
        float frequency = pow(2,i);
        float amplitude = pow(persistence,i);
        total = total + interpolated_noise(x * frequency, y * frequency) * amplitude;
    }
    return total;
}
 
float perlin::getNormal(float omax, float omin, float max, float min, float value) {
	return (max - min) / (omax - omin) * (value - omax) + max;
}

float perlin::getNormalNoise(float x, float z) {
	float noise = perlin::getPerlin(x,z);
	return perlin::getNormal(1,-1,1,0,noise);
}