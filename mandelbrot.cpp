#include <iostream>
#include <fstream>
#include <math.h>

#include <Eigen/Dense>

using Eigen::Matrix2d;
using Eigen::Vector2d;
using Eigen::Vector3d;

double mod(std::complex<double> z){
    double real = static_cast<double>(z.real());
    double imag = static_cast<double>(z.imag());
    double mod_squared = std::pow(real, 2) + std::pow(imag, 2);
    double modulus = sqrt(mod_squared);
    return modulus;
}

bool in_mandelbrot(std::complex<double> c){
    std::complex<double> z_n_1;
    std::complex<double> z_n = 0;
    bool in_set = true;

    for(int i; i < 10; ++i){
        z_n_1 = std::pow(z_n, 2) + c;
        z_n = z_n_1;
        double modulus = mod(z_n);
        if(modulus > 2){
            in_set = false;
        }
    }
    return in_set;
}


void mandelbrot_data(){
/*

current idea is to make a matrix roughly from (tl, tr, br, bl) 
(-2+1.2i, 1.2+1.2i, 1.2-1.2i, -2-1.2i) which is very dense (lots
of points close together) and then iterate through all of them and if check if
they're in the mandelbrot set using "in_mandelbrot" function. Then get this 
data written to data.txt file which can be used to plot in python. This should
give a pretty decent visual of the mandelbrot set.

*/
}

int main(){
    std::ofstream myfile;
    myfile.open("data.txt");
    myfile << in_mandelbrot(1);
    myfile.close();
    std::cout << in_mandelbrot(1);
    return 0;
}

