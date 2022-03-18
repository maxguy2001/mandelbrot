#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>


#include <Eigen/Dense>

using Eigen::Matrix;

double mod(std::complex<double> z){
    double real = z.real();
    double imag = z.imag();
    double mod_squared = std::pow(real, 2) + std::pow(imag, 2);
    double modulus = sqrt(mod_squared);
    return modulus;
}

bool in_mandelbrot(std::complex<double> c){
    std::complex<double> z_n_1;
    std::complex<double> z_n = 0;
    bool in_set = false;
    int modulus;

    for(size_t i=0; i < 1000; ++i){
        z_n_1 = std::pow(z_n, 2) + c;
        z_n = z_n_1;
        modulus = mod(z_n);

        //std::cout << modulus << std::endl;
        if(modulus > 5){
            break;
        }
    }
    if (modulus <= 2){
        in_set = true;
    }

    return in_set;
}

void mandelbrot_data(int points_along_real, int points_along_imag){

//define area we want to focus on and step size
    std::complex<double> top_left = -2+1.2j;
    std::complex<double> bottom_right = 1.2-1.2j;
    std::complex<double> delta_real = std::abs(top_left.real() - bottom_right.real()) /
                        static_cast<std::complex<double>>(points_along_real);
    double delta_imag_real = std::abs(top_left.imag() - bottom_right.imag()) /
                        static_cast<double>(points_along_imag);
    std::complex<double> imaginary_number = 1j;
    std::complex<double> delta_imag = delta_imag_real*imaginary_number;

    std::cout << delta_imag;

//vectors for storing datapoints and if they are in the mandelbrot set
    std::vector<std::complex<double>> datapoints;    

    std::complex<double> current_point = top_left;

    for(size_t re = 1; re < points_along_real; ++re){
        std::complex<double> real_multiplier = static_cast<std::complex<double>>(re);
        current_point = top_left + real_multiplier*delta_real;

        for(size_t im = 1; im < points_along_imag; ++im){

            // finding points and adding to points vector
            std::complex<double> imag_multiplier = static_cast<std::complex<double>>(im);
            std::complex<double> imag_current_point = current_point - delta_imag*imag_multiplier;
            datapoints.push_back(imag_current_point);
        }
    }

    std::vector<bool> bool_mask;
    for(size_t i = 0; i < datapoints.size()-1; ++i){
        std::complex<double> point = datapoints.at(i);
        bool in_set = in_mandelbrot(point);
        bool_mask.push_back(in_set);
    }

    //get the sum of the boolen matrix
    int sum = 0;
    for(size_t i = 0; i < bool_mask.size()-1;++i){
        sum += bool_mask.at(i);
    }
    //print sum of bool to file
    std::ofstream sumfile;
    sumfile.open("check.txt");
    sumfile << sum;
    sumfile.close();


    //print matrix of points to file
    std::ofstream datfile;
    datfile.open("data_matrix.txt");
    for(int i = 0; i < datapoints.size()-1; ++i){
        datfile << datapoints.at(i) << "\n";
    }
    datfile.close();

    //print matrix of bool to file
    std::ofstream boolfile;
    boolfile.open("bool_matrix.txt");
    for(size_t i = 0; i < bool_mask.size()-1; ++i){
        boolfile << bool_mask.at(i) << "\n";
    }
    boolfile.close();
}

int main(){
    mandelbrot_data(100, 100);
    //in_mandelbrot(1);
    //std::cout << in_mandelbrot(1);
    return 0;
}

