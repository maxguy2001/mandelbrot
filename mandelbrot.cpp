#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <complex>


/*
@brief returns the modulus of a complex number
*/
double mod(std::complex<double> z){
    double real = z.real();
    double imag = z.imag();
    double mod_squared = std::pow(real, 2) + std::pow(imag, 2);
    double modulus = sqrt(mod_squared);
    return modulus;
}

/*
@brief returns a boolean value of if a complex number is in the mandelbrot set.
*/
bool inMandelbrot(std::complex<double> c){
    std::complex<double> z_n_1;
    std::complex<double> z_n = 0;
    bool in_set = false;
    int modulus;

    for(size_t i=0; i < 200; ++i){
        z_n_1 = std::pow(z_n, 2) + c;
        z_n = z_n_1;
        modulus = mod(z_n);

        if(modulus > 5){
            break;
        }
    }
    if (modulus <= 2){
        in_set = true;
    }

    return in_set;
}

/*
@brief returns a progress bar in terminal
*/
std::string progressBar(double progress, double total){
    int num_progress_markers = static_cast<int>((progress/total)*50);
    int num_not_progress_markers = 50 - num_progress_markers;
    std::string progress_bar = "[";

    for(size_t i = 0; i < num_progress_markers; ++i){
        progress_bar.append("#");
    }
    for(size_t i = 0; i < num_not_progress_markers; ++i){
        progress_bar.append("."); 
    }
    progress_bar.append("]");
    std::cout << progress_bar << "\r";

    return "";
}

/*
@brief takes n by m evenly spaced points in an area of the complex plane
and wrties them along with a boolean mask of if they are in the mandelbrot
set to files.
*/
void mandelbrotData(int points_along_real, int points_along_imag){

//define area we want to focus on and step size
    std::complex<double> top_left = -2+1.2j;
    std::complex<double> bottom_right = 1.2-1.2j;


    std::complex<double> delta_real = std::abs(top_left.real() - bottom_right.real()) /
                        static_cast<std::complex<double>>(points_along_real/2);
    double delta_imag_real = std::abs(top_left.imag() - bottom_right.imag()) /
                        static_cast<double>(points_along_imag/2);
    std::complex<double> imaginary_number = 1j;
    std::complex<double> delta_imag = delta_imag_real*imaginary_number;

    //vectors for storing datapoints and if they are in the mandelbrot set
    std::vector<std::complex<double>> datapoints;    
    std::complex<double> current_point = top_left;

    for(size_t re = 1; re < points_along_real; ++re){
        std::complex<double> real_multiplier = static_cast<std::complex<double>>(re);
        current_point = top_left + real_multiplier*delta_real;

        if(re % 10 == 0){
            std::string progress_bar = progressBar(re, points_along_real);
            std::cout << progress_bar << "\r";
        }

        for(size_t im = 1; im < points_along_imag; ++im){

            // finding points and adding to points vector
            std::complex<double> imag_multiplier = static_cast<std::complex<double>>(im);
            std::complex<double> imag_current_point = current_point - delta_imag*imag_multiplier;

            std::ofstream datfile;
            datfile.open("data_matrix.txt", std::ios::app);
            datfile << imag_current_point << "\n";
            datfile.close();

            bool in_set = inMandelbrot(imag_current_point);
            std::ofstream boolfile;
            boolfile.open("bool_matrix.txt", std::ios::app);
            boolfile << in_set << "\n";
            boolfile.close();

        }
    }
}


int main(){
    mandelbrotData(1e5, 1e5);
    return 0;
}

