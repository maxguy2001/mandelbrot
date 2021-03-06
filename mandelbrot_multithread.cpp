#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <complex>
#include <string>
#include <thread>
#include <pthread.h>


/*
@brief returns the modulus of a complex number
*/
double mod(const std::complex<double> z){
    const double real = z.real();
    const double imag = z.imag();
    const double mod_squared = std::pow(real, 2) + std::pow(imag, 2);
    const double modulus = sqrt(mod_squared);
    return modulus;
}

/*
@brief returns a boolean value of if a complex number is in the mandelbrot set.
*/
bool inMandelbrot(const std::complex<double> c){
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
@brief takes n by m evenly spaced points in an area of the complex plane
and wrties them along with a boolean mask of if they are in the mandelbrot
set to files.
*/
void MandelbrotData_t1(const int points_along_real, const int points_along_imag){
    const std::complex<double> top_left = -1.5+1.2j;
    const std::complex<double> bottom_right = -0.45-0.05j;
    
    const std::complex<double> delta_real = std::abs(top_left.real() - bottom_right.real()) /
                        static_cast<std::complex<double>>(points_along_real);
    const double delta_imag_real = std::abs(top_left.imag() - bottom_right.imag()) /
                        static_cast<double>(points_along_imag);
    const std::complex<double> imaginary_number = 1j;
    const std::complex<double> delta_imag = delta_imag_real*imaginary_number;

    std::complex<double> current_point = top_left;

    for(size_t re = 1; re < points_along_real; ++re){
        std::complex<double> real_multiplier = static_cast<std::complex<double>>(re);
        current_point = top_left + real_multiplier*delta_real;

        for(size_t im = 1; im < points_along_imag; ++im){

            // finding points and adding to points vector
            std::complex<double> imag_multiplier = static_cast<std::complex<double>>(im);
            std::complex<double> imag_current_point = current_point - delta_imag*imag_multiplier;

            if(inMandelbrot(imag_current_point) == true){
                std::ofstream datfile_1;
                datfile_1.open("data_t1.txt", std::ios::app);
                datfile_1 << imag_current_point << "\n";
                datfile_1.close();
            }

        }
    }
}



/*
@brief takes n by m evenly spaced points in an area of the complex plane
and wrties them along with a boolean mask of if they are in the mandelbrot
set to files.
*/
void MandelbrotData_t2(const int points_along_real, const int points_along_imag){
    const std::complex<double> top_left = -0.55+1.2j;
    const std::complex<double> bottom_right = 0.5-0.05j;
    
    const std::complex<double> delta_real = std::abs(top_left.real() - bottom_right.real()) /
                        static_cast<std::complex<double>>(points_along_real);
    const double delta_imag_real = std::abs(top_left.imag() - bottom_right.imag()) /
                        static_cast<double>(points_along_imag);
    const std::complex<double> imaginary_number = 1j;
    const std::complex<double> delta_imag = delta_imag_real*imaginary_number;

    std::complex<double> current_point = top_left;

    for(size_t re = 1; re < points_along_real; ++re){
        std::complex<double> real_multiplier = static_cast<std::complex<double>>(re);
        current_point = top_left + real_multiplier*delta_real;

        for(size_t im = 1; im < points_along_imag; ++im){

            // finding points and adding to points vector
            std::complex<double> imag_multiplier = static_cast<std::complex<double>>(im);
            std::complex<double> imag_current_point = current_point - delta_imag*imag_multiplier;

            if(inMandelbrot(imag_current_point) == true){
                std::ofstream datfile_2;
                datfile_2.open("data_t2.txt", std::ios::app);
                datfile_2 << imag_current_point << "\n";
                datfile_2.close();
            }

        }
    }
}


/*
@brief takes n by m evenly spaced points in an area of the complex plane
and wrties them along with a boolean mask of if they are in the mandelbrot
set to files.
*/
void MandelbrotData_t3(const int points_along_real, const int points_along_imag){
    const std::complex<double> top_left = -1.5+0.05j;
    const std::complex<double> bottom_right = -0.45-1.2j;
    
    const std::complex<double> delta_real = std::abs(top_left.real() - bottom_right.real()) /
                        static_cast<std::complex<double>>(points_along_real);
    const double delta_imag_real = std::abs(top_left.imag() - bottom_right.imag()) /
                        static_cast<double>(points_along_imag);
    const std::complex<double> imaginary_number = 1j;
    const std::complex<double> delta_imag = delta_imag_real*imaginary_number;

    std::complex<double> current_point = top_left;

    for(size_t re = 1; re < points_along_real; ++re){
        std::complex<double> real_multiplier = static_cast<std::complex<double>>(re);
        current_point = top_left + real_multiplier*delta_real;

        for(size_t im = 1; im < points_along_imag; ++im){

            // finding points and adding to points vector
            std::complex<double> imag_multiplier = static_cast<std::complex<double>>(im);
            std::complex<double> imag_current_point = current_point - delta_imag*imag_multiplier;

            if(inMandelbrot(imag_current_point) == true){
                std::ofstream datfile_3;
                datfile_3.open("data_t3.txt", std::ios::app);
                datfile_3 << imag_current_point << "\n";
                datfile_3.close();
            }

        }
    }
}



/*
@brief takes n by m evenly spaced points in an area of the complex plane
and wrties them along with a boolean mask of if they are in the mandelbrot
set to files.
*/
void MandelbrotData_t4(int points_along_real, int points_along_imag){
    const std::complex<double> top_left = -0.55+0.05j;
    const std::complex<double> bottom_right = 0.5-1.2j;
    
    const std::complex<double> delta_real = std::abs(top_left.real() - bottom_right.real()) /
                        static_cast<std::complex<double>>(points_along_real);
    const double delta_imag_real = std::abs(top_left.imag() - bottom_right.imag()) /
                        static_cast<double>(points_along_imag);
    const std::complex<double> imaginary_number = 1j;
    const std::complex<double> delta_imag = delta_imag_real*imaginary_number;

    std::complex<double> current_point = top_left;

    for(size_t re = 1; re < points_along_real; ++re){
        std::complex<double> real_multiplier = static_cast<std::complex<double>>(re);
        current_point = top_left + real_multiplier*delta_real;

        for(size_t im = 1; im < points_along_imag; ++im){

            // finding points and adding to points vector
            std::complex<double> imag_multiplier = static_cast<std::complex<double>>(im);
            std::complex<double> imag_current_point = current_point - delta_imag*imag_multiplier;

            if(inMandelbrot(imag_current_point) == true){
                std::ofstream datfile_4;
                datfile_4.open("data_t4.txt", std::ios::app);
                datfile_4 << imag_current_point << "\n";
                datfile_4.close();
            }

        }
    }
}


int main(){
    const uint16_t num_points = 2500;
    std::thread t1(MandelbrotData_t1, num_points, num_points);
    std::thread t2(MandelbrotData_t2, num_points, num_points);
    std::thread t3(MandelbrotData_t3, num_points, num_points);
    std::thread t4(MandelbrotData_t4, num_points, num_points);

    t1.join();
    t2.join();
    t3.join();
    t4.join();
    return 0;
}