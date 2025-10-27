#include <fstream>
#include <string>
#include <iostream>
#include <cmath>
#include "json.hpp"


// are taken from config file
double omega = 0;
double dt = 0;
double x_0 = 0;
double v_0 = 0;
double simTime = 0;

double t = 0;
double x = -0.5;
double v = 0;

double v_t(double t, double x, double v){
    return -omega*omega*sin(x);
}

double x_t(double t, double x, double v){
    return v;
}

void stepEuler(double x_i, double v_i, double dt){
    v = v_i + dt * v_t(t, x_i, v_i);
    x = x_i + dt * x_t(t, x_i, v_i);
    t += dt;
}

void stepHeun(double x_i, double v_i, double dt){
    double x_ = x_i + dt * x_t(t, x_i, v_i);
    double v_ = v_i + dt * v_t(t, x_i, v_i);
    v = v_i + dt / 2 * (v_t(t, x_i, v_i) + v_t(t + dt, x_i, v_));
    x = x_i + dt / 2 * (x_t(t, x_i, v) + x_t(t + dt, x_, v));
    t += dt;
}

void stepRK(double x_i, double v_i, double dt){
    double vk1 = v_t(t, x_i, v_i);
    double xk1 = x_t(t, x_i, v_i);

    double vk2 = v_t(t + dt / 2, x_i + dt * xk1 / 2, v_i + dt * vk1 / 2);
    double xk2 = x_t(t + dt / 2, x_i + dt * xk1 / 2, v_i + dt * vk1 / 2);

    double xk3 = x_t(t + dt / 2, x_i + dt * xk2 / 2, v_i + dt * vk2 / 2);
    double vk3 = v_t(t + dt / 2, x_i + dt * xk2 / 2, v_i + dt * vk2 / 2);

    double vk4 = v_t(t + dt, x_i + dt * xk3, v_i + dt * vk3);
    double xk4 = x_t(t + dt, x_i + dt * xk3, v_i + dt * vk3);

    v = v_i + dt / 6 * (vk1 + 2 * vk2 + 2 * vk3 + vk4);
    x = x_i + dt / 6 * (xk1 + 2 * xk2 + 2 * xk3 + xk4);
    t += dt;
}

int main(int argc, char** argv){
    if(argc < 2){
        std::cout << "Usage: application_name config_file_name.json" << std::endl;
        return -1;
    }
    std::string configFileName = argv[1];
    std::ifstream configFile;
    configFile.open(configFileName);
    if(!configFile.is_open()){
        std::cout << "Cannot open config file" << std::endl;
        return -1;
    }
    nlohmann::json configJson;
    configFile >> configJson;

    std::string outputFileName = configJson["output_file_name"];
    std::string method = configJson["method"];


    omega = configJson["omega"];
    dt = configJson["dt"];
    x_0 = configJson["x_0"];
    v_0 = configJson["v_0"];
    simTime = configJson["simTime"];

    std::ofstream outFile;
    outFile.open(outputFileName);
    if(!outFile.is_open()){
        std::cout << "Cannot open output file" << std::endl;
        return -1;
    }



    for(int i = 0; i < int(simTime / dt); i++){
        outFile << t << "," << x << "," << v << std::endl;
        if(method == "E") stepEuler(x,v,dt);
        if(method == "H") stepHeun(x,v,dt);
        if(method == "RK") stepRK(x,v,dt);
    }

    return 0;
}