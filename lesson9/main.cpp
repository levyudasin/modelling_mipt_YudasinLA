#include <fstream>
#include <string>
#include <iostream>
#include <cmath>
#include <array>
#include <vector>
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

template<size_t N>
using State = std::array<double, N>;

template<size_t N>
using TState = std::pair<double, State<N>>;

// N - number of variables
template<size_t N>
class Model{
public:
    virtual State<N> RightPart(State<N> state, double t) = 0;
    State<N> startValues;
};

template<size_t N>
class Method
{
protected:
    State<N> virtual Step(State<N> state, Model<N> &m, double t, double dt) = 0;
public:
    std::vector<TState<N>> Solve(Model<N> &m, double dt, double simTime){
        std::vector<TState<N>> answer = {{0,m.startValues}};
        double t = 0;
        for(int i = 0; i < int(simTime / dt); i++){
            t += dt;
            answer.emplace_back(t, Step(answer.back().second, m, t, dt));
        }
        return answer;
    }

    virtual ~Method() = default;
};

template <size_t N>
class EulerMethod : public Method<N>{
protected:
    State<N> Step (State<N> state, Model<N> &m, double t, double dt) override {
        State<N> answer = state;
        for(int i = 0; i < N; i++){
            answer[i] = state[i] + dt * m.RightPart(state, t)[i];
        }
        return answer;
    }
};

template <size_t N>
class HeunMethod : public Method<N> {
protected:
    State<N> Step(State<N> state, Model<N>& m, double t, double dt) override {
        State<N> answer = state;
        State<N> pre_answer = state;
        for (int i = 0; i < N; i++) {
            pre_answer[i] = state[i] + dt * m.RightPart(state, t)[i];
            answer[i] = state[i] + dt / 2 * (m.RightPart(state, t)[i] + m.RightPart(pre_answer, t + dt)[i]);
        }
        return answer;
    }
};

class MathPendulum : public Model<2>
{
private:

    double v_t(double t, double x, double v) {
        return -omega * omega * x;
    }

    double x_t(double t, double x, double v) {
        return v;
    }
public:
    MathPendulum(double x, double v){
        startValues = {x, v};
    }
    
    State<2> RightPart(State<2> state, double t) override {
        return State<2>({ x_t(t, state[0], state[1]), v_t(t, state[0], state[1])});
    }
};



// void stepRK(double x_i, double v_i, double dt){
//     double vk1 = v_t(t, x_i, v_i);
//     double xk1 = x_t(t, x_i, v_i);

//     double vk2 = v_t(t + dt / 2, x_i + dt * xk1 / 2, v_i + dt * vk1 / 2);
//     double xk2 = x_t(t + dt / 2, x_i + dt * xk1 / 2, v_i + dt * vk1 / 2);

//     double xk3 = x_t(t + dt / 2, x_i + dt * xk2 / 2, v_i + dt * vk2 / 2);
//     double vk3 = v_t(t + dt / 2, x_i + dt * xk2 / 2, v_i + dt * vk2 / 2);

//     double vk4 = v_t(t + dt, x_i + dt * xk3, v_i + dt * vk3);
//     double xk4 = x_t(t + dt, x_i + dt * xk3, v_i + dt * vk3);

//     v = v_i + dt / 6 * (vk1 + 2 * vk2 + 2 * vk3 + vk4);
//     x = x_i + dt / 6 * (xk1 + 2 * xk2 + 2 * xk3 + xk4);
//     t += dt;
// }

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


    MathPendulum p(x_0, v_0);
    std::vector<TState<2>> solution;

    if(method == "E"){
        EulerMethod<2> method;
        solution = method.Solve(p, dt, simTime);
    }

    else if(method == "H"){
        HeunMethod<2> method;
        solution = method.Solve(p, dt, simTime);
    }

    // else if(method == "RK"){
    //     EulerMethod<2> method;
    //     solution = method.Solve(p, dt, simTime);
    // }
    
    for (auto &&ts : solution)
    {
        outFile << ts.first << "," << ts.second[0] << "," << ts.second[1] << std::endl;
    }
    



    // for(int i = 0; i < int(simTime / dt); i++){
    //     outFile << t << "," << x << "," << v << std::endl;
    //     //if(method == "E") stepEuler(x,v,dt);
    //     //if(method == "H") stepHeun(x,v,dt);
    //     if(method == "RK") stepRK(x,v,dt);
    // }

    return 0;
}