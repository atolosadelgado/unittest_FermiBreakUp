#ifndef INPUTPARSER_HH
#define INPUTPARSER_HH

#include <string>

class InputParser
{
public:
    InputParser(int argc, char* argv[]);

    int massNumberA;
    int atomicNumberZ;
    double excitationEnergyMeV;
    std::string fermiBreakUpModel;
    int nEvents;

private:
    void PrintUsage(const char* programName);
};

#include <cstdlib>
#include <iostream>
#include <stdexcept>

InputParser::InputParser(int argc, char* argv[])
{
    if (argc != 6) {
        PrintUsage(argv[0]);
        throw std::runtime_error("Invalid number of arguments");
    }

    try {
        massNumberA       = std::stoi(argv[1]);
        atomicNumberZ     = std::stoi(argv[2]);
        excitationEnergyMeV = std::stod(argv[3]);
        fermiBreakUpModel = argv[4];
        nEvents           = std::stoi(argv[5]);
    }
    catch (const std::exception&) {
        PrintUsage(argv[0]);
        throw std::runtime_error("Invalid input arguments");
    }

    if (massNumberA <= 0 || atomicNumberZ < 0 ||
        atomicNumberZ > massNumberA ||
        excitationEnergyMeV < 0.0 || nEvents <= 0) {
        throw std::runtime_error("Input values outside allowed ranges");
    }
}

void InputParser::PrintUsage(const char* programName)
{
    std::cerr << "Usage: " << programName
              << " <A> <Z> <ExcitationEnergyMeV>"
              << " <FermiBreakUpModel> <nEvents>\n";
}

#endif
