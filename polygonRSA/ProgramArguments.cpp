//
// Created by pkua on 07.04.19.
//

#include <fstream>
#include <iostream>

#include "ProgramArguments.h"

ProgramArguments::ProgramArguments(int argc, char **argv) {
    Expects(argc >= 1);     // valid program arguments have at least 1 arg - cmd
    cmd = argv[0];

    // Now valid program arguments, but we assume we always have mode indicator
    if (argc < 2)
        throw InvalidArgumentsException("Usage: " + cmd + " <mode> (additional parameters)");

    std::stringstream input;
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-f") {
            if (++i == argc)
                throw InvalidArgumentsException("Expected input file after -f. Aborting.");

            std::string inputFileName = argv[i];
            std::ifstream inputFile(inputFileName);
            if (!inputFile)
                throw InvalidArgumentsException("Cannot open input file " + inputFileName + "to read. Aborting.");

            input << inputFile.rdbuf() << std::endl;
        } else if (arg == "-i") {
            input << std::cin.rdbuf() << std::endl;
        } else if (arg.length() > 1 && arg[0] == '-') {
            input << arg.substr(1) << std::endl;
        } else {
            positionalArguments.push_back(arg);
        }
    }

    if (positionalArguments.empty())
        throw std::runtime_error("No mode parameter. Aborting.");
    mode = positionalArguments[0];
    positionalArguments.erase(positionalArguments.begin());

    parameters = Parameters(input);
}

std::string ProgramArguments::formatUsage(const std::string &additionalArgs) const {
    return "Usage: " + cmd + " " + mode + " " + additionalArgs;
}
