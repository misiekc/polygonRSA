//
// Created by pkua on 07.04.19.
//

#ifndef RSA3D_PROGRAMARGUMENTS_H
#define RSA3D_PROGRAMARGUMENTS_H

#include <string>
#include <vector>

#include "Parameters.h"

class InvalidArgumentsException : public std::runtime_error {
public:
    explicit InvalidArgumentsException(const std::string &msg) : runtime_error{msg} {}
};

class ProgramArguments {
private:
    Parameters parameters;
    std::string cmd;
    std::string mode;
    std::vector<std::string> positionalArguments;

public:
    ProgramArguments(int argc, char **argv);

    const std::string &getCmd() const { return cmd; }
    const Parameters &getParameters() const { return parameters; }
    const std::string &getMode() const { return mode; }
    const std::vector<std::string> &getPositionalArguments() const { return positionalArguments; }

    std::string formatUsage(const std::string &additionalArgs) const;
};

#endif //RSA3D_PROGRAMARGUMENTS_H
