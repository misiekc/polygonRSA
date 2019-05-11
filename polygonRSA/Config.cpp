//--------------------------------------------------------------------------------------------
// Simple configuration file parser
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include "Config.h"
#include "Utils.h"

#include <memory>


Config Config::parse(std::istream & in, char delim, bool allowRedefinition) {
    if (delim == '#')
        throw std::invalid_argument("delim == #");

    Config result;
    std::size_t line_num = 0;
    std::string line;
    while (std::getline(in, line)) {
        line_num++;
        stripComment(line);
        if (line.empty())
            continue;

        auto field = splitField(line, delim, line_num);
        if (result.hasParam(field.key) && !allowRedefinition)
            throw ConfigParseException("Redefinition of field \"" + field.key + "\" in line " + std::to_string(line_num));

        result.fieldMap[field.key] = field.value;
        if (!result.hasParam(field.key))
            result.keys.push_back(field.key);
    }

    return result;
}

Config::Field Config::splitField(const std::string &line, char delim, std::size_t line_num) {
    Field keyValue;
    std::size_t pos = line.find(delim);
    if (pos == std::string::npos)
        throw ConfigParseException("No '" + std::string(1, delim) + "' sign in line " + std::to_string(line_num));

    keyValue.key = line.substr(0, pos);
    trim(keyValue.key);
    keyValue.value = (pos == line.length() - 1) ? "" : line.substr(pos + 1);
    trim(keyValue.value);
    return keyValue;
}

void Config::stripComment(std::string &line) {
    std::size_t pos = line.find('#');
    if (pos != std::string::npos)
        line.erase(pos);
}

std::string Config::getString(const std::string & field) const {
    auto iter = this->fieldMap.find(field);
    if (iter == this->fieldMap.end())
        throw ConfigNoFieldException("No \"" + field + "\" field in config");
    return (*iter).second;
}

int Config::getInt(const std::string & field) const {
    return std::stoi(this->getString(field));
}

unsigned long Config::getUnsignedLong(const std::string & field) const {
    auto str = this->getString(field);
    if (std::stoi(str) < 0)
        throw std::invalid_argument("unsigned long field negative");
    return std::stoul(str);
}

double Config::getDouble(const std::string & field) const {
    return std::stod(this->getString(field));
}

float Config::getFloat(const std::string & field) const
{
    return std::stof(this->getString(field));
}

