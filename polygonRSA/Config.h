//--------------------------------------------------------------------------------------------
// Simple configuration file parser
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _CONFIG_H
    #define _CONFIG_H


#include <map>
#include <stdexcept>
#include <iosfwd>
#include <vector>
#include <algorithm>


// Parse exception class
class ConfigParseException : public std::runtime_error
{
public:
    explicit ConfigParseException(const std::string & _what) : std::runtime_error(_what)
    {}
};


// No field exception class
class ConfigNoFieldException : public std::runtime_error
{
public:
    explicit ConfigNoFieldException(const std::string & _what) : std::runtime_error(_what)
    {}
};


/**
 * @brief A key=value config file parser
 *
 * File format:
 * \code
 * key1=value1
 * # standalone comment
 * key2=value2 # end-line comment;
 *
 * # empty lines are omitted
 * key3 = value3 # whitespace is trimmed
 * # key3 = value3 - duplicate fields are forbidden
 * \endcode
 */
class Config
{
private:
    std::map<std::string, std::string>  fieldMap;
    std::vector<std::string>            keys;

    struct Field {
        std::string key;
        std::string value;
    };

    static void stripComment(std::string &line);
    static Field splitField(const std::string &line, char delim, std::size_t line_num);

    Config() = default;

public:

    /**
     * @brief Parses given stream.
     *
     * Format:
     * \code
     * key1=value1
     * # standalone comment
     * key2=value2 # end-line comment;
     *
     * # empty lines are omitted
     * key3 = value3 # whitespace is trimmed
     * # key3 = value3 - duplicate fields are forbidden
     * \endcode
     *
     * @param in stream to parse from
     * @param delim delimiter for key, value; defaults to '='
     * @throws ConfigParseException on parse error (no delimiter of a duplicate field)
     * @throws std::invalid_argument when delim = '#' (comment)
     * @return Config object to be deleted manualy after use
     */
    static Config parse(std::istream &in, char delim = '=', bool allowRedefinition = false);

    bool hasParam(const std::string &field) const { return std::find(keys.begin(), keys.end(), field) != keys.end(); };
    std::string getString(const std::string &field) const;
    int getInt(const std::string &field) const;
    unsigned long getUnsignedLong(const std::string &field) const;
    double getDouble(const std::string &field) const;
    float getFloat(const std::string &field) const;

    /**
     * Returns keys in a config, preserving order from an input
     * @return keys in a config, preserving order from an input
     */
    std::vector<std::string> getKeys() const { return this->keys; };
};

    
#endif // _CONFIG_H
