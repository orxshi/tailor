#ifndef LOGGER_H
#define	LOGGER_H

#include <fstream>
#include <string>

class Logger
{
    std::ofstream out;
    public:

    Logger(std::string file_name);

    void record(std::string s, int i);
};

#endif
