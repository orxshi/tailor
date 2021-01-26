#include "logger.h"

Logger::Logger(std::string file_name)
{
    out.open(file_name);
}

void Logger::record(std::string s, int i)
{
    out << s << " - " << i << std::endl;
}
