#include <iostream>
#include <fstream>
#include "logger.h"

static Logger::Ptr_t s_logger(new Logger());
static std::ofstream s_logger_file;

void Logger::info(const std::string& text) {
    if (s_logger_file && s_logger_file.good())
        s_logger_file << "INFO: " << text << std::endl;
    else
        std::cout << "INFO: " << text << std::endl;
}
void Logger::warn(const std::string& text) {
    if (s_logger_file && s_logger_file.good())
        s_logger_file << "WARN: " << text << std::endl;
    else
        std::cout << "WARN: " << text << std::endl;
}
void Logger::critical(const std::string& text) {
    if (s_logger_file && s_logger_file.good())
        s_logger_file << "CRITICAL: " << text << std::endl;
    else
        std::cout << "CRITICAL: " << text << std::endl;
}
void Logger::debug(const std::string& text) {
    if (s_logger_file && s_logger_file.good())
        s_logger_file << "DEBUG: " <<  text << std::endl;
    else
        std::cout << "DEBUG: " << text << std::endl;
}
void Logger::set_logger_file(const std::string& log_file) {
    s_logger_file = std::ofstream(log_file, std::ostream::app);
}
void Logger::flush() {
    if (s_logger_file && s_logger_file.good())
        s_logger_file << std::flush;
    else
        std::cout << std::flush;
}

// these are the static defined functions
namespace Log {
    void info(const std::string& text) {
        s_logger->info(text);
    }
    void warn(const std::string& text) {
        s_logger->warn(text);
    }
    void critical(const std::string& text) {
        s_logger->critical(text);
    }
    void debug(const std::string& text) {
        s_logger->debug(text);
    }
    void set_logger_file(const std::string& log_file) {
        s_logger->set_logger_file(log_file);
    }
    void set_logger(Logger* logger) {
        s_logger.reset(logger);
    }
    void flush() {
        s_logger->flush();
    }
}