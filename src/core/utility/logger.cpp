#include <iostream>
#include <fstream>
#include "logger.h"

static Logger::Ptr_t s_logger(new Logger());
static std::ofstream s_logger_file;

void Logger::info(const std::string& text) {
    if (s_logger_file && s_logger_file.good())
        s_logger_file << "INFO: " << text << "\n";
    else
        std::cout << "INFO: " << text << "\n";
}
void Logger::warn(const std::string& text) {
    if (s_logger_file && s_logger_file.good())
        s_logger_file << "WARN: " << text << "\n";
    else
        std::cout << "WARN: " << text << "\n";
}
void Logger::critical(const std::string& text) {
    if (s_logger_file && s_logger_file.good())
        s_logger_file << "CRITICAL: " << text << "\n";
    else
        std::cout << "CRITICAL: " << text << "\n";
}
void Logger::debug(const std::string& text) {
    if (s_logger_file && s_logger_file.good())
        s_logger_file << "DEBUG: " <<  text << "\n";
    else
        std::cout << "DEBUG: " << text << "\n";
}
void Logger::set_logger_file(const std::string& log_file) {
    s_logger_file = std::ofstream(log_file);
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
}