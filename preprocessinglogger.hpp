#ifndef __PREPROCESSING_LOGGER_HPP__
#define __PREPROCESSING_LOGGER_HPP__

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>

// Simple logger for the preprocessing server in CSV format
class PreprocessingLogger {
private:
  std::ofstream log_file;
  std::string timestamp;

  // Get current timestamp in ISO 8601 format
  std::string get_timestamp() {
    auto now = std::time(nullptr);
    auto tm = *std::localtime(&now);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
    return oss.str();

  }

  // Create logs directory if it doesn't exist
  void ensure_logs_directory() {
    struct stat st = {0};
    if (stat("logs", &st) == -1) {
      mkdir("logs", 0755);
    }
  }

public:
  PreprocessingLogger() {
    ensure_logs_directory();
    timestamp = get_timestamp();

    // Open log file
    // Filename format: logs/preprocessing_server_<TIMESTAMP>.csv
    // Using a sanitized timestamp for filename to avoid issues with colons
    std::string filename_timestamp = timestamp;
    for (char &c : filename_timestamp) {
      if (c == ':')
        c = '-';
    }

    std::string log_filename =
        "logs/preprocessing_" + filename_timestamp + ".csv";
    log_file.open(log_filename, std::ios::out);

    if (log_file.is_open()) {
      // CSV Header
      log_file << "Timestamp,Message,Duration(s)\n";
    }
  }

  ~PreprocessingLogger() {
    if (log_file.is_open()) {
      log_file.close();
    }
  }

  void log_info(const std::string &message) {
    std::string ts = get_timestamp();
    // Also print to stdout for immediate feedback
    std::cout << "[" << ts << "] " << message << "\n";

    // Write to CSV
    if (log_file.is_open()) {
      log_file << ts << "," << message << ",\n";
      log_file.flush();
    }
  }

  void log_duration(const std::string &operation, double duration_seconds) {
    std::string ts = get_timestamp();
    // Also print to stdout for immediate feedback
    std::cout << "[" << ts << "] " << operation << " took " << duration_seconds
              << "s\n";

    // Write to CSV
    if (log_file.is_open()) {
      log_file << ts << "," << operation << "," << std::fixed
               << std::setprecision(4) << duration_seconds << "\n";
      log_file.flush();
    }
  }
};

#endif
