#ifndef __LOGGER_HPP__
#define __LOGGER_HPP__

#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <iomanip>
#include <sys/stat.h>
#include "types.hpp"

// Performance logging utility for segment tree experiments
class PerformanceLogger {
private:
    std::ofstream csv_file;
    std::ofstream log_file;
    std::string experiment_id;
    size_t depth;
    size_t num_updates;
    size_t num_queries;
    std::string mode_str;
    
    // Get current timestamp in ISO 8601 format
    std::string get_timestamp() {
        auto now = std::time(nullptr);
        auto tm = *std::localtime(&now);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S");
        return oss.str();
    }
    
    // Create logs directory if it doesn't exist
    void ensure_logs_directory() {
        struct stat st = {0};
        if (stat("logs", &st) == -1) {
            mkdir("logs", 0755);
        }
    }
    
    // Convert ProcessingMode to string
    std::string mode_to_string(ProcessingMode mode) {
        switch (mode) {
            case MODE_ONLINE: return "online";
            case MODE_PREPROCESSING: return "preprocessing";
            case MODE_ONLINEONLY: return "onlineonly";
            default: return "unknown";
        }
    }
    
public:
    PerformanceLogger(const std::string& exp_id, size_t d, size_t n_u, size_t n_q, int player, ProcessingMode mode) 
        : experiment_id(exp_id), depth(d), num_updates(n_u), num_queries(n_q), mode_str(mode_to_string(mode)) {
        
        // Only create log files for player 0
        if (player != 0) {
            return;
        }
        
        ensure_logs_directory();
        
        // Open CSV file for performance metrics (include mode in filename)
        std::string csv_filename = "logs/performance_" + mode_str + "_" + experiment_id + ".csv";
        csv_file.open(csv_filename, std::ios::out);
        
        // Write CSV header (include mode column)
        csv_file << "timestamp,experiment_id,mode,depth,num_updates,num_queries,"
                 << "phase,operation_id,metric,value,unit\n";
        
        // Open log file for full execution output (include mode in filename)
        std::string log_filename = "logs/execution_" + mode_str + "_" + experiment_id + ".log";
        log_file.open(log_filename, std::ios::out);
        
        // Write log header
        log_file << "========================================\n";
        log_file << "Experiment ID: " << experiment_id << "\n";
        log_file << "Mode: " << mode_str << "\n";
        log_file << "Timestamp: " << get_timestamp() << "\n";
        log_file << "Depth: " << depth << ", Updates: " << num_updates 
                 << ", Queries: " << num_queries << "\n";
        log_file << "========================================\n\n";
    }
    
    ~PerformanceLogger() {
        if (csv_file.is_open()) csv_file.close();
        if (log_file.is_open()) log_file.close();
    }
    
    // Log a performance metric to CSV
    void log_metric(const std::string& phase, size_t operation_id, 
                   const std::string& metric, double value, 
                   const std::string& unit = "ms") {
        if (!csv_file.is_open()) return; // Skip if not player 0
        csv_file << get_timestamp() << ","
                 << experiment_id << ","
                 << mode_str << ","
                 << depth << ","
                 << num_updates << ","
                 << num_queries << ","
                 << phase << ","
                 << operation_id << ","
                 << metric << ","
                 << value << ","
                 << unit << "\n";
        csv_file.flush(); // Ensure immediate write
    }
    
    // Log general execution output
    void log_output(const std::string& message) {
        if (log_file.is_open()) {
            log_file << message;
            log_file.flush();
        }
        std::cout << message; // Always print to console
    }
    
    // Log a section header
    void log_section(const std::string& title) {
        std::string header = "\n===== " + title + " =====\n";
        if (log_file.is_open()) {
            log_file << header;
            log_file.flush();
        }
        std::cout << header;
    }
    
    // Log network/computation stats
    void log_stats(const std::string& phase, size_t messages_sent, 
                  size_t message_bytes, size_t lamport_clock,
                  size_t aes_ops, size_t wall_clock_ms, size_t mem_kb) {
        if (!csv_file.is_open()) return; // Skip if not player 0
        log_metric(phase, 0, "messages_sent", messages_sent, "count");
        log_metric(phase, 0, "message_bytes", message_bytes, "bytes");
        log_metric(phase, 0, "lamport_clock", lamport_clock, "latencies");
        log_metric(phase, 0, "aes_operations", aes_ops, "count");
        log_metric(phase, 0, "wall_clock_time", wall_clock_ms, "ms");
        log_metric(phase, 0, "memory_usage", mem_kb, "KiB");
    }
};

#endif
