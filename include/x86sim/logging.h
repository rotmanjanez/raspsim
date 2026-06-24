#ifndef PTLSIM_LOGGING_H
#define PTLSIM_LOGGING_H

#include <print>
#include <format>
#include <memory>
#include <cstdio>
#include <string>
#include <iterator>
#include <utility>
#include <filesystem>

namespace x86sim {
// Global enable flag for backward compatibility
extern bool logenable;

namespace logging {

// Python-compatible log levels
enum class Level : int {
  VERBOSE = 1,  // Maximum detail
  TRACE = 5,    // Below DEBUG
  DEBUG = 10,   // matches logging.DEBUG
  INFO = 20,    // matches logging.INFO (default)
  WARNING = 30, // matches logging.WARNING
  ERROR = 40,   // matches logging.ERROR
  CRITICAL = 50 // matches logging.CRITICAL
};

// For convenience in Level specification
inline constexpr Level VERBOSE = Level::VERBOSE;
inline constexpr Level TRACE = Level::TRACE;
inline constexpr Level DEBUG = Level::DEBUG;
inline constexpr Level INFO = Level::INFO;
inline constexpr Level WARNING = Level::WARNING;
inline constexpr Level ERROR = Level::ERROR;
inline constexpr Level CRITICAL = Level::CRITICAL;

// Abstract sink interface - pluggable backends
class LogSink {
public:
  virtual ~LogSink() = default;

  // Write a log message (level already checked by Logger)
  virtual void write(Level level, const std::string& message) = 0;

  // Flush buffered output
  virtual void flush() = 0;
};

// FILE* based sink (default implementation)
class FileLogSink : public LogSink {
private:
  FILE* fp_;
  bool owns_file_;

public:
  explicit FileLogSink(FILE* fp, bool owns_file = false) : fp_(fp), owns_file_(owns_file) {}

  explicit FileLogSink(const std::filesystem::path& filename) : fp_(nullptr), owns_file_(true) {
    std::string p = filename.string();
    fp_ = std::fopen(p.c_str(), "w");
    if (!fp_) {
      fp_ = stderr;
      owns_file_ = false;
    }
  }

  ~FileLogSink() override {
    if (owns_file_ && fp_) {
      std::fclose(fp_);
    }
  }

  void write(Level level, const std::string& message) override {
    if (fp_) {
      std::print(fp_, "{}", message);
    }
  }

  void flush() override {
    if (fp_) {
      std::fflush(fp_);
    }
  }
};

// Forward declaration
class Logger;

// Output iterator for std::format_to support
// Uses shared buffer to handle copying correctly
class log_output_iterator {
public:
  struct buffer_state {
    std::string buffer;
    Logger* logger;
    Level level;

    ~buffer_state(); // Defined after Logger is complete
  };

private:
  std::shared_ptr<buffer_state> state_;

public:
  using iterator_category = std::output_iterator_tag;
  using value_type = void;
  using difference_type = std::ptrdiff_t;
  using pointer = void;
  using reference = void;

  log_output_iterator(Logger* logger, Level level) : state_(std::make_shared<buffer_state>()) {
    state_->logger = logger;
    state_->level = level;
  }

  // Copyable (shares buffer state)
  log_output_iterator(const log_output_iterator&) = default;
  log_output_iterator& operator=(const log_output_iterator&) = default;

  // Output iterator interface
  log_output_iterator& operator=(char c) {
    if (state_) {
      state_->buffer += c;
    }
    return *this;
  }

  log_output_iterator& operator*() { return *this; }
  log_output_iterator& operator++() { return *this; }
  log_output_iterator operator++(int) { return *this; }

  // Manual flush (defined after Logger is complete)
  void flush();
};

// Singleton Logger - manages sink and level checking
class Logger {
private:
  friend class log_output_iterator;

  std::unique_ptr<LogSink> sink_;
  int current_level_;

  Logger() : sink_(std::make_unique<FileLogSink>(stderr)), current_level_(static_cast<int>(INFO)) {}

  void write_to_sink(Level level, const std::string& message) {
    if (sink_) {
      sink_->write(level, message);
    }
  }

  void flush_sink() {
    if (sink_) {
      sink_->flush();
    }
  }

public:
  static Logger& instance() {
    static Logger logger;
    return logger;
  }

  Logger(const Logger&) = delete;
  Logger& operator=(const Logger&) = delete;

  int current_level() const { return current_level_; }

  void set_level(Level level) { current_level_ = static_cast<int>(level); }

  void set_level(int level) { current_level_ = level; }

  void set_sink(std::unique_ptr<LogSink> sink) {
    if (sink) {
      sink_ = std::move(sink);
    }
  }

  void set_file_sink(const std::filesystem::path& filename) { set_sink(std::make_unique<FileLogSink>(filename)); }
  void set_file_sink(const char* filename) { set_file_sink(std::filesystem::path(filename)); }
  void set_file_sink(FILE* fp, bool owns_file = false) { set_sink(std::make_unique<FileLogSink>(fp, owns_file)); }

  // Return output iterator for std::format_to
  log_output_iterator output(Level level) { return log_output_iterator{this, level}; }

  // Explicit flush
  void flush() { flush_sink(); }
};

// Define buffer_state destructor after Logger is complete
inline log_output_iterator::buffer_state::~buffer_state() {
  if (!buffer.empty() && logger) {
    logger->write_to_sink(level, buffer);
  }
}

// Define flush() after Logger is complete
inline void log_output_iterator::flush() {
  if (state_ && !state_->buffer.empty() && state_->logger) {
    state_->logger->write_to_sink(state_->level, state_->buffer);
    state_->buffer.clear();
  }
}

// Fast inline level check - NO virtual call on hot path
inline bool logable_at(Level level) {
  return logenable && (static_cast<int>(level) >= Logger::instance().current_level());
}

// Convenience function: get output iterator for formatting
inline log_output_iterator output(Level level) {
  return Logger::instance().output(level);
}

// Templated print functions for easy migration from old API
template<typename... Args>
inline void print(Level lvl, std::format_string<Args...> fmt, Args&&... args) {
  if (logable_at(lvl)) {
    auto out = output(lvl);
    std::format_to(out, fmt, std::forward<Args>(args)...);
  }
}

template<typename... Args>
inline void println(Level lvl, std::format_string<Args...> fmt, Args&&... args) {
  if (logable_at(lvl)) {
    auto out = output(lvl);
    std::format_to(out, fmt, std::forward<Args>(args)...);
    std::format_to(out, "\n");
  }
}

// Default to INFO level
template<typename... Args>
inline void print(std::format_string<Args...> fmt, Args&&... args) {
  print(INFO, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void println(std::format_string<Args...> fmt, Args&&... args) {
  println(INFO, fmt, std::forward<Args>(args)...);
}

// Flush output
inline void flush() {
  Logger::instance().flush();
}

// Convenience functions for API compatibility
inline void set_level(Level level) {
  Logger::instance().set_level(level);
}

inline void set_level(int level) {
  Logger::instance().set_level(level);
}

inline void set_file_sink(const std::filesystem::path& filename) { Logger::instance().set_file_sink(filename); }
inline void set_file_sink(const char* filename) { Logger::instance().set_file_sink(filename); }
inline void set_file_sink(FILE* fp, bool owns_file = false) { Logger::instance().set_file_sink(fp, owns_file); }

template<typename... Args>
inline void eprint(std::format_string<Args...> fmt, Args&&... args) {
  std::print(stderr, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void eprintln(std::format_string<Args...> fmt, Args&&... args) {
  std::println(stderr, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void oprint(std::format_string<Args...> fmt, Args&&... args) {
  std::print(stdout, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void oprintln(std::format_string<Args...> fmt, Args&&... args) {
  std::println(stdout, fmt, std::forward<Args>(args)...);
}

} // namespace logging
} // namespace x86sim

#endif // PTLSIM_LOGGING_H
