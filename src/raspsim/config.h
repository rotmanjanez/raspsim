// -*- c++ -*-
//
// Copyright 2005-2008 Matt T. Yourst <yourst@yourst.com>
//
// This program is free software; it is licensed under the
// GNU General Public License, Version 2.
//
#ifndef _CONFIG_H_
#define _CONFIG_H_

#include "globals.h"
#include "superstl.h"
#include "x86sim/x86sim.hpp"
#include <stdarg.h>

static const x86sim::W64 infinity = std::numeric_limits<x86sim::W64s>::max();

struct ConfigurationOption {
  const char* name;
  const char* description;
  int type;
  int fieldsize;
  x86sim::Waddr offset;

  ConfigurationOption* next;

  ConfigurationOption() {}

  ConfigurationOption(const char* name, const char* description, int type, x86sim::Waddr offset, int fieldsize = 0) {
    this->name = name;
    this->description = description;
    this->type = type;
    this->offset = offset;
    this->fieldsize = fieldsize;
    this->next = nullptr;
  }
};

enum {
  OPTION_TYPE_NONE = 0,
  OPTION_TYPE_W64 = 1,
  OPTION_TYPE_FLOAT = 2,
  OPTION_TYPE_STRING = 3,
  OPTION_TYPE_TRAILER = 4,
  OPTION_TYPE_BOOL = 5,
  OPTION_TYPE_CORE_MODEL = 6,
  OPTION_TYPE_PATH = 7,
  OPTION_TYPE_SECTION = -1
};

struct ConfigurationParserBase {
  ConfigurationOption* options;
  ConfigurationOption* lastoption;

  void addentry(void* baseptr, void* field, int type, const char* name, const char* description) {
    x86sim::Waddr offset = ((x86sim::Waddr)field) - ((x86sim::Waddr)baseptr);
    ConfigurationOption* option = new ConfigurationOption(name, description, type, offset);
    if (lastoption)
      lastoption->next = option;
    if (!options)
      options = option;
    lastoption = option;
  }

  ConfigurationParserBase() {
    options = nullptr;
    lastoption = nullptr;
  }

  int parse(void* baseptr, int argc, char* argv[]);
  std::string format_to_string_usage(const void* baseptr) const;
  std::string format_to_string_config(const void* baseptr) const;
};

template<typename T>
struct ConfigurationParser : public T {
  ConfigurationParserBase options;

  ConfigurationParser() {}

  void add(x86sim::W64& field, const char* name, const char* description) {
    options.addentry(this, &field, OPTION_TYPE_W64, name, description);
  }

  void add(double& field, const char* name, const char* description) {
    options.addentry(this, &field, OPTION_TYPE_FLOAT, name, description);
  }

  void add(bool& field, const char* name, const char* description) {
    options.addentry(this, &field, OPTION_TYPE_BOOL, name, description);
  }

  void add(x86sim::CoreModel& field, const char* name, const char* description) {
    options.addentry(this, &field, OPTION_TYPE_CORE_MODEL, name, description);
  }

  void add(std::string& field, const char* name, const char* description) {
    options.addentry(this, &field, OPTION_TYPE_STRING, name, description);
  }

  void add(std::filesystem::path& field, const char* name, const char* description) {
    options.addentry(this, &field, OPTION_TYPE_PATH, name, description);
  }

  void section(const char* name) { options.addentry(this, nullptr, OPTION_TYPE_SECTION, name, name); }

  int parse(T& config, int argc, char* argv[]) { return options.parse(&config, argc, argv); }

  // Provided by user:
  void setup();
};

#endif // _CONFIG_H_
