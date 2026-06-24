// -*- c++ -*-
//
// Branch Prediction
//
// Copyright 2003-2008 Matt T. Yourst <yourst@yourst.com>
//
// This program is free software; it is licensed under the
// GNU General Public License, Version 2.
// TODO: change to alexis

#ifndef _ADDRSPACE_H_
#define _ADDRSPACE_H_

#include <unordered_map>
#include <memory>

#include "x86sim/x86sim.hpp"
#include "globals.h"
#include "superstl.h"

#include "ptlsim.h"

//
// Address space management
//
namespace x86sim {

// Each chunk covers 2 GB of virtual address space:
#define SPAT_TOPLEVEL_CHUNK_BITS 17
#define SPAT_PAGES_PER_CHUNK_BITS 19
#define SPAT_TOPLEVEL_CHUNKS (1 << SPAT_TOPLEVEL_CHUNK_BITS)  // 262144
#define SPAT_PAGES_PER_CHUNK (1 << SPAT_PAGES_PER_CHUNK_BITS) // 524288
#define SPAT_BYTES_PER_CHUNK (SPAT_PAGES_PER_CHUNK / 8)       // 65536
#define ADDRESS_SPACE_BITS (48)
#define ADDRESS_SPACE_SIZE (1LL << ADDRESS_SPACE_BITS)


class AddressSpace {
public:
  AddressSpace();
  ~AddressSpace();

public:
  using Page = std::array<W8, PAGE_SIZE>;
  std::unordered_map<Waddr, std::unique_ptr<Page>> mapped_mem;

  void map(Waddr start, Waddr length, Protection prot) {
    start = floor(start, PAGE_SIZE);
    length = ceil(length, PAGE_SIZE);
    Waddr num_pages = length / PAGE_SIZE;
    foreach (i, num_pages) {
      mapped_mem.insert_or_assign(start + i * PAGE_SIZE, std::make_unique<Page>());
    }
    setattr((byte*)start, length, prot);
  }
  void unmap(Waddr start, Waddr length) {
    start = floor(start, PAGE_SIZE);
    length = ceil(length, PAGE_SIZE);
    Waddr num_pages = length / PAGE_SIZE;
    foreach (i, num_pages) {
      mapped_mem.erase(start + i * PAGE_SIZE);
    }
    setattr((byte*)start, length, Protection::none);
  }

  void* page_virt_to_mapped(Waddr addr) {
    auto it = mapped_mem.find(floor(addr, PAGE_SIZE));
    if (it == mapped_mem.end())
      return nullptr;

    W8* base = it->second.get()->data();
    return base + lowbits(addr, 12);
  }

  //
  // Shadow page attribute table
  //
  typedef byte SPATChunk[SPAT_BYTES_PER_CHUNK];
  typedef SPATChunk** spat_t;

  spat_t readmap = nullptr;
  spat_t writemap = nullptr;
  spat_t execmap = nullptr;
  spat_t dirtymap = nullptr;

  spat_t allocmap();
  void freemap(spat_t top);

  byte& pageid_to_map_byte(spat_t top, Waddr pageid);
  void make_accessible(void* address, Waddr size, spat_t top);
  void make_inaccessible(void* address, Waddr size, spat_t top);

  Waddr pageid(void* address) const { return ((W64)lowbits((W64)address, ADDRESS_SPACE_BITS)) >> log2(PAGE_SIZE); }

  Waddr pageid(Waddr address) const { return pageid((void*)address); }

  void make_page_accessible(void* address, spat_t top) {
    setbit(pageid_to_map_byte(top, pageid(address)), lowbits(pageid(address), 3));
  }

  void make_page_inaccessible(void* address, spat_t top) {
    clearbit(pageid_to_map_byte(top, pageid(address)), lowbits(pageid(address), 3));
  }

  void allow_read(void* address, Waddr size) { make_accessible(address, size, readmap); }
  void disallow_read(void* address, Waddr size) { make_inaccessible(address, size, readmap); }
  void allow_write(void* address, Waddr size) { make_accessible(address, size, writemap); }
  void disallow_write(void* address, Waddr size) { make_inaccessible(address, size, writemap); }
  void allow_exec(void* address, Waddr size) { make_accessible(address, size, execmap); }
  void disallow_exec(void* address, Waddr size) { make_inaccessible(address, size, execmap); }

public:
  //
  // Memory management passthroughs
  //
  void setattr(void* start, Waddr length, Protection prot);
  Protection getattr(void* start);

  bool fastcheck(Waddr addr, spat_t top) const {

    // Is it outside of userspace address range?
    // Check disabled to allow access to VDSO in kernel space.
    if unlikely (addr >> 48)
      return 0;

    W64 chunkid = pageid(addr) >> log2(SPAT_PAGES_PER_CHUNK);

    if unlikely (!top[chunkid])
      return false;

    AddressSpace::SPATChunk& chunk = *top[chunkid];
    Waddr byteid = bits(pageid(addr), 3, log2(SPAT_BYTES_PER_CHUNK));
    return bit(chunk[byteid], lowbits(pageid(addr), 3));
  }

  bool fastcheck(void* addr, spat_t top) const { return fastcheck((Waddr)addr, top); }

  bool check(void* p, Protection prot) const {
    if unlikely (has_protection(prot, Protection::read) && (!fastcheck(p, readmap)))
      return false;

    if unlikely (has_protection(prot, Protection::write) && (!fastcheck(p, writemap)))
      return false;

    if unlikely (has_protection(prot, Protection::execute) && (!fastcheck(p, execmap)))
      return false;

    return true;
  }

  bool isdirty(Waddr mfn) { return fastcheck(mfn << 12, dirtymap); }
  void setdirty(Waddr mfn) { make_page_accessible((void*)(mfn << 12), dirtymap); }
  void cleardirty(Waddr mfn) { make_page_inaccessible((void*)(mfn << 12), dirtymap); }

  void resync_with_process_maps();
};

} // namespace x86sim

#endif
