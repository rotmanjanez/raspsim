// -*- c++ -*-
//
// Super Standard Template Library
//
// Faster and more optimized than stock STL implementation,
// plus includes various customized features
//
// Copyright 1997-2008 Matt T. Yourst <yourst@yourst.com>
//
// This program is free software; it is licensed under the
// GNU General Public License, Version 2.
//

#ifndef _SUPERSTL_H_
#define _SUPERSTL_H_

#include <algorithm>
#include <bitset>
#include <bit>
#include <format>
#include <print>
#include <string>
#include <string_view>
#include <vector>

namespace x86sim {

//
// Formatting
//
#define FMT_ZEROPAD 1  /* pad with zero */
#define FMT_SIGN 2     /* unsigned/signed long */
#define FMT_PLUS 4     /* show plus */
#define FMT_SPACE 8    /* space if plus */
#define FMT_LEFT 16    /* left justified */
#define FMT_SPECIAL 32 /* 0x */
#define FMT_LARGE 64   /* use 'ABCDEF' instead of 'abcdef' */

char* format_number(char* buf, char* end, W64 num, int base, int size, int precision, int type);
int format_integer(char* buf, int bufsize, W64s v, int size = 0, int flags = 0, int base = 10, int precision = 0);
int format_float(char* buf, int bufsize, double v, int precision = 6, int pad = 0);

int current_vcpuid();

extern bool force_synchronous_streams;

namespace superstl {
// Print bits as a string:
struct bitstring {
  W64 bits;
  int n;
  bool reverse;

  bitstring() {}

  bitstring(const W64 bits, const int n, bool reverse = false) {
    assert(n <= 64);
    this->bits = bits;
    this->n = n;
    this->reverse = reverse;
  }
};

struct hexstring {
  W64 value;
  int n;

  hexstring() {}

  hexstring(const W64 value, const int n) {
    this->value = value;
    this->n = n;
  }
};

struct bytemaskstring {
  const byte* bytes;
  W64 mask;
  int n;
  int splitat;

  bytemaskstring() {}

  bytemaskstring(const byte* bytes, W64 mask, int n, int splitat = 16) {
    assert(n <= 64);
    this->bytes = bytes;
    this->mask = mask;
    this->n = n;
    this->splitat = splitat;
  }
};

struct intstring {
  W64s value;
  int width;

  intstring() {}

  intstring(W64s value, int width) {
    this->value = value;
    this->width = width;
  }
};

struct floatstring {
  double value;
  int width;
  int precision;

  floatstring() {}

  floatstring(double value, int width = 0, int precision = 6) {
    this->value = value;
    this->width = width;
    this->precision = precision;
  }
};

struct padstring {
  const char* value;
  int width;
  char pad;

  padstring() {}

  padstring(const char* value, int width, char pad = ' ') {
    this->value = value;
    this->width = width;
    this->pad = pad;
  }
};

struct percentstring {
  double fraction;
  int width;

  percentstring() {}

  percentstring(W64s value, W64s total, int width = 7) {
    fraction = (total) ? (double(value) / double(total)) : 0;
    this->width = width;
  }
};

struct substring {
  const char* str;
  int length;

  substring() {}

  substring(const char* str, int start, int length) {
    int r = strlen(str);
    this->length = std::min(length, r - start);
    this->str = str + std::min(start, r);
  }
};

/*
   * CRC32
   */
struct CRC32 {
  static const W32 crctable[256];
  W32 crc;

  inline W32 update(byte value) {
    crc = crctable[byte(crc ^ value)] ^ (crc >> 8);
    return crc;
  }

  inline W32 update(byte* data, unsigned count) {
    foreach (i, count) {
      update(data[i]);
    }
    return crc;
  }

  CRC32() { reset(); }

  CRC32(W32 newcrc) { reset(newcrc); }

  inline void reset(W32 newcrc = 0xffffffff) { crc = newcrc; }

  operator W32() const { return crc; }
};

template<typename T>
static inline CRC32& operator<<(CRC32& crc, const T& t) {
  crc.update((byte*)&t, sizeof(T));
  return crc;
}

template<class T>
static inline CRC32& operator,(CRC32& crc, const T& v) {
  return crc << v;
}

//
// selflistlink class
// Double linked list without pointer: useful as root
// of inheritance hierarchy for another class to save
// space, since object pointed to is implied
//
class selflistlink {
public:
  selflistlink* next;
  selflistlink* prev;

public:
  void reset() {
    next = null;
    prev = null;
  }
  selflistlink() { reset(); }

  selflistlink* unlink() {
    if likely (prev)
      prev->next = next;
    if likely (next)
      next->prev = prev;
    prev = null;
    next = null;
    return this;
  }

  selflistlink* replacewith(selflistlink* newlink) {
    if likely (prev)
      prev->next = newlink;
    if likely (next)
      next->prev = newlink;
    newlink->prev = prev;
    newlink->next = next;
    return newlink;
  }

  void addto(selflistlink*& root) {
    // THIS <-> root <-> a <-> b <-> c
    this->prev = (selflistlink*)&root;
    this->next = root;
    if likely (root)
      root->prev = this;
    // Do not touch root->next since it might not even exist
    root = this;
  }

  bool linked() const { return (next || prev); }

  bool unlinked() const { return !linked(); }
};

class selfqueuelink {
public:
  selfqueuelink* next;
  selfqueuelink* prev;

public:
  void reset() {
    next = this;
    prev = this;
  }
  selfqueuelink() {}

  selfqueuelink& unlink() {
    // No effect if next = prev = this (i.e., unlinked)
    next->prev = prev;
    prev->next = next;
    prev = this;
    next = this;
    return *this;
  }

  void addhead(selfqueuelink& root) { addlink(&root, root.next); }

  void addhead(selfqueuelink* root) { addhead(*root); }

  void addto(selfqueuelink& root) { addhead(root); }

  void addto(selfqueuelink* root) { addto(*root); }

  void addtail(selfqueuelink& root) { addlink(root.prev, &root); }

  void addtail(selfqueuelink* root) { addtail(*root); }

  selfqueuelink* removehead() {
    if unlikely (empty())
      return null;
    selfqueuelink* link = next;
    link->unlink();
    return link;
  }

  selfqueuelink* removetail() {
    if unlikely (empty())
      return null;
    selfqueuelink* link = prev;
    link->unlink();
    return link;
  }

  selfqueuelink* head() const { return next; }

  selfqueuelink* tail() const { return prev; }

  bool empty() const { return (next == this); }

  bool unlinked() const { return ((!prev && !next) || ((prev == this) && (next == this))); }

  bool linked() const { return !unlinked(); }

  operator bool() const { return (!empty()); }

protected:
  void addlink(selfqueuelink* prev, selfqueuelink* next) {
    next->prev = this;
    this->next = next;
    this->prev = prev;
    prev->next = this;
  }
};

//
// Default link manager for objects in which the
// very first member (or superclass) is selflistlink.
//
template<typename T>
struct ObjectLinkManager {
  static inline T* objof(selflistlink* link) { return (T*)link; }
  static inline selflistlink* linkof(T* obj) { return (selflistlink*)obj; }
  //
  // Example:
  //
  // T* objof(selflistlink* link) {
  //   return baseof(T, hashlink, link); // a.k.a. (T*)((byte*)link) - offsetof(T, hashlink);
  // }
  //
  // selflistlink* linkof(T* obj) {
  //   return &obj->link;
  // }
  //
};

template<class T>
class queuelink {
public:
  queuelink<T>* next;
  queuelink<T>* prev;
  T* data;

public:
  void reset() {
    next = this;
    prev = this;
    data = null;
  }
  queuelink() { reset(); }
  queuelink(const T& t) {
    reset();
    data = &t;
  }
  queuelink(const T* t) {
    reset();
    data = t;
  }
  queuelink<T>& operator()(T* t) {
    reset();
    data = t;
    return *this;
  }

  T& unlink() {
    // No effect if next = prev = this (i.e., unlinked)
    next->prev = prev;
    prev->next = next;
    prev = this;
    next = this;
    return *data;
  }

  void add_to_head(queuelink<T>& root) { addlink(&root, root.next); }

  void addto(queuelink<T>& root) { addhead(root); }

  void add_to_tail(queuelink<T>& root) { addlink(root.prev, &root); }

  queuelink<T>* remove_head() {
    queuelink<T>* link = next;
    link->unlink();
    return link;
  }

  queuelink<T>* remove_tail() {
    queuelink<T>* link = prev;
    link->unlink();
    return link;
  }

  queuelink<T>* head() const { return next; }

  queuelink<T>* tail() const { return prev; }

  bool empty() const { return (next == this); }

  bool unlinked() const { return ((!prev && !next) || ((prev == this) && (next == this))); }

  bool linked() const { return !unlinked(); }

  operator bool() const { return (!empty()); }

  T* operator->() const { return data; }
  operator T*() const { return data; }
  operator T&() const { return *data; }

protected:
  void addlink(queuelink<T>* prev, queuelink<T>* next) {
    next->prev = this;
    this->next = next;
    this->prev = prev;
    prev->next = this;
  }
};

template<typename T, typename LM = superstl::ObjectLinkManager<T>>
class queue : selflistlink {
public:
  void reset() {
    next = this;
    prev = this;
  }
  queue() { reset(); }

  void add_to_head(selflistlink* link) { addlink(this, link, next); }
  void add_to_head(T& obj) { add_to_head(LM::linkof(&obj)); }
  void add_to_head(T* obj) { add_to_head(LM::linkof(obj)); }

  void add_to_tail(selflistlink* link) { addlink(prev, link, this); }
  void add_to_tail(T& obj) { add_to_tail(LM::linkof(&obj)); }
  void add_to_tail(T* obj) { add_to_tail(LM::linkof(obj)); }

  T* remove_head() {
    if unlikely (empty())
      return null;
    selflistlink* link = next;
    link->unlink();
    return LM::objof(link);
  }

  T* remove_tail() {
    if unlikely (empty())
      return null;
    selflistlink* link = prev;
    link->unlink();
    return LM::objof(link);
  }

  void enqueue(T* obj) { add_to_tail(obj); }
  T* dequeue() { return remove_head(); }

  void push(T* obj) { add_to_tail(obj); }
  void pop(T* obj) { remove_tail(); }

  T* head() const { return (unlikely(empty())) ? null : next; }

  T* tail() const { return (unlikely(empty())) ? null : prev; }

  bool empty() const { return (next == this); }

  operator bool() const { return (!empty()); }

protected:
  void addlink(selflistlink* prev, selflistlink* link, selflistlink* next) {
    next->prev = link;
    link->next = next;
    link->prev = prev;
    prev->next = link;
  }
};


//
// Index References (indexrefs) work exactly like pointers but always
// index into a specific structure. This saves considerable space and
// can allow aliasing optimizations not possible with pointers.
//

template<typename T, typename P = W32, Waddr base = 0, int granularity = 1>
struct shortptr {
  P p;

  shortptr() {}

  shortptr(const T& obj) { *this = obj; }

  shortptr(const T* obj) { *this = obj; }

  shortptr<T, P, base, granularity>& operator=(const T& obj) {
    p = (P)((((Waddr)&obj) - base) / granularity);
    return *this;
  }

  shortptr<T, P, base, granularity>& operator=(const T* obj) {
    p = (P)((((Waddr)obj) - base) / granularity);
    return *this;
  }

  T* get() const { return (T*)((p * granularity) + base); }

  T* operator->() const { return get(); }

  T& operator*() const { return *get(); }

  operator T*() const { return get(); }

  shortptr<T, P, base, granularity>& operator++() {
    (*this) = (get() + 1);
    return *this;
  }

  shortptr<T, P, base, granularity>& operator--() {
    (*this) = (get() - 1);
    return *this;
  }
};


// Index of the least significant set bit (N if no bit is set)
template<std::size_t N>
inline int lsb(const std::bitset<N>& bits) {
  return std::countr_zero(bits.to_ullong());
}

// Remove the bit at the given index, shifting all higher bits down by one
template<std::size_t N>
inline std::bitset<N> remove_bit(const std::bitset<N>& bits, int index) {
  const unsigned long long v = bits.to_ullong();
  const unsigned long long lo = v & ((1ULL << index) - 1);
  const unsigned long long hi = (v >> index) >> 1;
  return std::bitset<N>(lo | (hi << index));
}

//
// Convenient list iterator
//
#define foreachlink(list, type, iter)                                                                                  \
  for (type* iter = (type*)((list)->first); (iter != NULL); __builtin_prefetch(iter->next), iter = (type*)(iter->next))

template<typename K, typename T>
struct KeyValuePair {
  T value;
  K key;
};

template<typename K, int setcount>
struct HashtableKeyManager {
  static inline int hash(const K& key);
  static inline bool equal(const K& a, const K& b);
  static inline K dup(const K& key);
  static inline void free(K& key);
};

template<int setcount>
struct HashtableKeyManager<W64, setcount> {
  static inline int hash(W64 key) { return foldbits<log2(setcount)>(key); }

  static inline bool equal(W64 a, W64 b) { return (a == b); }
  static inline W64 dup(W64 key) { return key; }
  static inline void free(W64 key) {}
};

template<int setcount>
struct HashtableKeyManager<const char*, setcount> {
  static inline int hash(const char* key) {
    int len = strlen(key);
    CRC32 h;
    foreach (i, len) {
      h << key[i];
    }
    return h;
  }

  static inline bool equal(const char* a, const char* b) { return (strcmp(a, b) == 0); }

  static inline const char* dup(const char* key) { return strdup(key); }

  static inline void free(const char* key) { ::free((void*)key); }
};

template<typename T, typename K>
struct HashtableLinkManager {
  static inline T* objof(selflistlink* link);
  static inline K& keyof(T* obj);
  static inline selflistlink* linkof(T* obj);
  //
  // Example:
  //
  // T* objof(selflistlink* link) {
  //   return baseof(T, hashlink, link); // a.k.a. *(T*)((byte*)link) - offsetof(T, hashlink);
  // }
  //
};

template<typename K, typename T, int setcount = 64, typename LM = ObjectLinkManager<T>,
         typename KM = HashtableKeyManager<K, setcount>>
struct SelfHashtable {
protected:
  selflistlink* sets[setcount];

public:
  int count;

  T* get(const K& key) {
    selflistlink* tlink = sets[lowbits(KM::hash(key), log2(setcount))];
    while (tlink) {
      T* obj = LM::objof(tlink);
      if likely (KM::equal(LM::keyof(obj), key))
        return obj;
      tlink = tlink->next;
    }

    return null;
  }

  struct Iterator {
    SelfHashtable<K, T, setcount, LM, KM>* ht;
    selflistlink* link;
    int slot;

    Iterator() {}

    Iterator(SelfHashtable<K, T, setcount, LM, KM>* ht) { reset(ht); }

    Iterator(SelfHashtable<K, T, setcount, LM, KM>& ht) { reset(ht); }

    void reset(SelfHashtable<K, T, setcount, LM, KM>* ht) {
      this->ht = ht;
      slot = 0;
      link = ht->sets[slot];
    }

    void reset(SelfHashtable<K, T, setcount, LM, KM>& ht) { reset(&ht); }

    T* next() {
      for (;;) {
        if unlikely (slot >= setcount)
          return null;

        if unlikely (!link) {
          // End of chain: advance to next chain
          slot++;
          if unlikely (slot >= setcount)
            return null;
          link = ht->sets[slot];
          continue;
        }

        T* obj = LM::objof(link);
        link = link->next;
        __builtin_prefetch(link);
        return obj;
      }
    }
  };

  std::vector<T*> getentries() {
    std::vector<T*> result(count);
    int n = 0;
    Iterator iter(this);
    T* t;
    while ((t = iter.next())) {
      assert(n < count);
      result[n++] = t;
    }
    return result;
  }

  SelfHashtable() { reset(); }

  void reset() {
    count = 0;
    foreach (i, setcount) {
      sets[i] = null;
    }
  }

  void clear(bool free_after_remove = false) {
    foreach (i, setcount) {
      selflistlink* tlink = sets[i];
      while (tlink) {
        selflistlink* tnext = tlink->next;
        tlink->unlink();
        if unlikely (free_after_remove) {
          T* obj = LM::objof(tlink);
          delete obj;
        }
        tlink = tnext;
      }
      sets[i] = null;
    }
    count = 0;
  }

  void clear_and_free() { clear(true); }

  T* operator()(const K& key) { return get(key); }

  T* add(T* obj) {
    T* oldobj = get(LM::keyof(obj));
    if unlikely (oldobj) {
      remove(oldobj);
    }

    if (LM::linkof(obj)->linked())
      return obj;

    LM::linkof(obj)->addto(sets[lowbits(KM::hash(LM::keyof(obj)), log2(setcount))]);
    count++;
    return obj;
  }

  T& add(T& obj) { return *add(&obj); }

  T* remove(T* obj) {
    selflistlink* link = LM::linkof(obj);
    if (!link->linked())
      return obj;
    link->unlink();
    count--;
    return obj;
  }

  T& remove(T& obj) { return *remove(&obj); }
};


template<typename K, typename T, typename KM>
struct ObjectHashtableEntry : public KeyValuePair<K, T> {
  typedef KeyValuePair<K, T> base_t;
  selflistlink hashlink;

  ObjectHashtableEntry() {}

  ObjectHashtableEntry(const K& key, const T& value) {
    this->value = value;
    this->key = KM::dup(key);
  }

  ~ObjectHashtableEntry() {
    hashlink.unlink();
    KM::free(this->key);
  }
};

template<typename K, typename T, typename KM>
struct ObjectHashtableLinkManager {
  typedef ObjectHashtableEntry<K, T, KM> entry_t;

  static inline entry_t* objof(selflistlink* link) { return baseof(entry_t, hashlink, link); }

  static inline K& keyof(entry_t* obj) { return obj->key; }

  static inline selflistlink* linkof(entry_t* obj) { return &obj->hashlink; }
};

template<typename T, int N>
struct ChunkList {
  struct Chunk;

  struct Chunk {
    selflistlink link;
    std::bitset<N> freemap;

    // Formula: (CHUNK_SIZE - sizeof(ChunkHeader<T>)) / sizeof(T);
    T data[N];

    Chunk() {
      link.reset();
      freemap.set();
    }

    bool full() const { return (freemap.none()); }
    bool empty() const { return freemap.all(); }

    int add(const T& entry) {
      if unlikely (full())
        return -1;
      int idx = lsb(freemap);
      freemap[idx] = 0;
      data[idx] = entry;
      return idx;
    }

    bool contains(T* entry) const {
      int idx = entry - data;
      return ((idx >= 0) & (idx < lengthof(data)));
    }

    bool remove(int idx) {
      data[idx] = 0;
      freemap[idx] = 1;

      return empty();
    }

    struct Iterator {
      Chunk* chunk;
      size_t i;

      Iterator() {}

      Iterator(Chunk* chunk_) { reset(chunk_); }

      void reset(Chunk* chunk_) {
        this->chunk = chunk_;
        i = 0;
      }

      T* next() {
        for (;;) {
          if unlikely (i >= lengthof(chunk->data))
            return null;
          if unlikely (chunk->freemap[i]) {
            i++;
            continue;
          }
          return &chunk->data[i++];
        }
      }
    };

    int getentries(T* a, int limit) {
      Iterator iter(this);
      T* entry;
      int n = 0;
      while ((entry = iter.next())) {
        if unlikely (n >= limit)
          return n;
        a[n++] = *entry;
      }

      return n;
    }
  };

  struct Locator {
    Chunk* chunk;
    int index;

    void reset() {
      chunk = null;
      index = 0;
    }
  };

  selflistlink* head;
  int elemcount;

  ChunkList() {
    head = null;
    elemcount = 0;
  }

  bool add(const T& entry, Locator& hint) {
    Chunk* chunk = (Chunk*)head;

    while (chunk) {
      __builtin_prefetch(chunk->link.next);
      int index = chunk->add(entry);
      if likely (index >= 0) {
        hint.chunk = chunk;
        hint.index = index;
        elemcount++;
        return true;
      }
      chunk = (Chunk*)chunk->link.next;
    }

    Chunk* newchunk = new Chunk();
    newchunk->link.addto(head);

    int index = newchunk->add(entry);
    assert(index >= 0);

    hint.chunk = newchunk;
    hint.index = index;
    elemcount++;

    return true;
  }

  bool remove(const Locator& locator) {
    locator.chunk->remove(locator.index);
    elemcount--;

    if (locator.chunk->empty()) {
      locator.chunk->link.unlink();
      delete locator.chunk;
    }

    return empty();
  }

  void clear() {
    Chunk* chunk = (Chunk*)head;

    while (chunk) {
      Chunk* next = (Chunk*)chunk->link.next;
      __builtin_prefetch(next);
      delete chunk;
      chunk = next;
    }

    elemcount = 0;
    head = null;
  }

  int count() { return elemcount; }

  bool empty() { return (elemcount == 0); }

  ~ChunkList() { clear(); }

  struct Iterator {
    Chunk* chunk;
    Chunk* nextchunk;
    int i;

    Iterator() {}

    Iterator(ChunkList<T, N>* chunklist) { reset(chunklist); }

    void reset(ChunkList<T, N>* chunklist) {
      chunk = (Chunk*)chunklist->head;
      nextchunk = (chunk) ? (Chunk*)chunk->link.next : null;
      i = 0;
    }

    T* next() {
      for (;;) {
        if unlikely (!chunk)
          return null;

        if unlikely (i >= lengthof(chunk->data)) {
          chunk = nextchunk;
          if unlikely (!chunk)
            return null;
          nextchunk = (Chunk*)chunk->link.next;
          __builtin_prefetch(nextchunk);
          i = 0;
        }

        if unlikely (chunk->freemap[i]) {
          i++;
          continue;
        }

        return &chunk->data[i++];
      }
    }
  };

  int getentries(T* a, int limit) {
    Iterator iter(this);
    T* entry;
    int n = 0;
    while ((entry = iter.next())) {
      if unlikely (n >= limit)
        return n;
      a[n++] = *entry;
    }

    return n;
  }
};

// Fast vectorized method: empty only if all slots are literally zero
static inline bool bytes_are_all_zero(const void* ptr, size_t bytes) {
  // Fast vectorized method: empty only if all slots are literally zero
  const W64* p = (const W64*)ptr;
  W64 v = 0;
  foreach (i, (bytes / 8))
    v |= p[i];
  if unlikely (v % 8) {
    v |= (p[(bytes / 8)] & bitmask((bytes % 8) * 8));
  }
  return (v == 0);
}

//
// sort - sort an array of elements
// @p: pointer to data to sort
// @n: number of elements
//
// This function does a heapsort on the given array. You may provide a
// comparison function optimized to your element type.
//
// Sorting time is O(n log n) both on average and worst-case. While
// qsort is about 20% faster on average, it suffers from exploitable
// O(n*n) worst-case behavior and extra memory requirements that make
// it less suitable for kernel use.
//
template<typename T>
struct DefaultComparator {
  int operator()(const T& a, const T& b) const {
    int r = (a < b) ? -1 : +1;
    if (a == b)
      r = 0;
    return r;
  }
};

template<typename T, bool backwards = 0>
struct SortPrecomputedIndexListComparator {
  const T* v;

  SortPrecomputedIndexListComparator(const T* values) : v(values) {}

  int operator()(unsigned long a, unsigned long b) const {
    //
    // This strange construction helps the compiler do better peephole
    // optimization tricks using conditional moves when using integers.
    //
    if (backwards) {
      int r = (v[a] > v[b]) ? -1 : +1;
      if (v[a] == v[b])
        r = 0;
      return r;
    } else {
      int r = (v[a] < v[b]) ? -1 : +1;
      if (v[a] == v[b])
        r = 0;
      return r;
    }
  }
};

template<typename T, bool backwards = 0>
struct PointerSortComparator {
  PointerSortComparator() {}

  int operator()(T* a, T* b) const {
    //
    // This strange construction helps the compiler do better peephole
    // optimization tricks using conditional moves when using integers.
    //
    const T& aa = *a;
    const T& bb = *b;
    if (backwards) {
      int r = (aa > bb) ? -1 : +1;
      if (aa == bb)
        r = 0;
      return r;
    } else {
      int r = (aa < bb) ? -1 : +1;
      if (aa == bb)
        r = 0;
      return r;
    }
  }
};

template<typename T, typename Comparator>
void sort(T* p, size_t n, const Comparator& compare = DefaultComparator<T>()) {
  int c;

  // heapify
  for (int i = (n / 2); i >= 0; i--) {
    for (int r = i; r * 2 < n; r = c) {
      c = r * 2;
      c += ((c < (n - 1)) && (compare(p[c], p[c + 1]) < 0));
      if (compare(p[r], p[c]) >= 0)
        break;
      swap(p[r], p[c]);
    }
  }

  // sort
  for (int i = n - 1; i >= 0; i--) {
    swap(p[0], p[i]);
    for (int r = 0; r * 2 < i; r = c) {
      c = r * 2;
      c += ((c < i - 1) && (compare(p[c], p[c + 1]) < 0));
      if (compare(p[r], p[c]) >= 0)
        break;
      swap(p[r], p[c]);
    }
  }
}

//
// Safe divide and remainder functions that return true iff operation did not generate an exception:
//
template<typename T>
bool div_rem(T& quotient, T& remainder, T dividend_hi, T dividend_lo, T divisor);
template<typename T>
bool div_rem_s(T& quotient, T& remainder, T dividend_hi, T dividend_lo, T divisor);

} // namespace superstl

//
// std::formatter specializations
//
} // namespace x86sim

namespace std {

template<typename K, typename T, typename LM, int setcount, typename KM>
struct formatter<x86sim::superstl::SelfHashtable<K, T, setcount, LM, KM>> {
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

  auto format(const x86sim::superstl::SelfHashtable<K, T, setcount, LM, KM>& ht, format_context& ctx) const {
    auto out = ctx.out();
    out = format_to(out, "Hashtable of {} sets containing {} entries:\n", setcount, ht.count);
    for (int i = 0; i < setcount; i++) {
      x86sim::superstl::selflistlink* tlink = ht.sets[i];
      if (!tlink)
        continue;
      out = format_to(out, "  Set {}:\n", i);
      while (tlink) {
        T* obj = LM::objof(tlink);
        out = format_to(out, "    {} -> {}\n", LM::keyof(obj), *obj);
        tlink = tlink->next;
      }
    }
    return out;
  }
};

template<>
struct formatter<x86sim::superstl::bitstring> {
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

  auto format(const x86sim::superstl::bitstring& bs, format_context& ctx) const {
    auto out = ctx.out();
    if (bs.reverse) {
      for (int i = 0; i < bs.n; i++)
        out = format_to(out, "{}", (char)(((bs.bits >> i) & 1) + '0'));
    } else {
      for (int i = bs.n - 1; i >= 0; i--)
        out = format_to(out, "{}", (char)(((bs.bits >> i) & 1) + '0'));
    }
    return out;
  }
};

template<>
struct formatter<x86sim::superstl::hexstring> {
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

  auto format(const x86sim::superstl::hexstring& hs, format_context& ctx) const {
    int n = (hs.n + 3) / 4;
    const x86sim::W64 mask = (hs.n >= 64) ? ~0ULL : ((1ULL << hs.n) - 1);
    return format_to(ctx.out(), "{:0{}x}", hs.value & mask, n);
  }
};

template<>
struct formatter<x86sim::superstl::bytemaskstring> {
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

  auto format(const x86sim::superstl::bytemaskstring& bs, format_context& ctx) const {
    auto out = ctx.out();
    for (int i = 0; i < bs.n; i++) {
      if ((bs.mask >> i) & 1)
        out = format_to(out, "{:02x}", bs.bytes[i]);
      else
        out = format_to(out, "xx");
      if (((i % bs.splitat) == (bs.splitat - 1)) && (i != bs.n - 1))
        out = format_to(out, "\n");
      else if (i != bs.n - 1)
        out = format_to(out, " ");
    }
    return out;
  }
};

} // namespace std

#endif // _SUPERSTL_H_
