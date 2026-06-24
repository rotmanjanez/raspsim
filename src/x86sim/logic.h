// -*- c++ -*-
//
// Sequential Logic Primitives for C++
//
// Copyright 1999-2008 Matt T. Yourst <yourst@yourst.com>
//

#ifndef _LOGIC_H_
#define _LOGIC_H_

#include <array>
#include "globals.h"
#include "superstl.h"

namespace x86sim {

using namespace superstl;

//
// Queue
//

// Iterate forward through queue from head to tail
#define foreach_forward(Q, i) for (int i = (Q).head; i != (Q).tail; i = add_index_modulo(i, +1, (Q).size))

// Iterate forward through queue from the specified entry until the tail
#define foreach_forward_from(Q, E, i) for (int i = E->index(); i != (Q).tail; i = add_index_modulo(i, +1, (Q).size))

// Iterate forward through queue from the entry after the specified entry until the tail
#define foreach_forward_after(Q, E, i)                                                                                 \
  for (int i = add_index_modulo(E->index(), +1, (Q).size); i != (Q).tail; i = add_index_modulo(i, +1, (Q).size))

// Iterate backward through queue from tail to head
#define foreach_backward(Q, i)                                                                                         \
  for (int i = add_index_modulo((Q).tail, -1, (Q).size); i != add_index_modulo((Q).head, -1, (Q).size);                \
       i = add_index_modulo(i, -1, (Q).size))

// Iterate backward through queue from the specified entry until the tail
#define foreach_backward_from(Q, E, i)                                                                                 \
  for (int i = E->index(); i != add_index_modulo((Q).head, -1, (Q).size); i = add_index_modulo(i, -1, (Q).size))

// Iterate backward through queue from the entry before the specified entry until the head
#define foreach_backward_before(Q, E, i)                                                                               \
  for (int i = add_index_modulo(E->index(), -1, (Q).size);                                                             \
       ((i != add_index_modulo((Q).head, -1, (Q).size)) && (E->index() != (Q).head));                                  \
       i = add_index_modulo(i, -1, (Q).size))

template<class T, int SIZE>
struct FixedQueue : public std::array<T, SIZE> {
  int head;  // used for allocation
  int tail;  // used for deallocation
  int count; // count of entries

  static const int size = SIZE;

  FixedQueue() { reset(); }

  void flush() { head = tail = count = 0; }

  void reset() { head = tail = count = 0; }

  int remaining() const { return std::max((SIZE - count) - 1, 0); }

  bool empty() const { return (!count); }

  bool full() const { return (!remaining()); }

  T* alloc() {
    if (!remaining())
      return null;

    T* entry = &(*this)[tail];

    tail = add_index_modulo(tail, +1, SIZE);
    count++;

    return entry;
  }

  T* push() { return alloc(); }

  T* push(const T& data) {
    T* slot = push();
    if (!slot)
      return null;
    *slot = data;
    return slot;
  }

  T* enqueue(const T& data) { return push(data); }

  void commit(T& entry) {
    assert(entry.index() == head);
    count--;
    head = add_index_modulo(head, +1, SIZE);
  }

  void annul(T& entry) {
    // assert(entry.index() == add_index_modulo(tail, -1, SIZE));
    count--;
    tail = add_index_modulo(tail, -1, SIZE);
  }

  T* pop() {
    if (empty())
      return null;
    tail = add_index_modulo(tail, -1, SIZE);
    count--;
    return &(*this)[tail];
  }

  T* peek() {
    if (empty())
      return null;
    return &(*this)[head];
  }

  T* dequeue() {
    if (empty())
      return null;
    count--;
    T* entry = &(*this)[head];
    head = add_index_modulo(head, +1, SIZE);
    return entry;
  }

  void commit(T* entry) { commit(*entry); }
  void annul(T* entry) { annul(*entry); }

  T* pushhead() {
    if (full())
      return null;
    head = add_index_modulo(head, -1, SIZE);
    count++;
    return &(*this)[head];
  }

  T* pophead() {
    if (empty())
      return null;
    T* p = &(*this)[head];
    count--;
    head = add_index_modulo(head, +1, SIZE);
    return p;
  }

  T* peekhead() {
    if (empty())
      return null;
    return &(*this)[head];
  }

  T* peektail() {
    if (empty())
      return null;
    int t = add_index_modulo(tail, -1, SIZE);
    return &(*this)[t];
  }

  T& operator()(int index) {
    index = add_index_modulo(head, index, SIZE);
    return (*this)[index];
  }

  const T& operator()(int index) const {
    index = add_index_modulo(head, index, SIZE);
    return (*this)[index];
  }
};

} // namespace x86sim

template<class T, int SIZE>
struct std::formatter<x86sim::FixedQueue<T, SIZE>> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const x86sim::FixedQueue<T, SIZE>& q, std::format_context& ctx) const {
    using namespace x86sim;
    auto out = ctx.out();
    out = std::format_to(out, "Queue<{}>: head {} to tail {} ({} entries):\n", SIZE, q.head, q.tail, q.count);

    for (int i = q.head; i != q.tail; i = add_index_modulo(i, +1, SIZE)) {
      const T& entry = q[i];
      out = std::format_to(out, "  slot {:>3}: {}\n", i, entry);
    }

    return out;
  }
};

namespace x86sim {

template<class T, int SIZE>
struct Queue : public FixedQueue<T, SIZE> {
  typedef FixedQueue<T, SIZE> base_t;

  Queue() { reset(); }

  void reset() {
    base_t::reset();
    foreach (i, SIZE) {
      (*this)[i].init(i);
    }
  }

  T* alloc() {
    T* p = base_t::alloc();
    if likely (p)
      p->validate();
    return p;
  }
};

} // namespace x86sim

template<class T, int Size>
struct std::formatter<x86sim::Queue<T, Size>> : std::formatter<x86sim::FixedQueue<T, Size>> {
  // Inherits parse and format from x86sim::FixedQueue formatter
};

namespace x86sim {

//
// Fully Associative Arrays
//

template<typename T>
struct InvalidTag {
  static const T INVALID;
};
template<>
struct InvalidTag<W64> {
  static const W64 INVALID = 0xffffffffffffffffULL;
};
template<>
struct InvalidTag<W32> {
  static const W32 INVALID = 0xffffffff;
};
template<>
struct InvalidTag<W16> {
  static const W16 INVALID = 0xffff;
};
template<>
struct InvalidTag<W8> {
  static const W8 INVALID = 0xff;
};

//
// The replacement policy is pseudo-LRU using a most recently used
// bit vector (mLRU), as described in the paper "Performance Evaluation
// of Cache Replacement Policies for the SPEC CPU2000 Benchmark Suite"
// by Al-Zoubi et al. Essentially we maintain one MRU bit per way and
// set the bit for the way when that way is accessed. The way to evict
// is the first way without its MRU bit set. If all MRU bits become
// set, they are all reset and we start over. Surprisingly, this
// simple method performs as good as, if not better than, true LRU
// or tree-based hot sector LRU.
//

template<typename T, int ways>
struct FullyAssociativeTags {
  std::bitset<ways> evictmap;
  T tags[ways];

  static const T INVALID = InvalidTag<T>::INVALID;

  FullyAssociativeTags() { reset(); }

  void reset() {
    evictmap = 0;
    foreach (i, ways) {
      tags[i] = INVALID;
    }
  }

  void use(int way) {
    evictmap[way] = 1;
    // Performance is somewhat better with this off with higher associativity caches:
    // if (evictmap.all()) evictmap = 0;
  }

  //
  // This is a clever way of doing branch-free matching
  // with conditional moves and addition. It relies on
  // having at most one matching entry in the array;
  // otherwise the algorithm breaks:
  //
  int match(T target) {
    int way = 0;
    foreach (i, ways) {
      way += (tags[i] == target) ? (i + 1) : 0;
    }

    return way - 1;
  }

  int probe(T target) {
    int way = match(target);
    if (way < 0)
      return -1;

    use(way);
    return way;
  }

  int lru() const { return (evictmap.all()) ? 0 : lsb(~evictmap); }

  int select(T target, T& oldtag) {
    int way = probe(target);
    if (way < 0) {
      way = lru();
      if (evictmap.all())
        evictmap = 0;
      oldtag = tags[way];
      tags[way] = target;
    }
    use(way);
    return way;
  }

  int select(T target) {
    T dummy;
    return select(target, dummy);
  }

  void invalidate_way(int way) {
    tags[way] = INVALID;
    evictmap[way] = 0;
  }

  int invalidate(T target) {
    int way = probe(target);
    if (way < 0)
      return -1;
    invalidate_way(way);
    return way;
  }

  const T& operator[](int index) const { return tags[index]; }

  T& operator[](int index) { return tags[index]; }
  int operator()(T target) { return probe(target); }

};

} // namespace x86sim

template<typename T, int ways>
struct std::formatter<x86sim::FullyAssociativeTags<T, ways>> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const x86sim::FullyAssociativeTags<T, ways>& tags, std::format_context& ctx) const {
    using namespace x86sim;
    auto out = ctx.out();
    for (int i = 0; i < ways; i++) {
      out = std::format_to(out, "  way {:<2}: ", i);
      if (tags[i] != tags.INVALID) {
        out = std::format_to(out, "tag 0x{}", hexstring(tags[i], sizeof(T) * 8));
        if (tags.evictmap[i])
          out = std::format_to(out, " (MRU)");
      } else {
        out = std::format_to(out, "<invalid>");
      }
      out = std::format_to(out, "\n");
    }
    return out;
  }
};

namespace x86sim {


//
// Associative array
//
// Limitations:
//
// - Every tag in the array must be unique,
//   except for the invalid tag (all 1s)
//
// - <size> can be from 1 to 128.
//
// - <width> in bits can be from 1 to 64
//

template<int size, int width, int padsize = 0>
struct FullyAssociativeTagsNbitOneHot {
  typedef W64 base_t;

  static const int slices = (width + 7) / 8;
  static const int storewidth = slices * 8;
  static constexpr base_t invalid_tag = ~base_t{0};
  static constexpr base_t storemask = (storewidth >= 64) ? invalid_tag : ((base_t{1} << storewidth) - 1);

  std::array<base_t, size + padsize> tags;
  std::bitset<size> valid;
  std::bitset<size> evictmap;

  FullyAssociativeTagsNbitOneHot() { reset(); }

  void reset() {
    valid = 0;
    evictmap = 0;
    tags.fill(invalid_tag);
  }

  static base_t stored_tag(base_t tag) { return tag & storemask; }

  int match(base_t tag) const {
    const base_t target = stored_tag(tag);
    foreach (i, size) {
      if (valid[i] && (stored_tag(tags[i]) == target))
        return i;
    }
    return -1;
  }

  int search(base_t tag) const { return match(tag); }

  int operator()(base_t tag) const { return search(tag); }

  void update(int index, base_t tag) {
    tags[index] = tag;
    valid[index] = 1;
    evictmap[index] = 1;
  }

  class ref {
    friend class FullyAssociativeTagsNbitOneHot;

    FullyAssociativeTagsNbitOneHot<size, width, padsize>& tags;
    int index;

    ref();

  public:
    inline ref(FullyAssociativeTagsNbitOneHot& tags_, int index_) : tags(tags_), index(index_) {}

    inline ~ref() {}

    inline ref& operator=(base_t tag) {
      tags.update(index, tag);
      return *this;
    }

    inline ref& operator=(const ref& other) {
      tags.update(index, other.tags[other.index]);
      return *this;
    }
  };

  friend class ref;

  ref operator[](int index) { return ref(*this, index); }
  base_t operator[](int index) const { return tags[index]; }

  bool isvalid(int index) { return valid[index]; }

  int insertslot(int idx, base_t tag) {
    valid[idx] = 1;
    (*this)[idx] = tag;
    return idx;
  }

  int insert(base_t tag) {
    foreach (idx, size) {
      if (!valid[idx])
        return insertslot(idx, tag);
    }
    return -1;
  }

  void invalidateslot(int index) {
    valid[index] = 0;
    tags[index] = invalid_tag;
  }

  void validateslot(int index) { valid[index] = 1; }

  int invalidate(base_t target) {
    int index = match(target);
    if (index < 0)
      return 0;
    invalidateslot(index);
    return 1;
  }

  std::bitset<size> masked_match(base_t targettag, base_t tagmask) {
    std::bitset<size> m;

    foreach (i, size) {
      base_t tag = tags[i];
      m[i] = valid[i] && ((tag & tagmask) == targettag);
    }

    return m;
  }

  void masked_invalidate(const std::bitset<size>& slotmask) {
    foreach (i, size) {
      if unlikely (slotmask[i])
        invalidateslot(i);
    }
  }

  void use(int way) {
    evictmap[way] = 1;

    if (evictmap.all()) {
      evictmap = 0;
      evictmap[way] = 1;
    }
  }

  int probe(base_t target) {
    int way = match(target);
    if (way < 0)
      return way;
    use(way);
    return way;
  }

  int lru() const {
    foreach (i, size) {
      if (!evictmap[i])
        return i;
    }
    return 0;
  }

  int select(base_t target, base_t& oldtag) {
    int way = probe(target);
    if (way < 0) {
      way = lru();
      if (evictmap.all())
        evictmap = 0;
      oldtag = tags[way];
      update(way, target);
    }
    use(way);
    return way;
  }

  int select(base_t target) {
    base_t dummy;
    return select(target, dummy);
  }

  std::string slotid(int slot) const {
    base_t tag = (*this)[slot];
    std::string result = std::format("{:>3}: {:016x} ", slot, tag);
    for (int i = 0; i < slices; i++) {
      const byte b = (byte)(tag >> (i * 8));
      result += std::format(" {:02x}", b);
    }
    if (!valid[slot])
      result += " <invalid>";
    return result;
  }
};

} // namespace x86sim

template<int Size, int Width, int PadSize>
struct std::formatter<x86sim::FullyAssociativeTagsNbitOneHot<Size, Width, PadSize>> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const x86sim::FullyAssociativeTagsNbitOneHot<Size, Width, PadSize>& tags, std::format_context& ctx) const {
    using namespace x86sim;
    auto out = ctx.out();
    for (int i = 0; i < Size; i++) {
      out = std::format_to(out, "{}\n", tags.slotid(i));
    }
    return out;
  }
};

namespace x86sim {

template<typename T, typename V>
struct NullAssociativeArrayStatisticsCollector {
  static void inserted(V& elem, T newtag, int way) {}
  static void replaced(V& elem, T oldtag, T newtag, int way) {}
  static void probed(V& elem, T tag, int way, bool hit) {}
  static void overflow(T tag) {}
  static void locked(V& slot, T tag, int way) {}
  static void unlocked(V& slot, T tag, int way) {}
  static void invalidated(V& elem, T oldtag, int way) {}
};

template<typename T, typename V, int ways, typename stats = NullAssociativeArrayStatisticsCollector<T, V>>
struct FullyAssociativeArray {
  FullyAssociativeTags<T, ways> tags;
  V data[ways];

  FullyAssociativeArray() { reset(); }

  void reset() {
    tags.reset();
    foreach (i, ways) {
      data[i].reset();
    }
  }

  V* probe(T tag) {
    int way = tags.probe(tag);
    stats::probed((way < 0) ? data[0] : data[way], tag, way, (way >= 0));
    return (way < 0) ? null : &data[way];
  }

  V* select(T tag, T& oldtag) {
    int way = tags.select(tag, oldtag);

    V& slot = data[way];

    if ((way >= 0) & (tag == oldtag)) {
      stats::probed(slot, tag, way, 1);
    } else {
      if (oldtag == tags.INVALID)
        stats::inserted(slot, tag, way);
      else
        stats::replaced(slot, oldtag, tag, way);
    }

    return &slot;
  }

  V* select(T tag) {
    T dummy;
    return select(tag, dummy);
  }

  int wayof(const V* line) const {
    int way = (line - (const V*)&data);
#if 0
    assert(inrange(way, 0, ways-1));
#endif
    return way;
  }

  T tagof(V* line) {
    int way = wayof(line);
    return tags.tags[way];
  }

  void invalidate_way(int way) {
    stats::invalidated(data[way], tags[way], way);
    tags.invalidate_way(way);
    data[way].reset();
  }

  void invalidate_line(V* line) { invalidate_way(wayof(line)); }

  int invalidate(T tag) {
    int way = tags.probe(tag);
    if (way < 0)
      return -1;
    invalidate_way(way);
    return way;
  }

  V& operator[](int way) { return data[way]; }

  V* operator()(T tag) { return select(tag); }
};

} // namespace x86sim

template<typename T, typename V, int ways, typename stats>
struct std::formatter<x86sim::FullyAssociativeArray<T, V, ways, stats>> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const x86sim::FullyAssociativeArray<T, V, ways, stats>& arr, std::format_context& ctx) const {
    using namespace x86sim;
    auto out = ctx.out();
    for (int i = 0; i < ways; i++) {
      std::string waystr = std::format("  way {:<2}: ", i);
      if (arr.tags[i] != arr.tags.INVALID) {
        waystr += std::format("tag 0x{}", hexstring(arr.tags[i], sizeof(T) * 8));
        if (arr.tags.evictmap[i])
          waystr += " (MRU)";
      } else {
        waystr += "<invalid>";
      }
      out = std::format_to(out, "{:<40} -> ", waystr);
      out = std::format_to(out, "{}\n", arr.data[i]);
    }
    return out;
  }
};

namespace x86sim {

template<typename T, typename V, int setcount, int waycount, int linesize,
         typename stats = NullAssociativeArrayStatisticsCollector<T, V>>
struct AssociativeArray {
  typedef FullyAssociativeArray<T, V, waycount, stats> Set;
  Set sets[setcount];

  AssociativeArray() { reset(); }

  void reset() {
    foreach (set, setcount) {
      sets[set].reset();
    }
  }

  static int setof(T addr) { return bits(addr, log2(linesize), log2(setcount)); }

  static T tagof(T addr) { return floor(addr, linesize); }

  V* probe(T addr) { return sets[setof(addr)].probe(tagof(addr)); }

  V* select(T addr, T& oldaddr) { return sets[setof(addr)].select(tagof(addr), oldaddr); }

  V* select(T addr) {
    T dummy;
    return sets[setof(addr)].select(tagof(addr), dummy);
  }

  void invalidate(T addr) { sets[setof(addr)].invalidate(tagof(addr)); }
};

} // namespace x86sim

template<typename T, typename V, int setcount, int waycount, int linesize, typename stats>
struct std::formatter<x86sim::AssociativeArray<T, V, setcount, waycount, linesize, stats>> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const x86sim::AssociativeArray<T, V, setcount, waycount, linesize, stats>& arr, std::format_context& ctx) const {
    using namespace x86sim;
    auto out = ctx.out();
    out = std::format_to(out, "x86sim::AssociativeArray<{} sets, {} ways, {}-byte lines>:\n", setcount, waycount, linesize);
    for (int set = 0; set < setcount; set++) {
      out = std::format_to(out, "  Set {}:\n", set);
      out = std::format_to(out, "{}", arr.sets[set]);
    }
    return out;
  }
};

namespace x86sim {

//
// Lockable version of associative arrays:
//

template<typename T, int ways>
struct LockableFullyAssociativeTags {
  std::bitset<ways> evictmap;
  std::bitset<ways> unlockedmap;
  T tags[ways];

  static const T INVALID = InvalidTag<T>::INVALID;

  LockableFullyAssociativeTags() { reset(); }

  void reset() {
    evictmap = 0;
    unlockedmap.set();
    foreach (i, ways) {
      tags[i] = INVALID;
    }
  }

  void use(int way) {
    evictmap[way] = 1;
    // Performance is somewhat better with this off with higher associativity caches:
    // if (evictmap.all()) evictmap = 0;
  }

  //
  // This is a clever way of doing branch-free matching
  // with conditional moves and addition. It relies on
  // having at most one matching entry in the array;
  // otherwise the algorithm breaks:
  //
  int match(T target) {
    int way = 0;
    foreach (i, ways) {
      way += (tags[i] == target) ? (i + 1) : 0;
    }

    return way - 1;
  }

  int probe(T target) {
    int way = match(target);
    if (way < 0)
      return -1;

    use(way);
    return way;
  }

  int lru() const {
    if (unlockedmap.none())
      return -1;
    auto w = (~evictmap) & unlockedmap;
    return w.any() ? lsb(w) : 0;
  }

  int select(T target, T& oldtag) {
    int way = probe(target);
    if (way < 0) {
      way = lru();
      if (way < 0)
        return -1;
      if (evictmap.all())
        evictmap = 0;
      oldtag = tags[way];
      tags[way] = target;
    }
    use(way);
    return way;
  }

  int select(T target) {
    T dummy;
    return select(target, dummy);
  }

  int select_and_lock(T tag, bool& firstlock, T& oldtag) {
    int way = select(tag, oldtag);
    if (way < 0)
      return way;
    firstlock = unlockedmap[way];
    lock(way);
    return way;
  }

  int select_and_lock(T tag, bool& firstlock) {
    T dummy;
    return select_and_lock(tag, firstlock, dummy);
  }

  int select_and_lock(T target) {
    bool dummy;
    return select_and_lock(target, dummy);
  }

  void invalidate_way(int way) {
    tags[way] = INVALID;
    evictmap[way] = 0;
    unlockedmap[way] = 1;
  }

  int invalidate(T target) {
    int way = probe(target);
    if (way < 0)
      return -1;
    invalidate_way(way);
  }

  bool islocked(int way) const { return !unlockedmap[way]; }

  void lock(int way) { unlockedmap[way] = 0; }
  void unlock(int way) { unlockedmap[way] = 1; }

  const T& operator[](int index) const { return tags[index]; }

  T& operator[](int index) { return tags[index]; }
  int operator()(T target) { return probe(target); }

};

} // namespace x86sim

template<typename T, int ways>
struct std::formatter<x86sim::LockableFullyAssociativeTags<T, ways>> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const x86sim::LockableFullyAssociativeTags<T, ways>& tags, std::format_context& ctx) const {
    using namespace x86sim;
    auto out = ctx.out();
    for (int i = 0; i < ways; i++) {
      out = std::format_to(out, "  way {:<2}: ", i);
      if (tags[i] != tags.INVALID) {
        out = std::format_to(out, "tag 0x{}", hexstring(tags[i], sizeof(T) * 8));
        if (tags.evictmap[i])
          out = std::format_to(out, " (MRU)");
        if (!tags.unlockedmap[i])
          out = std::format_to(out, " (locked)");
      } else {
        out = std::format_to(out, "<invalid>");
      }
      out = std::format_to(out, "\n");
    }
    return out;
  }
};

namespace x86sim {

template<typename T, typename V, int ways, typename stats = NullAssociativeArrayStatisticsCollector<T, V>>
struct LockableFullyAssociativeArray {
  LockableFullyAssociativeTags<T, ways> tags;
  V data[ways];

  LockableFullyAssociativeArray() { reset(); }

  void reset() {
    tags.reset();
    foreach (i, ways) {
      data[i].reset();
    }
  }

  V* probe(T tag) {
    int way = tags.probe(tag);
    stats::probed((way < 0) ? data[0] : data[way], tag, way, (way >= 0));
    return (way < 0) ? null : &data[way];
  }

  V* select(T tag, T& oldtag) {
    int way = tags.select(tag, oldtag);

    if (way < 0) {
      stats::overflow(tag);
      return null;
    }

    V& slot = data[way];

    if ((way >= 0) & (tag == oldtag)) {
      stats::probed(slot, tag, way, 1);
    } else {
      if (oldtag == tags.INVALID)
        stats::inserted(slot, tag, way);
      else
        stats::replaced(slot, oldtag, tag, way);
    }

    return &slot;
  }

  V* select(T tag) {
    T dummy;
    return select(tag, dummy);
  }

  V* select_and_lock(T tag, bool& firstlock, T& oldtag) {
    int way = tags.select_and_lock(tag, firstlock, oldtag);

    if (way < 0) {
      stats::overflow(tag);
      return null;
    }

    V& slot = data[way];

    if (tag == oldtag) {
      stats::probed(slot, tag, way, 1);
    } else {
      if (oldtag == tags.INVALID)
        stats::inserted(slot, tag, way);
      else
        stats::replaced(slot, oldtag, tag, way);
      stats::locked(slot, tag, way);
    }

    return &slot;
  }

  V* select_and_lock(T tag, bool& firstlock) {
    T dummy;
    return select_and_lock(tag, firstlock, dummy);
  }

  V* select_and_lock(T tag) {
    bool dummy;
    return select_and_lock(tag, dummy);
  }

  int wayof(const V* line) const {
    int way = (line - (const V*)&data);
#if 0
    assert(inrange(way, 0, ways-1));
#endif
    return way;
  }

  T tagof(V* line) {
    int way = wayof(line);
    return tags.tags[way];
  }

  void invalidate_way(int way) {
    unlock_way(way);
    stats::invalidated(data[way], tags[way], way);
    tags.invalidate_way(way);
    data[way].reset();
  }

  void invalidate_line(V* line) { invalidate_way(wayof(line)); }

  int invalidate(T tag) {
    int way = tags.probe(tag);
    if (way < 0)
      return -1;
    invalidate_way(way);
    return way;
  }

  void unlock_way(int way) {
    stats::unlocked(data[way], tags[way], way);
    tags.unlock(way);
  }

  void unlock_line(V* line) { unlock_way(wayof(line)); }

  int unlock(T tag) {
    int way = tags.probe(tag);
    if (way < 0)
      return 0;
    unlock_way(way);
    if (tags.islocked(way))
      stats::unlocked(data[way], tags[way], way);
    return way;
  }

  V& operator[](int way) { return data[way]; }

  V* operator()(T tag) { return select(tag); }
};

} // namespace x86sim

template<typename T, typename V, int ways, typename stats>
struct std::formatter<x86sim::LockableFullyAssociativeArray<T, V, ways, stats>> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const x86sim::LockableFullyAssociativeArray<T, V, ways, stats>& arr, std::format_context& ctx) const {
    using namespace x86sim;
    auto out = ctx.out();
    for (int i = 0; i < ways; i++) {
      std::string waystr = std::format("  way {:<2}: ", i);
      if (arr.tags[i] != arr.tags.INVALID) {
        waystr += std::format("tag 0x{}", hexstring(arr.tags[i], sizeof(T) * 8));
        if (arr.tags.evictmap[i])
          waystr += " (MRU)";
        if (!arr.tags.unlockedmap[i])
          waystr += " (locked)";
      } else {
        waystr += "<invalid>";
      }
      out = std::format_to(out, "{:<40} -> ", waystr);
      out = std::format_to(out, "{}\n", arr.data[i]);
    }
    return out;
  }
};

namespace x86sim {

template<typename T, typename V, int setcount, int waycount, int linesize,
         typename stats = NullAssociativeArrayStatisticsCollector<T, V>>
struct LockableAssociativeArray {
  typedef LockableFullyAssociativeArray<T, V, waycount, stats> Set;
  Set sets[setcount];

  LockableAssociativeArray() { reset(); }

  void reset() {
    foreach (set, setcount) {
      sets[set].reset();
    }
  }

  static int setof(T addr) { return bits(addr, log2(linesize), log2(setcount)); }

  static T tagof(T addr) { return floor(addr, linesize); }

  V* probe(T addr) { return sets[setof(addr)].probe(tagof(addr)); }

  V* select(T addr, T& oldaddr) { return sets[setof(addr)].select(tagof(addr), oldaddr); }

  V* select(T addr) {
    T dummy;
    return select(addr, dummy);
  }

  void invalidate(T addr) { sets[setof(addr)].invalidate(tagof(addr)); }

  V* select_and_lock(T addr, bool& firstlock, T& oldtag) {
    V* line = sets[setof(addr)].select_and_lock(tagof(addr), firstlock, oldtag);
    return line;
  }

  V* select_and_lock(T addr, bool& firstlock) {
    W64 dummy;
    return select_and_lock(addr, firstlock, dummy);
  }

  V* select_and_lock(T addr) {
    bool dummy;
    return select_and_lock(addr, dummy);
  }
};

} // namespace x86sim

template<typename T, typename V, int setcount, int waycount, int linesize, typename stats>
struct std::formatter<x86sim::LockableAssociativeArray<T, V, setcount, waycount, linesize, stats>> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const x86sim::LockableAssociativeArray<T, V, setcount, waycount, linesize, stats>& arr,
              std::format_context& ctx) const {
    using namespace x86sim;
    auto out = ctx.out();
    out = std::format_to(out, "x86sim::LockableAssociativeArray<{} sets, {} ways, {}-byte lines>:\n", setcount, waycount,
                         linesize);
    for (int set = 0; set < setcount; set++) {
      out = std::format_to(out, "  Set {}:\n", set);
      out = std::format_to(out, "{}", arr.sets[set]);
    }
    return out;
  }
};

namespace x86sim {

template<int size, int width, int padsize = 0>
struct FullyAssociativeTagsNbit {
  static_assert((width == 8) || (width == 16), "FullyAssociativeTagsNbit supports 8-bit and 16-bit tags");

  using base_t = std::conditional_t<width == 8, byte, W16>;
  static constexpr int tags_per_chunk = 16 / sizeof(base_t);
  using tag_array_t = std::array<base_t, tags_per_chunk>;
  using vec_t = tag_array_t;

  std::array<base_t, size + padsize> tags;
  std::bitset<size> valid;

  W64 getvalid() { return valid.to_ullong(); }

  FullyAssociativeTagsNbit() { reset(); }

  base_t operator[](int i) const { return tags[i]; }

  base_t& operator[](int i) { return tags[i]; }

  bool isvalid(int index) { return valid[index]; }

  void reset() {
    valid = 0;
    tags.fill((base_t)-1);
  }

  static tag_array_t prep(base_t tag) {
    tag_array_t result;
    result.fill(tag);
    return result;
  }

  int insertslot(int idx, base_t tag) {
    valid[idx] = 1;
    (*this)[idx] = tag;
    return idx;
  }

  int insert(base_t tag) {
    if (valid.all())
      return -1;
    int idx = lsb(~valid);
    return insertslot(idx, tag);
  }

  std::bitset<size> match(const tag_array_t& target) const {
    std::bitset<size> m = 0;

    foreach (i, size)
      if (tags[i] == target[i % tags_per_chunk])
        m[i] = 1;

    return m & valid;
  }

  std::bitset<size> match(base_t target) const { return match(prep(target)); }

  std::bitset<size> matchany(const tag_array_t& target) const {
    std::bitset<size> m = 0;

    foreach (i, size)
      if ((tags[i] & target[i % tags_per_chunk]) != 0)
        m[i] = 1;

    return m & valid;
  }

  std::bitset<size> matchany(base_t target) const { return matchany(prep(target)); }

  int search(const tag_array_t& target) const {
    std::bitset<size> bitmap = match(target);
    int idx = lsb(bitmap);
    if (bitmap.none())
      idx = -1;
    return idx;
  }

  int extract(const tag_array_t& target) {
    int idx = search(target);
    if (idx >= 0)
      valid[idx] = 0;
    return idx;
  }

  int search(base_t tag) const { return search(prep(tag)); }

  std::bitset<size> extract(base_t tag) { return extract(prep(tag)); }

  void invalidateslot(int index) { valid[index] = 0; }

  const std::bitset<size>& invalidatemask(const std::bitset<size>& mask) {
    valid &= ~mask;
    return mask;
  }

  std::bitset<size> invalidate(const tag_array_t& target) { return invalidatemask(match(target)); }

  std::bitset<size> invalidate(base_t target) { return invalidate(prep(target)); }

  void collapse(int index) {
    for (int i = index; i < (size + padsize) - 1; i++)
      tags[i] = tags[i + 1];
    tags[(size + padsize) - 1] = (base_t)-1;

    for (int i = index; i < size - 1; i++)
      valid[i] = valid[i + 1];
    valid.reset(size - 1);
  }

  void decrement(base_t amount = 1) {
    foreach (i, size)
      tags[i] = (tags[i] < amount) ? 0 : tags[i] - amount;
  }

  void increment(base_t amount = 1) {
    const base_t max = (base_t)-1;
    foreach (i, size)
      tags[i] = (tags[i] > max - amount) ? max : tags[i] + amount;
  }

  std::string slotid(int slot) const {
    int tag = (*this)[slot];
    if (valid[slot])
      return std::format("{:>3}", tag);
    else
      return "???";
  }
};

template<int size, int padsize = 0>
using FullyAssociativeTags8bit = FullyAssociativeTagsNbit<size, 8, padsize>;

template<int size, int padsize = 0>
using FullyAssociativeTags16bit = FullyAssociativeTagsNbit<size, 16, padsize>;

} // namespace x86sim

template<int Size, int Width, int PadSize>
struct std::formatter<x86sim::FullyAssociativeTagsNbit<Size, Width, PadSize>> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const x86sim::FullyAssociativeTagsNbit<Size, Width, PadSize>& tags, std::format_context& ctx) const {
    using namespace x86sim;
    auto out = ctx.out();
    for (int i = 0; i < Size; i++) {
      out = std::format_to(out, "{} ", tags.slotid(i));
    }
    return out;
  }
};

#endif // _LOGIC_H_
