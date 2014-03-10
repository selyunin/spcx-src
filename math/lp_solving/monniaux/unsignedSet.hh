#include <assert.h>
#include <iostream>

class unsignedSet {
  const unsigned maxSize;
  bool * const tab;

public:
  friend class iterator;

  class const_iterator {
    const unsignedSet& set;
    unsigned index;

    const_iterator(const unsignedSet& set0, unsigned index0) :
      set(set0), index(index0) {
    }

    void advance() {
      while(index < set.maxSize && !set.tab[index]) index++;
    }

  public:
    bool operator!=(const const_iterator& x) const {
      return index!=x.index;
    }

    bool operator==(const const_iterator& x) const {
      return index==x.index;
    }

    const_iterator& operator++() { // prefix
      index++;
      advance();
      return *this;
    }

    const_iterator operator++(int) { //postfix
      const_iterator r=*this;
      index++;
      advance();
      return r;
    }

    unsigned operator*() {
      return index;
    }

    friend class unsignedSet;
  };

  class iterator : public const_iterator {
    friend class unsignedSet;

    iterator(unsignedSet& set, unsigned index) : const_iterator(set, index) {
    }

  public:
    iterator& operator++() { // prefix
      index++;
      advance();
      return *this;
    }

    iterator operator++(int) { //postfix
      iterator r=*this;
      index++;
      advance();
      return r;
    }
  };

  unsignedSet(unsigned maxSize0) :
    maxSize(maxSize0), tab(new bool[maxSize0]) {
    for(unsigned i=0; i<maxSize0; i++) {
      tab[i] = false;
    }
  }

  ~unsignedSet() {
    delete[] tab;
  }

  const_iterator find(unsigned i) const {
    return const_iterator(*this, tab[i] ? i : maxSize);
  }

  void insert(unsigned i) {
    assert(i < maxSize);
    tab[i] = true;
  }

  void erase(unsigned i) {
    assert(i < maxSize);
    tab[i] = false;
  }

  void insert(iterator i) {
    assert(*i < maxSize);
    tab[*i] = true;
  }

  void erase(iterator i) {
    assert(*i < maxSize);
    tab[*i] = false;
  }

  iterator begin() {
    iterator r(*this, 0);
    r.advance();
    return r;
  }

  iterator end() {
    return iterator(*this, maxSize);
  }

  const_iterator begin() const {
    const_iterator r(*this, 0);
    r.advance();
    return r;
  }

  const_iterator end() const {
    return const_iterator(*this, maxSize);
  }

  unsigned size() const {
    unsigned count=0;
    for(unsigned i=0; i<maxSize; i++) {
      if (tab[i]) count++;
    }
    return count;
  }

  void clear() {
    for(unsigned i=0; i<maxSize; i++) tab[i]=false;
  }
};
