#include <stdlib.h>
#include <malloc.h>

typedef double double2 __attribute__((vector_size(16)));

template <> struct WithInfinitesimal<double> {
  union {
    struct {
      double standardValue, epsilonCoefficient;
    };
    double2 vec;
  };

  WithInfinitesimal() {
  }

  WithInfinitesimal(const double& v) :
    standardValue(v), epsilonCoefficient(0) {
  }

  WithInfinitesimal(const double& v, const double& eps) :
    standardValue(v), epsilonCoefficient(eps) {
  }

private:
  WithInfinitesimal(const double2& vec0) :
    vec(vec0) {
  }

public:
  WithInfinitesimal<double> operator+(const WithInfinitesimal <double>& v) const {
    return WithInfinitesimal<double>(vec + v.vec);
  }

  WithInfinitesimal<double> operator-(const WithInfinitesimal <double>& v) const {
    return WithInfinitesimal<double>(vec - v.vec);
  }

  WithInfinitesimal<double> operator*(const double s) const {
    return WithInfinitesimal<double>(vec * ((double2) {s, s}));
  }

  WithInfinitesimal<double> operator/(const double s) const {
    return WithInfinitesimal<double>(vec / ((double2) {s, s}));
  }

  WithInfinitesimal<double>& operator*=(const double s) {
    vec *= (double2) {s, s};
    return *this;
  }

  WithInfinitesimal<double>& operator+=(const WithInfinitesimal<double>& s) {
    vec += s.vec;
    return *this;
  }

  WithInfinitesimal<double>& operator-=(const WithInfinitesimal<double>& s) {
    vec -= s.vec;
    return *this;
  }

  WithInfinitesimal<double> operator-() const {
    return WithInfinitesimal<double>(-vec);
  }

  const WithInfinitesimal<double>& getFiniteValue() const {
    return *this;
  }

/* Workaround for what I consider a bug. */
  void *operator new(size_t size) {
    return memalign(__alignof__(WithInfinitesimal<double>), size);
  }

  void operator delete(void *block) {
    free(block);
  }

  void *operator new[](size_t size) {
    return memalign(__alignof__(WithInfinitesimal<double>), size);
  }

  void operator delete[](void *block) {
    free(block);
  }
};

namespace std {
template<> struct allocator<WithInfinitesimal<double> > {
  typedef WithInfinitesimal<double> value_type;
  typedef value_type* pointer;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;
  typedef const value_type* const_pointer;
  typedef value_type& reference;
  typedef const value_type& const_reference;

  pointer allocate(size_type n, const void* = 0) {
    return (pointer) memalign(__alignof__(value_type), n * sizeof(value_type));
  }

  void deallocate(pointer p, size_type) {
    free(p);
  }

  template<typename tp1> struct rebind {
    typedef allocator<tp1> other;
  };
};
}
