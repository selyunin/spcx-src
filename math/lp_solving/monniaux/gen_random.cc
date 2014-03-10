#include <iostream>
#include <gmpxx.h>
#include <stdlib.h>

mpz_class random_z(gmp_randstate_t state, mpz_class max) {
  mpz_class r;
  mpz_urandomm(r.get_mpz_t(), state, mpz_class(2*max+1).get_mpz_t());
  r -= max;
  return r;
}

unsigned random_ui(gmp_randstate_t state, unsigned n) {
  return gmp_urandomm_ui(state, n);
}

void printNumber(std::ostream& out, mpz_class nr) {
  if (nr < 0) {
    out << "(~ " << -nr << ")";
  } else {
    out << nr;
  }
}

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "usage: gen_random vars conjuncts max [seed]" << std::endl;
    return 1;
  }

  unsigned vars = atoi(argv[1]), conjuncts=atoi(argv[2]);
  mpz_class max(argv[3]);

  gmp_randstate_t randgen;
  gmp_randinit_default(randgen);
  gmp_randseed(randgen, mpz_class(argv[4]).get_mpz_t());

  std::ostream& out = std::cout;

  const char* operators[] = {"<", ">", "<=", ">=" };

  out << "(benchmark lemma" << std::endl;
  out << ":logic QF_LRA" << std::endl;

  for(unsigned var=0; var<vars; var++) {
    out << ":extrafuns ((v" << var << " Real))" << std::endl;
  }

  out << ":formula (and" << std::endl;

  for(unsigned conjunct=0; conjunct<conjuncts; conjunct++) {
    out << "(" << operators[random_ui(randgen, 4)] << " (+ ";

    for(unsigned var=0; var<vars; var++) {
      out << "(* ";
      printNumber(out, random_z(randgen, max));
      out << " v" << var << ")";
    }
    out << ") ";
    printNumber(out, random_z(randgen, max));
    out << ")" << std::endl;
  }

  out << "))" << std::endl;
}
