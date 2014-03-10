#ifndef OPERATOR_ENUMS_H_
#define OPERATOR_ENUMS_H_

#include <stdexcept>

typedef enum {ADD,SUB,MUL,DIV,NEG, SQRT, SIN, LOG, EXP, POW, TAN, COS} arithmetic_operator;
typedef enum {AND,OR,NOT} boolean_operator;
typedef enum {LT,LE,GT,GE,EQ} comparison_operator;

inline bool is_strict(comparison_operator s) {return s==LT || s==GT;};
inline bool is_EQ(comparison_operator s) {return s==EQ;};
inline bool is_LT_or_LE(comparison_operator s) {return (s==LE || s==LT);};
inline bool is_LT_or_LE_or_EQ(comparison_operator s) {return (s==LE || s==LT || s==EQ);};
inline bool is_LE_or_EQ(comparison_operator s) {return (s==LE || s==EQ);};
inline bool is_GT_or_GE(comparison_operator s) {return (s==GE || s==GT);};
inline bool is_GT_or_GE_or_EQ(comparison_operator s) {return (s==GE || s==GT || s==EQ);};
inline bool is_GE_or_EQ(comparison_operator s) {return (s==GE || s==EQ);};

inline comparison_operator get_weaker_sign(comparison_operator s1,comparison_operator s2) {
	if (s1==LT && is_LT_or_LE_or_EQ(s2)) return LT;
	if (s2==LT && is_LT_or_LE_or_EQ(s1)) return LT;
	if (s1==LE && is_LE_or_EQ(s2)) return LE;
	if (s2==LE && is_LE_or_EQ(s1)) return LE;
	if (s1==GT && is_GT_or_GE_or_EQ(s2)) return GT;
	if (s2==GT && is_GT_or_GE_or_EQ(s1)) return GT;
	if (s1==GE && is_GE_or_EQ(s2)) return GE;
	if (s2==GE && is_GE_or_EQ(s1)) return GE;
	if (s1==EQ && s2==EQ) return EQ;
	// otherwise signs are opposite
	throw std::runtime_error("can't combine opposite signs");
	return EQ; // arbitrary return value so the compiler doesn't complain
};

inline comparison_operator get_flipped_sign(comparison_operator s) {
	if (s==LT) return GT;
	if (s==LE) return GE;
	if (s==GT) return LT;
	if (s==GE) return LE;
	// if (s==EQ) 
	return EQ;
};

inline comparison_operator get_closed_sign(comparison_operator s) {
	if (s==LT) return LE;
	if (s==GT) return GE;
	return s;
};

inline std::string sign_string(comparison_operator s) {
	switch (s) {
	case LT:
		return "<";
	case LE:
		return "<=";
	case EQ:
		return "==";
	case GE:
		return ">=";
	case GT:
		return ">";
	default:
		throw std::runtime_error("sign_string: unknown sign");
		return "";
	}
}
;

/** Returns true if the operator is commutative */
inline bool is_commutative(arithmetic_operator op) {
	return (op==ADD || op == MUL);
};

/** Returns true if the operators are commutative */
inline bool is_commutative(arithmetic_operator op1,arithmetic_operator op2) {
	return is_commutative(op1) && (op1==op2);
};

#endif /*OPERATOR_ENUMS_H_*/
