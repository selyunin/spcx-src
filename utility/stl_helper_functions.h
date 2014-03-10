#ifndef STL_HELPER_FUNCTIONS_H_
#define STL_HELPER_FUNCTIONS_H_

#include <string>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm> /* for std::includes -wsc */

#include "infix_ostream_iterator.h"

// -----------------------------------------------------------------------------------------
// String Functions
// -----------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------------
// Convert to string
// -----------------------------------------------------------------------------------------

/** Convert an object t of class T to string using stream output.
 *
 * Adapted from http://www.codeguru.com/forum/showthread.php?t=231056, by Gabriel Fleseriu
 */
template<class T> std::string to_string(const T& t) {
	std::ostringstream oss;
	if (!(oss << t))
		throw std::runtime_error("Could not convert object to stream.");
	return oss.str();
}

/** Convert an object t of class T to string using stream output, using a format specifier f.
 * Example: std::cout<<to_string<long>(123456, std::hex)<<std::endl;
 *
 * Adapted from http://www.codeguru.com/forum/showthread.php?t=231056, by Gabriel Fleseriu
 */
template<class T> std::string to_string(const T& t, std::ios_base & (*f)(std::ios_base&)) {
	std::ostringstream oss;
	if (!(oss << f << t))
		throw std::runtime_error("Could not convert object to stream.");
	return oss.str();
}

/** Convert integer i to string, but only up to 256 digits! */
std::string int2string(const int i);

/** Convert double d to string, but only up to 256 digits! */
std::string double2string(const double d);

// -----------------------------------------------------------------------------------------
// Convert from string
// -----------------------------------------------------------------------------------------

/** Convert a string to an object t of class T using stream input.
 *
 * Adapted from http://www.parashift.com/c++-faq-lite/misc-technical-issues.html#faq-39.2
 */
template<class T> T from_string(const std::string& s) {
	std::istringstream iss(s);
	T x;
	if (!(iss >> x))
		throw std::runtime_error("Could not convert string \"" + s + "\" to desired type.");
	return x;
}

/** Convert a string to a bool.
 *
 * Adapted from http://www.parashift.com/c++-faq-lite/misc-technical-issues.html#faq-39.2
 */
template<>
inline bool from_string<bool> (const std::string& s) {
	std::istringstream iss(s);
	bool x;
	if (!(iss >> std::boolalpha >> x))
		throw std::runtime_error("Could not convert string " + s + " to object.");
	return x;
}

// -----------------------------------------------------------------------------------------
// Split strings
// -----------------------------------------------------------------------------------------

/** Return all chars before encountering c.
 * @Note: By definition, str is returned if c is empty or c is not found.
 *       This is so that string_before + c + string_after == str.
 */
std::string string_before(const std::string& str, const std::string& c);

/** Return all chars after encountering c.
 * @Note: By definition, nothing is returned if c is empty or c is not found.
 *       This is so that string_before + c + string_after == str
 */
std::string string_after(const std::string& str, const std::string& c);

/** Return the substrings of str obtained by splitting the string
 * at every occurrence of the character delim.
 */
std::vector<std::string> split_string(const std::string& str, char delim);

/** Replace in str any occurrence of search_str by repl_str.
 *
 * taken from http://snipplr.com/view/1055/find-and-replace-one-string-with-another/
 * */
void replace(std::string& str, std::string search_str, std::string repl_str);

/** Replace in str the first occurrence of search_str by repl_str.
 * */
void replace_first(std::string& str, std::string search_str, std::string repl_str);

/** Strip the path from a filename string. */
std::string strip_path(const std::string& str);

/** Get the filename extension */
std::string get_file_extension(const std::string& fname);

/** Replaces special caracters in a stream by their XML equivqlent.
 *
 * Performs the following replacements:
 * &amp;	&
 *	&lt;	<
 *	&gt;	>
 *  &quot;	"
 *	&apos;	'
 *
 * */
std::string string_to_xml(const std::string& str);

/** Remove quotes from string
 *
 * Removes one set of double quotes from beginning and end of s.
 * A different pair of characters can be specified as cbeg and cend. Removal only
 * takes place if both beginning and end match.
 * Replacement is not recursive.
 * A single character is not changed.
 */
void trim_quotes(std::string& s, char cbeg='\"', char cend='\"');

/** Replace any characters in txt between bmark and emark (including) with repl. */
void delimited_replace(std::string& txt, const std::string& bmark,
		const std::string& emark, const std::string repl);

// -----------------------------------------------------------------------------------------
// Compare strings
// -----------------------------------------------------------------------------------------

/** Return true if the string str matches the string wild, where wild may contain
 * the wildcard characters ? (matches any one character) and $ (matches any number of characters).
 *
 * Taken from http://www.codeproject.com/string/wildcmp.asp
 * which is adapted from code by Jack Handy - jakkhandy@hotmail.com
 */
bool wildcmp(const std::string& wild, const std::string& str);

// -----------------------------------------------------------------------------------------
// Stream Functions
// -----------------------------------------------------------------------------------------

/** A redirector class for redirecting streams
 *
 * Stream source will be redirected to stream target as long as the redirector is alive.
 * The original state will be restored upon destruction.
 *
 * @Note This functionality is wrapped in a class to ensure that streams are properly
 * restored in case of an exception.
 * Recall that C++ encourages RAII (ressource acquisition is initialization) rather than
 * things like try/catch/finally. */
class stream_redirector {
public:
	/** The constructor initiates redirection. */
	stream_redirector(std::ios& source, std::ios& target) : my_source(source) {
		orig_stream = my_source.rdbuf(target.rdbuf());
	}
	;
	/** The destructor restores the original state. */
	virtual ~stream_redirector() {
		my_source.rdbuf(orig_stream);
	}
private:
	std::ios& my_source;
	std::streambuf * orig_stream;
};


/** Prints an STL vector in the format:
 * [0:element0,1:element1,...,n-1:element(n-1)]
 */
template<typename T> std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
	os << "[";
	for (unsigned int i = 0; i < v.size(); ++i) {
		if (i > 0)
			os << ",";
		os << i << ":" << v[i];
	}
	os << "]";
	return os;
}

/** Prints an STL set in the format:
 * {element1,element2,...,elementn}
 */
template<typename T> std::ostream& operator<<(std::ostream& os, const std::set<T>& v) {
	os << "{";
	for (typename std::set<T>::const_iterator it = v.begin(); it != v.end(); ++it) {
		if (it != v.begin())
			os << ",";
		os << (T) (*it);
	}
	os << "}";
	return os;
}

/** Prints an STL list in the format:
 * <element1,element2,...,elementn>
 */
template<typename T> std::ostream& operator<<(std::ostream& os, const std::list<T>& v) {
	os << "<";
	for (typename std::list<T>::const_iterator it = v.begin(); it != v.end(); ++it) {
		if (it != v.begin())
			os << ",";
		os << (T) (*it);
	}
	os << ">";
	return os;
}

/** Outputs an STL vector in the format:
 * [key1->value1,key2->value2,...,keyn->valuen]
 */
template<typename T, typename T2, typename T3, typename T4> std::ostream& operator<<(
		std::ostream& os, const std::map<T, T2, T3, T4>& v) {
	os << "[";
	for (typename std::map<T, T2, T3, T4>::const_iterator it = v.begin(); it != v.end(); ++it) {
		if (it != v.begin())
			os << ",";
		os << it->first << "->" << it->second;
	}
	os << "]";
	return os;
}

// -----------------------------------------------------------------------------------------
// Set Functions
// -----------------------------------------------------------------------------------------

/** Returns true iff s1 and s2 have no common element. */
template<typename T> bool set_is_disjoint(const std::set<T>& s1,
		const std::set<T>& s2) {
	std::set<T> out_set;
	std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
			std::inserter(out_set, out_set.begin()));
	return out_set.empty();
}
;

/** Return true iff s1 contains s2. */
template<typename T> bool set_contains(const std::set<T>& s1,
		const std::set<T>& s2) {
	return std::includes(s1.begin(), s1.end(), s2.begin(), s2.end());
}
;

/** Removes from s1 all the elements of s2. */
template<typename T> void set_difference_assign(std::set<T>& s1,
		const std::set<T>& s2) {
	for (typename std::set<T>::const_iterator it = s2.begin(); it != s2.end(); ++it) {
		s1.erase(*it);
	}
}
;

/** Adds to s1 all the elements of s2. */
template<typename T> void set_union_assign(std::set<T>& s1,
		const std::set<T>& s2) {
	s1.insert(s2.begin(), s2.end());
}
;

/** Constructs a set from a list, removing doubles. */
template<typename T> std::set<T> list_to_set(const std::list<T>& l) {
	std::set<T> s;
	for (typename std::list<T>::const_iterator it = l.begin(); it != l.end(); ++it) {
		s.insert(*it);
	}
	return s;
}
;

#endif /*STL_HELPER_FUNCTIONS_H_*/
