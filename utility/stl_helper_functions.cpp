#include "utility/stl_helper_functions.h"
#include <iostream>
#include <cstdio>

using namespace std;

string int2string(const int i)
{
  char chr[256];
  string str;
  sprintf(chr,"%i",i);
  str=chr;
  return str;
}

string double2string(const double d)
{
  char chr[256];
  string str;
  sprintf(chr,"%g",d);
  str=chr;
  return str;
}

string string_before(const string& str, const string& c)
{
	string out("");
	bool ret_all=false;
	
	if (c.size()>0)
	{
		string::size_type i=str.find(c);
		if (i!=string::npos)
		{
			if (i>0)
				out=str.substr(0,i);
		}
		else
			ret_all=true;
	}
	else
		ret_all=true;
	
	if (ret_all) 
		return str;
	else
		return out;
}

string string_after(const string& str, const string& c)
{
	string out("");
	if (c.size()>0)
	{
		string::size_type i=str.find(c);
		if (i!=string::npos && i+1<str.size())
			out=str.substr(i+c.size(),str.size()-i-c.size());
	};
						
	return out;
}


bool wildcmp(const string& wild, const string& str) 
// taken from http://www.codeproject.com/string/wildcmp.asp
// which is adapted from code by Jack Handy - jakkhandy@hotmail.com
{
	unsigned int cp=0, mp=0;

	unsigned int i=0;
	unsigned int j=0;
	while (i < str.length() && j < wild.length() && wild[j] != '$')
	{
		if ((wild[j] != str[i]) && (wild[j] != '?')) 
		{
			return false;
		}
		i++;
		j++;
	}
	
	while (i<str.length()) 
	{
		if (j<wild.length() && wild[j] == '$') 
		{
			if ((j++)>=wild.length()) 
			{
				return true;
			}
			mp = j;
			cp = i+1;
		} 
		else if (j<wild.length() && (wild[j] == str[i] || wild[j] == '?')) 
		{
			j++;
			i++;
		} 
		else 
		{
			j = mp;
			i = cp++;
		}
	}
	
	while (j<wild.length() && wild[j] == '$')
	{
		j++;
	}
	return j>=wild.length();
}

vector<string> split_string(const string& str, char delim)
{
	istringstream instr;
	instr.str(str);
	string buf;
	vector<string> vstr;
	while (getline(instr,buf,delim))
		vstr.push_back(buf);
	// if the last char is the delimiter, add an empty string
	if (str.size()>1 && str[str.size()-1]==delim)
		vstr.push_back("");
	return vstr;
}

void replace(std::string& tInput, std::string tFind, std::string tReplace) {
	// taken from http://snipplr.com/view/1055/find-and-replace-one-string-with-another/

	size_t uPos = 0;
	size_t uFindLen = tFind.length();
	size_t uReplaceLen = tReplace.length();

	if (uFindLen != 0) {
		for (; (uPos = tInput.find(tFind, uPos)) != std::string::npos;) {
			tInput.replace(uPos, uFindLen, tReplace);
			uPos += uReplaceLen;
		}
	}
}

void replace_first(std::string& tInput, std::string tFind, std::string tReplace) {
	size_t uPos = 0;
	size_t uFindLen = tFind.length();
	size_t uReplaceLen = tReplace.length();

	if (uFindLen != 0) {
		if ((uPos = tInput.find(tFind, uPos)) != std::string::npos) {
			tInput.replace(uPos, uFindLen, tReplace);
			uPos += uReplaceLen;
		}
	}
}

std::string string_to_xml(const std::string& str) {
	std::string res(str);
	replace(res, "&", "&amp;");
	replace(res, "<", "&lt;");
	replace(res, ">", "&gt;");
	replace(res, "\"", "&quot;");
	replace(res, "'", "&apos;");
	return res;
}

std::string strip_path(const std::string& str) {
	string out("");

	string::size_type i=str.find("/");
	string::size_type ilast=0;
	while (i!=string::npos)
	{
		ilast=i+1;
		i=str.find("/",i+1);
	}
	out=str.substr(ilast,str.length()-ilast);

	return out;
}


void trim_quotes(std::string& s, char cbeg, char cend) {
	if (s.size() > 1) {
		unsigned int npos = s.size() - 1;
		if (s[0] == cbeg && s[npos] == cend) {
			// erase leading quote
			s.erase(0, 1);
			// erase ending quote
			if (npos>0)
				s.erase(npos-1, 1);
		}
	}
}

void delimited_replace(std::string& txt, const std::string& bmark,
		const std::string& emark, const std::string repl) {
	//			const std::tr1::regex pattern("\\[[^\\[\\]]*\\]");
	//			std::string replace = "";
	//			text = std::tr1::regex_replace(text, pattern, replace);
	std::size_t epos = 0;
	std::size_t bpos = txt.find(bmark, epos);
	while (bpos != std::string::npos) {
		epos = txt.find(emark, bpos);
		if (epos != std::string::npos) {
			// Replace everything from bpos to epos+emark.length()-1.
			// This has length (epos+emark.length()-1-bpos)+1
			//    = epos-bpos+emark.length()
			txt.replace(bpos, epos - bpos + emark.length(), repl);

			// Look for the next occurene, starting with the next
			// character, which is at bpos+repl.length()
			bpos = txt.find(bmark, bpos + repl.length());
		} else
			bpos = std::string::npos;
	}
}

std::string get_file_extension(const std::string& fname)
{
    if(fname.find_last_of(".") != std::string::npos)
        return fname.substr(fname.find_last_of(".")+1);
    return "";
}
