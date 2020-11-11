#pragma once
/*
Utilities for I/O
*/

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

namespace bats {
namespace util {
namespace io {

template<typename T>
T parse_argv(
	const int argc,
	char** argv,
	const std::string &&token,
	const T default_return
) {
	int i = 1; // skip first entry
	while (i < argc) {
		std::string entry(argv[i]);
		std::istringstream iss(entry);
		std::string prefix, postfix;
		getline(iss, prefix, '=');
		if (token.compare(prefix) == 0) {
			getline(iss, postfix);
			return T(std::stod(postfix));
		}
		i++;
	}
	return default_return;
}

template<>
std::string parse_argv(
	const int argc,
	char** argv,
	const std::string &&token,
	const std::string default_return
) {
	int i = 1; // skip first entry
	while (i < argc) {
		std::string entry(argv[i]);
		std::istringstream iss(entry);
		std::string prefix, postfix;
		getline(iss, prefix, '=');
		if (token.compare(prefix) == 0) {
			getline(iss, postfix);
			return postfix;
		}
		i++;
	}
	return default_return;
}


} // namespace io
} // namespace util
} // namespace bats
