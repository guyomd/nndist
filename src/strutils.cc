#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>

std::string trim(const std::string& str,
                 const std::string& whitespace = " \t")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}


std::vector<std::string> splitString(const std::string& text, char delimiter) {
    std::vector<std::string> tokens;
    std::istringstream iss(text);
    std::string token;
    while (std::getline(iss, token, delimiter)) {
        token = trim(token);
        if (!token.empty()) {
            tokens.push_back(token);
        }
    }
    return tokens;
}


std::string formatAlpha0(double alpha0) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << alpha0;
    std::string text = oss.str();
    if (text.find('.') != std::string::npos) {
        while (!text.empty() && text.back() == '0') {
            text.pop_back();
        }
        if (!text.empty() && text.back() == '.') {
            text.pop_back();
        }
    }
    if (text.empty()) {
        text = "0";
    }
    std::replace(text.begin(), text.end(), '.', 'p');
    return text;
}
