#include <string>
#include <sstream>
#include <vector>


//http://stackoverflow.com/questions/236129/splitting-a-string-in-c by Evan Teran
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

bool to_bool(std::string const& s) {
     return s != "0";
}
