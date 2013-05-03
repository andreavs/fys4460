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


//// array::fill example
//#include <iostream>
//#include <fstream>
//using namespace std;

//int main (int args, char *argv[]) {
//    int N = atoi(argv[1]);

//    double A[N];
//    for(int i=0;i<N;i++) {
//        A[i] = i+i/(float)N;
//        cout << A[i] << endl;
//    }

//    // char buffer[100];
//    ofstream myFile ("data.bin", ios::out | ios::binary);
//    myFile.write (reinterpret_cast<char*>(&N), sizeof(int)); // First write how many doubles we will write
//    myFile.write (reinterpret_cast<char*>(A), N*sizeof(double));
//    myFile.close();

//    return 0;
//}

//// array::fill example
//#include <iostream>
//#include <fstream>
//using namespace std;

//int main () {
//	int N;

//	ifstream myFile ("data.bin", ios::in | ios::binary);
//	myFile.read (reinterpret_cast<char*>(&N), sizeof(int)); // Read the number of doubles
//	double *A = new double[N];
//    myFile.read (reinterpret_cast<char*>(A), N*sizeof(double));

//    for(int i=0;i<N;i++) {
//    	cout << A[i] << endl;
//    }

//	return 0;
//}
