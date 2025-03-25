#include "Functions.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <sys/resource.h>


// formats time from seconds into mm:ss
string format_time(double seconds) {

    // get minutes and seconds
    int minutes = static_cast<int>(seconds) / 60;
    int sec = static_cast<int>(seconds) % 60;
    
    // put them into string
    stringstream ss;
    ss << setfill('0') << setw(2) << minutes << ":" << setfill('0') << setw(2) << sec;
    return ss.str();
}


// MEMORY: function to print out maximum memory usage
long long get_maxrss_memory() {
    
    // declare a rusage structure to store resource usage information
    struct rusage usage;

    // get resource usage statistics for the current process
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
    
        // print the RSS memory size in megabytes
        long long rssmax = usage.ru_maxrss;
        cout << "max RSS memory size: " << rssmax/1024.0/1024.0 << " MB" << endl;
        return rssmax;

    } else {

        cerr << "Error getting resource usage." << endl;

        return 0;
    }
}