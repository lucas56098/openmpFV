#include "Profiler.h"

unordered_map<string, chrono::high_resolution_clock::time_point> Profiler::m_StartTimes;
unordered_map<string, long long> Profiler::m_Timings;

void Profiler::StartTimer(const string& name) {
    m_StartTimes[name] = chrono::high_resolution_clock::now();
}

void Profiler::EndTimer(const string& name) {
    auto endTime = chrono::high_resolution_clock::now();
    auto startTime = m_StartTimes[name];
    m_Timings[name] += chrono::duration_cast<chrono::microseconds>(endTime - startTime).count();
}

void Profiler::PrintResults() {
    cout << "\n=== Profiling Results ===\n";
    long long totalRuntime = 0;
    long long parallelizedTime = 0;

    for (const auto& entry : m_Timings) {
        cout << "[PROFILE] " << entry.first << " took " << entry.second << " microseconds in total.\n";
        
        if (entry.first.find("(par)") != string::npos) {
            parallelizedTime += entry.second;
        }
        
        if (entry.first == "TOTAL_RUNTIME") {
            totalRuntime = entry.second;
        }
    }

    double parallelFraction = 0.0;
    if (totalRuntime > 0) {
        parallelFraction = static_cast<double>(parallelizedTime) / totalRuntime;
    }
    cout << "\nTOTAL_RUNTIME: " << totalRuntime << " microseconds\n";
    cout << "PARALLELIZED_TIME: " << parallelizedTime << " microseconds\n";
    cout << "PARALLEL_FRACTION: " << parallelFraction * 100.0 << " %\n";
    cout << "=========================\n";
}
