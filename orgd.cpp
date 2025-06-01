
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cctype>
#include <chrono>
#include <thread>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <sys/stat.h>

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/times.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>
#endif

class ORGD {
public:
    class Position {
    public:
        int startinRef;
        int endinRef;
    };

    static std::vector<Position> L_list;
    static std::vector<Position> N_list;
    static std::string meta_data;
    static int length;
    static int line_length;

    static std::string readSeq(const std::string& sequenceFileName) {
        std::ifstream file(sequenceFileName);
        if (!file.is_open()) {
            throw std::runtime_error("IOException: Cannot open file");
        }
        
        std::string line;
        std::stringstream stringbuilder;
        
        std::getline(file, line);
        
        while (std::getline(file, line)) {
            stringbuilder << line;
        }
        
        file.close();
        return stringbuilder.str();
    }

    static std::string readrefSeq(const std::string& sequenceFileName) {
        std::ifstream file(sequenceFileName);
        if (!file.is_open()) {
            throw std::runtime_error("IOException: Cannot open file");
        }
        
        std::string line;
        int line_length;
        char temp_ch;
        
        std::stringstream stringbuilder;
        
        std::getline(file, line);
        
        while (std::getline(file, line)) {
            line_length = line.length();
            for (int i = 0; i < line_length; i++) {
                temp_ch = line[i];
                if (!std::isupper(temp_ch)) {
                    temp_ch = std::toupper(temp_ch);
                }
                if (temp_ch != 'N') {
                    stringbuilder << temp_ch;
                }
            }
        }
        
        file.close();
        return stringbuilder.str();
    }

    static void use7zip(const std::string& filename, const std::string& Dfilename) {
        struct stat buffer;
        if (stat(filename.c_str(), &buffer) != 0) {
            return;
        }
        
        std::string exec = "./7za e " + filename + " -o" + Dfilename + " -aos";
        try {
            int result = std::system(exec.c_str());
            (void)result; 
        } catch (const std::exception& e) {
            std::cerr << "Exception: " << e.what() << std::endl;
        }
    }

    static std::stringstream reconstruct(const std::string& inFileName, const std::string& reference) {
        std::ifstream file(inFileName);
        if (!file.is_open()) {
            throw std::runtime_error("IOException: Cannot open file");
        }
        
        std::string line;
        std::stringstream stringbuilder;
        
        int prev_end = 0, index = 0;
        
        while (std::getline(file, line)) {
            if (index == 4) {
                if (line.find(",") != std::string::npos) {
                    size_t comma_pos = line.find(",");
                    std::string first_part = line.substr(0, comma_pos);
                    std::string second_part = line.substr(comma_pos + 1);
                    
                    int begin = std::stoi(first_part);
                    int end = std::stoi(second_part);
                    begin += prev_end;
                    end += begin;
                    
                    try {
                        std::string text = reference.substr(begin, end - begin + 1);
                        stringbuilder << text;
                        prev_end = end;
                    } catch (const std::exception& x) {
                        std::cout << "come" << prev_end << begin << end << std::endl;
                    }
                } else if (line.length() > 0) {
                    stringbuilder << line;
                }
                continue;
            } else if (index < 4) {
                index++;
                continue;
            }
        }
        
        file.close();
        return stringbuilder;
    }

    static void write(const std::string& fileName, const std::string& text) {
        std::ofstream file(fileName);
        if (!file.is_open()) {
            throw std::runtime_error("IOException: Cannot create file");
        }
        
        file << text;
        file.flush();
        file.close();
    }

    static long long getCpuTime() {
#ifdef _WIN32
        FILETIME creationTime, exitTime, kernelTime, userTime;
        if (GetThreadTimes(GetCurrentThread(), &creationTime, &exitTime, &kernelTime, &userTime)) {
            ULARGE_INTEGER userTimeInt, kernelTimeInt;
            userTimeInt.LowPart = userTime.dwLowDateTime;
            userTimeInt.HighPart = userTime.dwHighDateTime;
            kernelTimeInt.LowPart = kernelTime.dwLowDateTime;
            kernelTimeInt.HighPart = kernelTime.dwHighDateTime;
            // Convert from 100-nanosecond units to nanoseconds
            return (userTimeInt.QuadPart + kernelTimeInt.QuadPart) * 100;
        }
        return 0;
#else
        struct rusage usage;
        if (getrusage(RUSAGE_THREAD, &usage) == 0) {
            // Convert to nanoseconds
            return (usage.ru_utime.tv_sec + usage.ru_stime.tv_sec) * 1000000000LL +
                   (usage.ru_utime.tv_usec + usage.ru_stime.tv_usec) * 1000LL;
        }
        return 0;
#endif
    }

    static void test(const std::string& ref, const std::string& tar) {
        for (size_t i = 0; i < ref.length(); i++) {
            if (ref[i] != tar[i]) {
                std::cout << "error:" << i << std::endl;
            }
        }
    }

    static int main(int argc, char* argv[]) {
        try {
            if (argc != 4) {
                std::cout << "Make sure you have inputted 3 arguments." << std::endl;
                std::exit(0);
            }
            
            std::string final_file;
            std::string reference_file;
            std::string unzip_folder = std::string(argv[3]) + "/unzip";
            std::string output;
            
            std::string chrName[] = {
                "chr1.fa", "chr2.fa", "chr3.fa", "chr4.fa", "chr5.fa", "chr6.fa", "chr7.fa", "chr8.fa", "chr9.fa",
                "chr10.fa", "chr11.fa", "chr12.fa", "chr13.fa", "chr14.fa", "chr15.fa", "chr16.fa", "chr17.fa", "chr18.fa", "chr19.fa", "chr20.fa",
                "chr21.fa", "chr22.fa", "chrX.fa", "chrY.fa"
            };
            
            use7zip(argv[2], unzip_folder);
            
            auto startDate = std::chrono::system_clock::now();
            long long startCpuTimeNano = getCpuTime();
            std::time_t start_time = std::chrono::system_clock::to_time_t(startDate);
            std::cout << "Start time: " << std::ctime(&start_time);
            std::cout << argv[2] << " is decompressing..." << std::endl;
            
            for (int i = 0; i < 24; i++) {
                reference_file = std::string(argv[1]) + "/" + chrName[i];
                final_file = unzip_folder + "/" + chrName[i];
                output = std::string(argv[3]) + "/" + chrName[i];
                
                std::remove(output.c_str());
                
                std::ifstream file(final_file);
                if (!file.is_open()) {
                    throw std::runtime_error("IOException: Cannot open file");
                }
                
                std::string line;
                int Pindex = 0;
                
                L_list = std::vector<Position>();
                N_list = std::vector<Position>();
                
                while (std::getline(file, line)) {
                    if (Pindex == 1) {
                        line_length = std::stoi(line);
                        Pindex++;
                        continue;
                    } else if (Pindex == 2) {
                        if (line.empty() || line.length() <= 0) {
                            Pindex++;
                            continue;
                        }
                        std::istringstream iss(line);
                        std::string token;
                        std::vector<std::string> strings;
                        while (iss >> token) {
                            strings.push_back(token);
                        }
                        
                        int j = 0, prev = 0;
                        while (j < static_cast<int>(strings.size())) {
                            Position position;
                            position.startinRef = prev + std::stoi(strings[j]);
                            j++;
                            position.endinRef = position.startinRef + std::stoi(strings[j]);
                            j++;
                            L_list.push_back(position);
                            prev = position.endinRef;
                        }
                        Pindex++;
                        continue;
                    } else if (Pindex == 3) {
                        if (line.empty() || line.length() <= 0) {
                            Pindex++;
                            continue;
                        }
                        std::istringstream iss(line);
                        std::string token;
                        std::vector<std::string> strings;
                        while (iss >> token) {
                            strings.push_back(token);
                        }
                        
                        int j = 0, prev = 0;
                        while (j < static_cast<int>(strings.size())) {
                            Position position;
                            position.startinRef = prev + std::stoi(strings[j]);
                            j++;
                            position.endinRef = position.startinRef + std::stoi(strings[j]);
                            j++;
                            N_list.push_back(position);
                            prev = position.endinRef;
                        }
                        Pindex++;
                        continue;
                    } else if (Pindex == 0) {
                        meta_data = line + "\n";
                        Pindex++;
                        continue;
                    }
                }
                file.close();
                
                std::string reference = "";
                
                if (N_list.size() > 0) {
                    reference = readrefSeq(reference_file);
                } else if (N_list.size() <= 0) {
                    reference = readSeq(reference_file);
                    std::transform(reference.begin(), reference.end(), reference.begin(), ::toupper);
                }
                
                std::stringstream target_string_stream = reconstruct(final_file, reference);
                std::string target_string = target_string_stream.str();
                
                std::stringstream interim;
                int index = 0, iterator = 0;
                bool accept = false;
                for (const Position& position : N_list) {
                    while (true) {
                        if (iterator >= position.startinRef && iterator < position.endinRef) {
                            interim << 'N';
                        } else if (iterator < position.startinRef && !accept) {
                            interim << target_string[index];
                            index = index + 1;
                            if (index >= static_cast<int>(target_string.length())) {
                                accept = true;
                            }
                        }
                        iterator++;
                        
                        if (iterator >= position.endinRef) {
                            break;
                        }
                    }
                }
                if (N_list.size() > 0) {
                    while (index < static_cast<int>(target_string.length())) {
                        interim << target_string[index];
                        index++;
                        iterator++;
                    }
                }
                
                if (interim.str().length() > 0) {
                    target_string = interim.str();
                } else {
                }
                
                for (const Position& position : L_list) {
                    for (int j = position.startinRef; j <= position.endinRef; j++) {
                        if (j == static_cast<int>(target_string.length())) {
                            break;
                        }
                        target_string[j] = std::tolower(target_string[j]);
                    }
                }
                
                std::string final_string = target_string;
                
                std::stringstream target_string_builder;
                target_string_builder << final_string[0];
                
                for (int t = 1; t < static_cast<int>(final_string.length()); t++) {
                    if (t % line_length == 0) {
                        target_string_builder << "\n";
                    }
                    target_string_builder << final_string[t];
                }
                
                final_string = meta_data + target_string_builder.str() + "\n";
                
                write(output, final_string);
                
                std::cout << "Decompressed file is " << argv[3] << "/" << chrName[i] << std::endl;
            }
            
            long long taskCpuTimeNano = getCpuTime() - startCpuTimeNano;
            std::cout << "Decompressed time: " << static_cast<double>(taskCpuTimeNano) / 1000000000.0 << " seconds." << std::endl;
            std::cout << "Done" << std::endl << "----------------------------------------------------------------------------" << std::endl;
            std::system(("rm -rf " + std::string(argv[3]) + "/unzip").c_str());
            
        } catch (const std::exception& e) {
            std::cerr << "Exception: " << e.what() << std::endl;
            return 1;
        }
        
        return 0;
    }
};

std::vector<ORGD::Position> ORGD::L_list;
std::vector<ORGD::Position> ORGD::N_list;
std::string ORGD::meta_data;
int ORGD::length;
int ORGD::line_length;

int main(int argc, char* argv[]) {
    return ORGD::main(argc, argv);
}

