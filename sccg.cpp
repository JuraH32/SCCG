
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <chrono>
#include <thread>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>
#include <cctype>
#include <cstring>
#include <sys/stat.h>

class SCCGD {
public:
    class Position {
    public:
        int startinRef;
        int endinRef;
        
        Position() : startinRef(0), endinRef(0) {}
    };

    static std::vector<Position> L_list;
    static std::vector<Position> N_list;
    static std::string meta_data;
    static int length;
    static int line_length;

    static std::string readSeq(const std::string& sequenceFileName) {
        std::ifstream file(sequenceFileName);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + sequenceFileName);
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
            throw std::runtime_error("Cannot open file: " + sequenceFileName);
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
        } catch (const std::exception& e) {
            std::cerr << "Exception: " << e.what() << std::endl;
        }
    }

    static std::string reconstruct(const std::string& inFileName, const std::string& reference) {
        std::ifstream file(inFileName);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + inFileName);
        }
        
        std::string line;
        std::stringstream stringbuilder;
        
        int prev_end = 0, index = 0;
        
        while (std::getline(file, line)) {
            if (index == 4) {
                if (line.find(",") != std::string::npos) {
                    size_t comma_pos = line.find(",");
                    std::string begin_str = line.substr(0, comma_pos);
                    std::string end_str = line.substr(comma_pos + 1);
                    
                    int begin = std::stoi(begin_str);
                    int end = std::stoi(end_str);
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
        return stringbuilder.str();
    }

    static void write(const std::string& fileName, const std::string& text) {
        std::ofstream file(fileName);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file for writing: " + fileName);
        }
        
        file << text;
        file.flush();
        file.close();
    }

    static long getCpuTime() {
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = now.time_since_epoch();
        return std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();
    }

    static void test(const std::string& ref, const std::string& tar) {
        for (size_t i = 0; i < ref.length(); i++) {
            if (ref[i] != tar[i]) {
                std::cout << "error:" << i << std::endl;
            }
        }
    }

    static std::vector<std::string> split(const std::string& str, const std::string& delimiter) {
        std::vector<std::string> tokens;
        size_t start = 0;
        size_t end = str.find(delimiter);
        
        while (end != std::string::npos) {
            if (end != start) {
                tokens.push_back(str.substr(start, end - start));
            }
            start = end + delimiter.length();
            end = str.find(delimiter, start);
        }
        
        if (start < str.length()) {
            tokens.push_back(str.substr(start));
        }
        
        return tokens;
    }

    static void deleteFile(const std::string& filename) {
        struct stat buffer;
        if (stat(filename.c_str(), &buffer) == 0) {
            std::remove(filename.c_str());
        }
    }

    static int main(int argc, char* argv[]) {
        try {
            if (argc != 4) { 
                std::cout << "Make sure you have inputted 3 arguments." << std::endl;
                std::exit(0);
            }
            
            std::string reference_file = argv[1];
            std::string in_file_7z = argv[2];
            
            std::vector<std::string> chrs = split(in_file_7z, "/");
            std::string final_file = std::string(argv[3]) + "/" + chrs[chrs.size()-1];
            
            size_t pos = final_file.find(".7z");
            if (pos != std::string::npos) {
                final_file = final_file.substr(0, pos) + final_file.substr(pos + 3);
            }
            
            auto startDate = std::chrono::system_clock::now();
            long startCpuTimeNano = getCpuTime();
            
            std::time_t start_time = std::chrono::system_clock::to_time_t(startDate);
            std::cout << "Start time: " << std::ctime(&start_time);
            
            std::cout << argv[2] << " is decompressing..." << std::endl;
            
            deleteFile(final_file);
            
            use7zip(in_file_7z, std::string(argv[3]) + "/");
            
            std::ifstream file(final_file);
            if (!file.is_open()) {
                throw std::runtime_error("Cannot open file: " + final_file);
            }
            
            std::string line;
            int Pindex = 0;
            
            L_list.clear();
            N_list.clear();
            
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
                    std::vector<std::string> strings = split(line, " ");
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
                    std::vector<std::string> strings = split(line, " ");
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
            
            std::string target_string = reconstruct(final_file, reference);
            
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
            
            std::stringstream target_builder;
            target_builder << final_string[0];
            
            for (int t = 1; t < static_cast<int>(final_string.length()); t++) {
                if (t % line_length == 0) {
                    target_builder << "\n";
                }
                target_builder << final_string[t];
            }
            
            final_string = meta_data + target_builder.str() + "\n";
            
            
            write(final_file, final_string);
            
            
            long taskCpuTimeNano = getCpuTime() - startCpuTimeNano;
            std::cout << "Decompressed time: " << static_cast<double>(taskCpuTimeNano) / 1000000000.0 << " seconds." << std::endl;
            std::cout << "Decompressed file is " << argv[3] << "/" << chrs[chrs.size()-1] << std::endl;
            std::cout << "Done\n" << "----------------------------------------------------------------------------" << std::endl;
            
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
        
        return 0;
    }
};

std::vector<SCCGD::Position> SCCGD::L_list;
std::vector<SCCGD::Position> SCCGD::N_list;
std::string SCCGD::meta_data;
int SCCGD::length;
int SCCGD::line_length;

int main(int argc, char* argv[]) {
    return SCCGD::main(argc, argv);
}

