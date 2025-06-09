#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

namespace Utils {
    int compressWith7zip(const std::string& inputFile) {
        std::string command = "./7za a " + inputFile + ".7z " + inputFile + " m0=PPMD";
        return system(command.c_str());
    }

    int decompressWith7zip(const std::string& inputFile, const std::string& outputDir) {
        std::string command = "../7za x " + inputFile + " -o" + outputDir + " -aos";
        return system(command.c_str());
    }

    std::vector<std::string> splitPath(const std::string& path) {
        std::vector<std::string> parts;
        size_t pos = 0;
        size_t found;
        while ((found = path.find_first_of("/\\", pos)) != std::string::npos) {
            if (found > pos) {
                parts.push_back(path.substr(pos, found - pos));
            }
            pos = found + 1;
        }
        if (pos < path.length()) {
            parts.push_back(path.substr(pos));
        }
        return parts;
    }

    void createFileIfNotExists(const std::string& filename) {
        std::ofstream file(filename, std::ios::app);
        if (!file) {
            std::cerr << "Error: Could not create file: " << filename << std::endl;
        } else {
            std::cout << "Created file: " << filename << std::endl;
        }
        file.close();
    }

    void removeFileIfExists(const std::string& filename, bool create_if_not_exists = false) {
        if (std::remove(filename.c_str()) != 0) {
            std::cerr << "Warning: Could not remove file: " << filename << std::endl;
        } else {
            std::cout << "Removed file: " << filename << std::endl;
        }
        if (create_if_not_exists) {
            createFileIfNotExists(filename);
        }
    }

    void savePositionsToFile(const std::string& filename, const std::vector<size_t>& positions) {
        std::ofstream outFile(filename, std::ios::app);
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
            return;
        }
        outFile << "LOWERCASE_POSITIONS_COUNT " << positions.size() << "\n";
        for (size_t pos : positions) {
            outFile << pos << "\n";
        }
        outFile.close();
    }
}