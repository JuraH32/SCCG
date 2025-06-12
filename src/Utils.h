#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

namespace Utils {
#ifdef __APPLE__
    const std::string sevenZipCommand = "7zz";
    const std::string compressionCommand = "-mm=ppmd";
#else
    const std::string sevenZipCommand = "./7za";
    const std::string compressionCommand = "m0=PPMD";
#endif

    inline int compressWith7zip(const std::string& inputFile) {
        std::string command = sevenZipCommand + " a " + inputFile + ".7z " + inputFile + " " + compressionCommand;
        return system(command.c_str());
    }

    inline int decompressWith7zip(const std::string& archivePath, const std::string& unzipDir) {
        std::string cmd = sevenZipCommand + " e \"" + archivePath + "\" -o\"" + unzipDir + "\" -aos";
        return system(cmd.c_str());
    }

    inline std::vector<std::string> splitPath(const std::string& path) {
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

    inline void createFileIfNotExists(const std::string& filename) {
        std::ofstream file(filename, std::ios::app);
        if (!file) {
            std::cerr << "Error: Could not create file: " << filename << std::endl;
        } else {
            std::cout << "Created file: " << filename << std::endl;
        }
        file.close();
    }

    inline void removeFileIfExists(const std::string& filename, bool create_if_not_exists = false) {
        if (std::remove(filename.c_str()) != 0) {
            std::cerr << "Warning: Could not remove file: " << filename << std::endl;
        } else {
            std::cout << "Removed file: " << filename << std::endl;
        }
        if (create_if_not_exists) {
            createFileIfNotExists(filename);
        }
    }
}