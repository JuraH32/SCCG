#include <cstdlib>
#include <string>

namespace Utils {
    // Basic example - add error checking and proper path handling!
    int compressWith7zip(const std::string& inputFile, const std::string& outputFile) {
        std::string command = "./7za a " + outputFile + " " + inputFile;
        return system(command.c_str());
    }

    int decompressWith7zip(const std::string& inputFile, const std::string& outputDir) {
        std::string command = "./7za x " + inputFile + " -o" + outputDir;
        return system(command.c_str());
    }
}