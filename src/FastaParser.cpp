#include "FastaParser.h"

std::string FastaParser::readSequence(const std::string& filename) {
    // Attempt to open the specified file for reading.
    std::ifstream inFile(filename);
    if (!inFile) {
        // If the file cannot be opened, print an error message to stderr and return an empty string.
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return "";
    }

    std::string sequence;      // Initialize an empty string to store the concatenated sequence.
    std::string line;          // String to hold each line read from the file.

    while (std::getline(inFile, line)) {
        // Skip empty lines or header lines (which start with '>').
        if (line.empty() || line[0] == '>') {
            continue;
        }
        // Remove carriage return characters, which can appear in files from different OS.
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        // Remove newline characters.
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        // Append the cleaned line to the main sequence string.
        sequence += line;
    }

    inFile.close();
    return sequence; // Return the complete sequence.
}
