#include "FastaParser.h"

// Implementation of the readSequence method for FastaParser.
std::string FastaParser::readSequence(const std::string& filename) {
    // Attempt to open the specified file for reading.
    std::ifstream inFile(filename);
    if (!inFile) {
        // If the file cannot be opened, print an error message to stderr and return an empty string.
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return "";
    }

    std::string sequence_data = ""; // Initialize an empty string to store the concatenated sequence.
    std::string line;          // String to hold each line read from the file.

    // Read the file line by line.
    while (std::getline(inFile, line)) {
        // Skip empty lines or header lines (which start with '>').
        if (line.empty() || line[0] == '>') {
            continue;
        }
        // Remove carriage return characters, which can appear in files from different OS.
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        // Remove newline characters (though getline usually handles this, it's a safeguard).
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        // Append the cleaned line to the main sequence string.
        sequence_data += line;
    }
    // Close the file stream. Although it's done automatically on destruction, explicit closing is good practice.
    inFile.close();
    return sequence_data; // Return the complete sequence.
}

// Implementation of the writeSequence method for FastaParser.
bool FastaParser::writeSequence(const std::string& filename, const std::string& header, const std::string& sequence, size_t line_width) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
        return false;
    }

    outfile << ">" << header << std::endl;

    if (line_width == 0) { // Write sequence in a single line
        outfile << sequence << std::endl;
    } else { // Write sequence with line breaks
        for (size_t i = 0; i < sequence.length(); i += line_width) {
            outfile << sequence.substr(i, line_width) << std::endl;
        }
    }

    outfile.close();
    return true;
}
