#ifndef FASTA_PARSER_H
#define FASTA_PARSER_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm> // Required for std::remove

// FastaParser class: Responsible for reading and writing biological sequences from/to FASTA formatted files.
class FastaParser {
public:
    /**
     * @brief Reads a biological sequence from a FASTA file.
     *
     * This function opens the specified FASTA file, reads its content,
     * and concatenates all sequence lines into a single string.
     * Header lines (starting with '>') and empty lines are ignored.
     * It also removes newline and carriage return characters from sequence lines.
     *
     * @param filename The path to the FASTA file.
     * @return A string containing the concatenated sequence. Returns an empty string
     * if the file cannot be opened or if an error occurs.
     */
    static std::string readSequence(const std::string& filename);

    /**
     * @brief Writes a biological sequence to a FASTA file.
     *
     * @param filename The path to the FASTA file to be created/overwritten.
     * @param header The header string (without the leading '>').
     * @param sequence The biological sequence string.
     * @param line_width The maximum number of characters per line for the sequence (0 for no line breaks).
     * @return True if writing was successful, false otherwise.
     */
    static bool writeSequence(const std::string& filename, const std::string& header, const std::string& sequence, size_t line_width = 60);
};

#endif // FASTA_PARSER_H
