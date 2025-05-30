#ifndef FASTA_PARSER_H
#define FASTA_PARSER_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

// FastaParser class: Responsible for reading biological sequences from FASTA formatted files.
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
};

#endif // FASTA_PARSER_H
