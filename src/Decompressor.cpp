#include "Decompressor.h"
#include "Position.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <algorithm>

namespace ORGD {
    std::string readSeq(const std::string& path) {
        std::ifstream fin(path);
        if (!fin) throw std::runtime_error("ORGD::readSeq: cannot open " + path);
        std::string line, seq;
        std::getline(fin, line);         
        while (std::getline(fin, line)) seq += line;
        return seq;
    }

    std::string readrefSeq(const std::string& path) {
        std::ifstream fin(path);
        if (!fin) throw std::runtime_error("ORGD::readrefSeq: cannot open " + path);
        std::string line, seq;
        std::getline(fin, line);         
        while (std::getline(fin, line)) {
            for (char c : line) {
                char C = std::toupper(static_cast<unsigned char>(c));
                if (C != 'N') seq.push_back(C);
            }
        }
        return seq;
    }

    std::stringstream reconstruct(const std::string& intermediateFile, const std::string& reference) {
        std::ifstream file(intermediateFile);
        if (!file.is_open()) throw std::runtime_error("ORGD::reconstruct: cannot open " + intermediateFile);

        std::string line;
        std::stringstream result;

        int prev_end = -1;
        int header_lines = 0;

        while (std::getline(file, line)) {
            if (header_lines < 4) {
                ++header_lines;
                continue;
            }

            size_t comma_pos = line.find(',');
            if (comma_pos != std::string::npos) {
                int delta = std::stoi(line.substr(0, comma_pos));
                int len = std::stoi(line.substr(comma_pos + 1));
                int begin = prev_end + delta + 1; // delta is from previous END
                int end = begin + len - 1;

                if (begin < 0 || end >= static_cast<int>(reference.size())) {
                    throw std::runtime_error("reconstruct: out-of-bounds access: begin=" + std::to_string(begin) + ", end=" + std::to_string(end));
                }

                result << reference.substr(begin, len);
                prev_end = end;
            } else {
                result << line;
            }
        }
        return result;
    }

    void write(const std::string& outFile, const std::string& text) {
        std::ofstream fout(outFile);
        if (!fout) throw std::runtime_error("ORGD::write: cannot open " + outFile);
        fout << text;
    }
}

void Decompressor::decompress(
    const std::string& referencePath,
    const std::string& archivePath,
    const std::string& outputFaPath
) {
    std::string unzipDir = outputFaPath + "_unzip";
    std::filesystem::remove_all(unzipDir);
    std::filesystem::create_directories(unzipDir);

    {
        std::string cmd = "./7za e \"" + archivePath + "\" -o\"" + unzipDir + "\" -aos";
        if (std::system(cmd.c_str()) != 0) {
            throw std::runtime_error("Decompressor: failed to extract archive: " + archivePath);
        }
    }

    std::string intermediate;
    for (auto &e : std::filesystem::directory_iterator(unzipDir)) {
        if (e.is_regular_file()) {
            intermediate = e.path().string();
            break;
        }
    }
    if (intermediate.empty()) {
        throw std::runtime_error("Decompressor: no file found in " + unzipDir);
    }

    // Step 1: Parse the header lines
    std::vector<Position> L_list, N_list;
    std::string meta_data;
    int line_length = 0;

    {
        std::ifstream fin(intermediate);
        if (!fin) throw std::runtime_error("Decompressor: cannot open " + intermediate);
        std::string line;
        for (int pass = 0; pass < 4 && std::getline(fin, line); ++pass) {
            if (pass == 0) {
                meta_data = line + "\n";
            } else if (pass == 1) {
                line_length = std::stoi(line);
            } else {
                std::istringstream iss(line);
                int delta, len, prev = 0;
                while (iss >> delta >> len) {
                    Position p;
                    p.setStartInTarget(prev + delta);
                    p.setEndInTarget(prev + delta + len);
                    prev = p.getEndInTarget();
                    if (pass == 2) L_list.push_back(p);
                    else           N_list.push_back(p);
                }
            }
        }
    }

    // Step 2: Read reference (with or without 'N')
    std::string reference = N_list.empty()
        ? ORGD::readSeq(referencePath)
        : ORGD::readrefSeq(referencePath);
    if (N_list.empty()) {
        std::transform(reference.begin(), reference.end(), reference.begin(), ::toupper);
    }

    // Step 3: Delta-decode from intermediate file
    std::stringstream core_ss = ORGD::reconstruct(intermediate, reference);
    std::string target = core_ss.str();

    // Step 4: Re-insert N-runs
    if (!N_list.empty()) {
        std::string updated;
        int t_idx = 0, ref_idx = 0;
        for (const auto& pos : N_list) {
            while (ref_idx < pos.getStartInTarget() && t_idx < (int)target.size()) {
                updated += target[t_idx++];
                ref_idx++;
            }
            for (int k = pos.getStartInTarget(); k < pos.getEndInTarget(); ++k) {
                updated += 'N';
                ref_idx++;
            }
        }
        while (t_idx < (int)target.size()) {
            updated += target[t_idx++];
        }
        target = updated;
    }

    // Step 5: Apply lowercase runs
    for (const auto& pos : L_list) {
        for (int i = pos.getStartInTarget(); i < pos.getEndInTarget() && i < (int)target.size(); ++i) {
            target[i] = std::tolower(target[i]);
        }
    }

    // Step 6: Wrap lines
    std::ostringstream wrapped;
    for (int i = 0; i < (int)target.size(); ++i) {
        if (i && i % line_length == 0) wrapped << "\n";
        wrapped << target[i];
    }
    wrapped << "\n";

    ORGD::write(outputFaPath, meta_data + wrapped.str());
    std::filesystem::remove_all(unzipDir);
    std::cout << "Decompressed FASTA written to " << outputFaPath << "\n";
}
