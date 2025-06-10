// Decompressor.cpp

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

//
// ORGD‐style helper definitions (no Compressor calls)
//
namespace ORGD {
    // Read an entire FASTA (preserving Ns & case), skipping the first line.
    std::string readSeq(const std::string& path) {
        std::ifstream fin(path);
        if (!fin) throw std::runtime_error("ORGD::readSeq: cannot open " + path);
        std::string line, seq;
        std::getline(fin, line);           // skip header
        while (std::getline(fin, line)) {
            seq += line;
        }
        return seq;
    }

    // Read a FASTA, uppercase everything and drop 'N'
    std::string readrefSeq(const std::string& path) {
        std::ifstream fin(path);
        if (!fin) throw std::runtime_error("ORGD::readrefSeq: cannot open " + path);
        std::string line, seq;
        std::getline(fin, line);           // skip header
        while (std::getline(fin, line)) {
            for (char c : line) {
                char C = std::toupper(static_cast<unsigned char>(c));
                if (C != 'N') seq.push_back(C);
            }
        }
        return seq;
    }

    std::stringstream reconstruct(const std::string& intermediateFile,
                                  const std::string& reference) {
        std::ifstream fin(intermediateFile);
        if (!fin) throw std::runtime_error("ORGD::reconstruct: cannot open " + intermediateFile);

        // Skip the first 4 header lines
        std::string line;
        for (int i = 0; i < 4 && std::getline(fin, line); ++i) {}

        std::string out;
        long long prev_end = -1;
        while (std::getline(fin, line)) {
            auto comma = line.find(',');
            if (comma != std::string::npos) {
                int delta  = std::stoi(line.substr(0, comma));
                int length = std::stoi(line.substr(comma + 1));
                long long begin = prev_end + delta;
                if (prev_end != -1 && (begin < 0 || begin + length > (long long)reference.size())) {
                    throw std::runtime_error(
                    "ORGD::reconstruct: computed segment out of range: begin="
                    + std::to_string(begin)
                    + " length=" + std::to_string(length)
                    + " ref_size=" + std::to_string(reference.size()));
                }
                out.append(reference.substr(fmax(begin, 0), length));
                prev_end = begin + length - 1;
            } else {
                out.append(line);
            }
        }
        return std::stringstream(std::move(out));
    }

    // Write the final FASTA header+sequence
    void write(const std::string& outFile, const std::string& text) {
        std::ofstream fout(outFile);
        if (!fout) throw std::runtime_error("ORGD::write: cannot open " + outFile);
        fout << text;
    }
}
// end ORGD helpers

void Decompressor::decompress(
    const std::string& referencePath,
    const std::string& archivePath,
    const std::string& outputFaPath
) {
    // 1) Extract .7z into a temp folder
//    std::string unzipDir = outputFaPath + "_unzip";
//    std::filesystem::remove_all(unzipDir);
//    std::filesystem::create_directories(unzipDir);
//    {
//        std::string sevenza = "./7za";
//        std::string cmd = sevenza + " e \"" + archivePath + "\" -o\"" + unzipDir + "\" -aos";
//        if (std::system(cmd.c_str()) != 0) {
//            throw std::runtime_error("Decompressor: failed to extract archive: " + archivePath);
//        }
//    }

    // 2) Locate the single intermediate file
//    std::string intermediate;
//    for (auto &e : std::filesystem::directory_iterator(unzipDir)) {
//        if (e.is_regular_file()) {
//            intermediate = e.path().string();
//            break;
//        }
//    }
//    if (intermediate.empty()) {
//        throw std::runtime_error("Decompressor: no file found in " + unzipDir);
//    }

    // 3) Parse header lines: 0=meta, 1=line_length, 2=LOWERCASE, 3=N_POSITIONS
    std::vector<Position> L_list, N_list;
    std::string meta_data;
    int line_length = 0;

    {
        std::ifstream fin(archivePath);
        if (!fin) throw std::runtime_error("Decompressor: cannot open " + archivePath);
        std::string line;
        for (int pass = 0; pass < 4 && std::getline(fin, line); ++pass) {
            if (pass == 0) {
                meta_data = line + "\n";
            } else if (pass == 1) {
                line_length = std::stoi(line);
            } else {
                std::istringstream iss(line);
                int delta, len, prev = 0;
                iss >> line; // Read the prefix (L: or N:)
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

    // 4) Load reference, stripping Ns if we have N_list
    std::string reference = N_list.empty()
        ? ORGD::readSeq(referencePath)
        : ORGD::readrefSeq(referencePath);
    if (N_list.empty()) {
        std::transform(reference.begin(), reference.end(),
                       reference.begin(), ::toupper);
    }

    // 5) Reconstruct the core sequence
    std::stringstream core_ss = ORGD::reconstruct(archivePath, reference);
    std::string target = core_ss.str();

    // 6) Re-insert N runs
    if (!N_list.empty()) {
        std::stringstream tmp;
        int ti = 0, ri = 0;
        for (auto &pos : N_list) {
            while (ri < pos.getStartInTarget() && ti < (int)target.size()) {
                tmp << target[ti++];
                ri++;
            }
            for (int k = pos.getStartInTarget(); k < pos.getEndInTarget(); ++k) {
                tmp << 'N';
                ri++;
            }
        }
        while (ti < (int)target.size()) {
            tmp << target[ti++];
        }
        target = tmp.str();
    }

    // 7) Restore lowercase runs
    for (auto &pos : L_list) {
        for (int i = pos.getStartInTarget();
             i < pos.getEndInTarget() && i < (int)target.size();
             ++i)
        {
            target[i] = std::tolower(target[i]);
        }
    }

    std::ostringstream wrapped;
    for (int i = 0; i < (int)target.size(); ++i) {
        if (i && i % line_length == 0) wrapped << "\n";
        wrapped << target[i];
    }

//    // Add newline every 50 characters
//    std::ostringstream final;
//    // Add meta data
//    final << meta_data;
//    for (size_t i = 0; i < wrapped.str().length(); i += 50) {
//        final << wrapped.str().substr(i, 50);
//        if (i + 50 < wrapped.str().length()) final << "\n";
//    }
//
//    wrapped << "\n";

    ORGD::write(outputFaPath, meta_data + wrapped.str() + "\n");

//    std::filesystem::remove_all(unzipDir);
    std::cout << "Decompressed FASTA written to " << outputFaPath << "\n";
}
