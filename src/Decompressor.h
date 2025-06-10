#include <string>

class Decompressor {
public:
    // Unpacks a single-chromosome .7z archive and reconstructs
    // the target FASTA from the reference.
    //
    // referencePath: path to the reference FASTA
    // archivePath:   path to the .7z archive holding one encoded chromosome
    // outputFaPath:  path where the decoded .fa should be written
    static void decompress(
        const std::string& referencePath,
        const std::string& archivePath,
        const std::string& outputFaPath
    );
};
