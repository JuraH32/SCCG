# SCCG - Sequence Compression by Comparison to a Genome

## Compilation

Compile the project using a C++17 compatible compiler (like g++):

```bash
g++ main.cpp src/*.cpp -o sccg -std=c++17 -O3
```
(Note: The executable name is changed to `sccg` from `main` for clarity. Adjust if you prefer `main`.)

## Usage

The program operates in three modes: `encode`, `decode`, and `check`.

### Encode Mode

Compresses a target FASTA file against a reference FASTA file.

**Syntax:**
```bash
./sccg encode <reference.fa> <target.fa> <output_dir> [k_mer k0 L_seg search_range_global T1 T2]
```

**Arguments:**
*   `<reference.fa>`: Path to the reference FASTA file.
*   `<target.fa>`: Path to the target FASTA file to be compressed.
*   `<output_dir>`: Directory where the compressed output (and intermediate files) will be stored. The final compressed file will be named `<target_filename>.7z` inside this directory.

**Behavior:**
*   If `<reference.fa>` and `<target.fa>` are single files, it compresses the target against the reference.
*   If `<reference.fa>` and `<target.fa>` are directories, it iterates through FASTA files in the target directory and attempts to compress each against a correspondingly named file in the reference directory.

**Example:**
```bash
./sccg encode data/hg18/chr1.fa data/hg17/chr1.fa ./compressed_output
./sccg encode data/hg18/ data/hg17/ ./compressed_output_genome
```

### Decode Mode

Decompresses a target file that was previously compressed by SCCG, using the original reference.

**Syntax:**
```bash
./sccg decode <reference.fa> <encoded_archive.7z> <decoded_output.fa>
```

**Arguments:**
*   `<reference.fa>`: Path to the original reference FASTA file used during encoding.
*   `<encoded_archive.7z>`: Path to the compressed `.7z` archive file.
*   `<decoded_output.fa>`: Path where the reconstructed FASTA file will be saved.

**Behavior:**
*   If `<reference.fa>` and `<encoded_archive.7z>` are single files, it decompresses the archive.
*   If `<reference.fa>` and `<encoded_archive.7z>` are directories, it iterates through `.7z` files in the archive directory and attempts to decompress each against a correspondingly named reference file (stripping `.7z` from the reference filename if present). The output directory for decoded files is specified by `<decoded_output.fa>` (which should be a directory path in this case).

**Example:**
```bash
./sccg decode data/hg18/chr1.fa compressed_output/chr1.fa.7z decoded_output/chr1_decoded.fa
./sccg decode data/hg18/ compressed_output_genome/ decoded_output_genome/
```

### Check Mode

Compares two FASTA files to see if their sequence contents are identical, ignoring newline characters.

**Syntax:**
```bash
./sccg check <sequence1.fa> <sequence2.fa>
```

**Arguments:**
*   `<sequence1.fa>`: Path to the first FASTA file.
*   `<sequence2.fa>`: Path to the second FASTA file.

**Behavior:**
*   If `<sequence1.fa>` and `<sequence2.fa>` are single files, it compares them.
*   If `<sequence1.fa>` and `<sequence2.fa>` are directories, it iterates through files in the first directory and compares each against a correspondingly named file in the second directory.

**Example:**
```bash
./sccg check data/hg18/chr1.fa decoded_output_genome/chr1.fa
./sccg check data/hg18/ decoded_output_genome
```
