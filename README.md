# SCCG

[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](./LICENSE)

## Original Paper  
This is a C++ implementation of the SCCG algorithm. You can find the paper here: [https://academic.oup.com/bioinformatics/article/35/12/2058/5165377], and the original Java implementation is available at: [https://github.com/jhchen5/SCCG].

## Course Information  
This project was completed as part of the *Bioinformatics* course at FER: [https://www.fer.unizg.hr/en/course/enbio1].


## Prepare data:
1. Download the data zip file from: https://drive.google.com/file/d/1x5dOx_hL22pVWZHQxDl6-V-DCAzistKj/view?usp=sharing

2. Move the data zip file to the repository with the name `data.zip`.

   The directory structure should look like this:
   ```
   SCCG
   ├── data.zip
   ├── prepare.sh
   ├── orgc.cpp
   ├── orgd.cpp
   |── README.md
   └── 7za
   ```

3. Run the following command to prepare the data:
    ```bash
    bash prepare.sh
    ```


## Compression:

1.  **Compile the code:**

    ```bash
    g++ orgc.cpp -o orgc -std=c++11
    ```

2.  **Run the compression tool:**

    ```bash
    ./orgc <path_to_ref_genome> <path_to_target> <output_dir>
    ```

    -   `<path_to_ref_genome>`:  Path to the directory containing the reference genome files (e.g., `./data/hg17`).
    -   `<path_to_target>`: Path to the directory containing the target genome files (e.g., `./data/hg38`).
    -   `<output_dir>`: Path to the directory where the compressed results will be stored (e.g., `result`).

    **Example:**

    ```bash
    ./orgc ./data/hg17 ./data/hg38 result
    ```

3.  **Locate the results:**

    The compressed output will be saved in the `<output_dir>/result` directory.

## Decompression:

1.  **Compile the code:**

    ```bash
    g++ orgd.cpp -o orgd -std=c++11
    ```

2.  **Run the decompression tool:**

    ```bash
    ./orgd <path_to_ref_genome> <path_to_target> <input_dir>
    ```

    -   `<path_to_ref_genome>`: Path to the directory containing the reference genome files (e.g., `./data/hg17`).
    -   `<path_to_target>`: Path to the directory containing the target genome files (e.g., `./data/hg38`).
    -   `<input_dir>`: Path to the directory containing the *compressed* results (e.g., `result`).

    **Example:**

    ```bash
    ./orgd ./data/hg17 ./data/hg18 result
    ```
