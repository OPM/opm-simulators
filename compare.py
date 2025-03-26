import re
import numpy as np
import sys
import glob
import argparse

def extract_vector(name, lines):
    """
    Extracts vector data from lines corresponding to the given matrix name.
    Returns a numpy array of the vector data.
    """
    vector_data = []
    capturing = False

    for line in lines:
        if line.startswith(name):
            capturing = True

        elif capturing and line.startswith("dune"):
            # Stop capturing when we hit the next matrix
            break

        elif capturing and len(line) > 0:
            numbers = re.findall(r'[-+]?\d*\.\d+e[-+]?\d+|[-+]?\d+\.\d+|[-+]?\d+', line)
            if len(numbers) > 0:
                vector_data.append([float(num) for num in numbers])

    return np.array(vector_data)


def extract_matrix(name, lines):
    """
    Extracts matrix data from lines corresponding to the given matrix name.
    Returns a numpy array of the matrix data.
    """
    matrix_data = []
    capturing = False

    for line in lines:
        if line.startswith(name):
            capturing = True

        elif capturing and line.startswith("dune"):
            # Stop capturing when we hit the next object
            break

        elif capturing and len(line) > 0:
            # Continue capturing matrix data
            numbersStrings = line.split('|')
            total = len(numbersStrings)
            even_numbers_reverse = list(range(total if total % 2 == 0 else total - 1, -1, -2))
            for i in even_numbers_reverse:
                numbersStrings.pop(i)
            for numbersString in numbersStrings:
                # Convert to floats if they are valid
                numbers = re.findall(r'[-+]?\d*\.\d+e[-+]?\d+|[-+]?\d+\.\d+|[-+]?\d+', numbersString)
                if len(numbers) > 0:
                    matrix_data.append([float(num) for num in numbers])

    return np.array(matrix_data)

def are_similar(obj1, obj2, verbose, precision=0.99, atol=5e-8): # todo: is 0.99 and 5e-8 ok??
    """
    Compares three matrix objects with a given precision.
    Returns True if all are similar, False otherwise.
    """
    if not (obj1.shape == obj2.shape):
        print("The objects do not have the same size!\n")
        return False
    if not np.allclose(obj1, obj2, rtol=1 - precision, atol=atol):
        if verbose:
            difference_mask = ~np.isclose(obj1, obj2, rtol=1 - precision, atol=atol)

            # Print or extract the locations of differences
            diff_indices = np.argwhere(difference_mask)
            print(f"Differences found at these indices:\n{diff_indices}")
            for idx in diff_indices:
                idx_tuple = tuple(idx)  # Convert to tuple for indexing
                print(f"At index {idx_tuple}: obj1 = {obj1[idx_tuple]}, obj2 = {obj2[idx_tuple]}")
        return False
    return True

def compare_across_files(file1_path, file2_path, verbose):
    """
    Extracts and compares the 'duneBGlobal_' matrix across three files.
    """
    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
        lines1 = file1.readlines()
        lines2 = file2.readlines()

    # Extract the matrices
    matrix_names = ['duneBGlobal_', 'duneCGlobal_', 'duneD_']
    for matrix_name in matrix_names :
        print(f"Will check matrix {matrix_name}...")
        matrix1 = extract_matrix(matrix_name, lines1)
        matrix2 = extract_matrix(matrix_name, lines2)

        # Compare matrices
        if are_similar(matrix1, matrix2, verbose):
            print(f"{matrix_name}: Matrices for {file1_path} and {file2_path} are similar.")
        else:
            print(f"\033[91m{matrix_name}: Matrices for {file1_path} and {file2_path} are different.\033[0m")

    vector_names = ['duneResidual']
    for vector_name in vector_names:
        print(f"Will check vector {vector_name}...")
        vector1 = extract_vector(vector_name, lines1)
        vector2 = extract_vector(vector_name, lines2)

        # Compare vectors
        if are_similar(matrix1, matrix2, verbose):
            print(f"{vector_name}: Vectors for {file1_path} and {file2_path} are similar.")
        else:
            print(f"\033[91m{vector_name}: Vectors for {file1_path} and {file2_path} are different.\033[0m")

def str_to_bool(s):
    return s.lower() in ['true', '1', 'yes', 'y']

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare matrices across multiple files. First argument is a verbosity flag, second argument are files.")
    parser.add_argument("verbosity", help="Verbosity flag")
    parser.add_argument("files", nargs='+', help="Paths to files (wildcards supported)")

    args = parser.parse_args()

    # Expand wildcards and get a sorted list of unique file paths
    file_list = sorted(set(file for pattern in args.files for file in glob.glob(pattern)))

    if len(file_list) < 2:
        print("Error: At least two files are required for comparison.")
        exit(1)

    verbosity = str_to_bool(args.verbosity)

    print(f"Comparing {len(file_list)} files: {file_list}")

    # Compare all unique pairs
    for i in range(len(file_list)):
        for j in range(i + 1, len(file_list)):
            print(f"\n\n======================================================== COMPARE FILES {file_list[i]} AND {file_list[j]} =======================================================")
            compare_across_files(file_list[i], file_list[j], verbosity)
