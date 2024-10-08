import re
import numpy as np
import sys

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

def are_similar(obj1, obj2, precision=0.99, atol=5e-8): # todo: is 0.99 and 5e-8 ok??
    """
    Compares three matrix objects with a given precision.
    Returns True if all are similar, False otherwise.
    """
    if not (obj1.shape == obj2.shape):
        print("The objrices to not have the same size!")
        return False
    if not np.allclose(obj1, obj2, rtol=1 - precision, atol=atol):
        difference_mask = ~np.isclose(obj1, obj2, rtol=1 - precision, atol=atol)

        # Print or extract the locations of differences
        diff_indices = np.argwhere(difference_mask)

        print(f"Differences found at these indices:\n{diff_indices}")
        for idx in diff_indices:
            idx_tuple = tuple(idx)  # Convert to tuple for indexing
            print(f"At index {idx_tuple}: obj1 = {obj1[idx_tuple]}, obj2 = {obj2[idx_tuple]}")
        return False
    return True

def compare_across_files(file1_path, file2_path, file3_path):
    """
    Extracts and compares the 'duneBGlobal_' matrix across three files.
    """
    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2, open(file3_path, 'r') as file3:
        lines1 = file1.readlines()
        lines2 = file2.readlines()
        lines3 = file3.readlines()

    # Extract the duneBGlobal_ matrix from all three files
    matrix_names = ['duneBGlobal_', 'duneCGlobal_', 'duneD_']
    for matrix_name in matrix_names :
        print(f"Will check matrix {matrix_name}...")
        matrix1 = extract_matrix(matrix_name, lines1)
        matrix2 = extract_matrix(matrix_name, lines2)
        matrix3 = extract_matrix(matrix_name, lines3)

        # Compare matrices across the three files
        if are_similar(matrix1, matrix2): 
            print(f"{matrix_name}: Matrices 1 and 2 are similar.")
        else:
            print(f"{matrix_name}: Matrices 1 and 2 are different.")
        if are_similar(matrix1, matrix3):
            print(f"{matrix_name}: Matrices 1 and 3 are similar.")
        else:
            print(f"{matrix_name}: Matrices 1 and 3 are different.")
        if are_similar(matrix2, matrix3):
            print(f"{matrix_name}: Matrices 2 and 3 are similar.\n")
        else:
            print(f"{matrix_name}: Matrices 2 and 3 are different.\n")
        
    vector_names = ['duneResidual']
    for vector_name in vector_names:
        print(f"Will check vector {vector_name}...")
        vector1 = extract_vector(vector_name, lines1)
        vector2 = extract_vector(vector_name, lines2)
        vector3 = extract_vector(vector_name, lines3)

        # Compare matrices across the three files
        if are_similar(vector1, vector2): 
            print(f"{vector_name}: Vectors 1 and 2 are similar.")
        else:
            print(f"{vector_name}: Vectors 1 and 2 are different.")
        if are_similar(vector1, vector3):
            print(f"{vector_name}: Vectors 1 and 3 are similar.")
        else:
            print(f"{vector_name}: Vectors 1 and 3 are different.")
        if are_similar(vector2, vector3):
            print(f"{vector_name}: Vectors 2 and 3 are similar.\n")
        else:
            print(f"{vector_name}: Vectors 2 and 3 are different.\n")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compare the duneBGlobal_ matrix in three files.")
    parser.add_argument("file1", help="Path to the first file")
    parser.add_argument("file2", help="Path to the second file")
    parser.add_argument("file3", help="Path to the third file")

    args = parser.parse_args()
    compare_across_files(args.file1, args.file2, args.file3)
