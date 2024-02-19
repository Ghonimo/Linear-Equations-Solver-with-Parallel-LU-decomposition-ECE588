# matrix_generator.py
# Author: Mohamed Ghonim
# Created: 02/18/2024
# Last Modified: 02/18/2024
# Function: This script generates a random matrix of a given dimension and saves it to a file in matrices/py_generated
# Usage: python matrix_generator.py <dimension>
# version: 0.1

import numpy as np
import sys
import os

def generate_matrix_and_save(dimension):

    # Define the target directory
    target_directory = os.path.join(os.path.dirname(__file__), "../matrices/py_generated")
    #target_directory = "../matrices"

    # Ensure the target directory exists
    os.makedirs(target_directory, exist_ok=True)
    
    # Generating a dimension x dimension matrix with random integers between -10 and 10
    matrix = np.random.randint(-10, 10, size=(dimension, dimension))
    
    # Generating the final row with random integers between -10 and 10
    final_row = np.random.randint(-10, 10, size=(dimension,))
    
    # Preparing the string representation and defining the file path
    matrix_str = f"{dimension}\n" + "\n".join(" ".join(map(str, row)) for row in matrix) + "\n" + " ".join(map(str, final_row))
    file_name = os.path.join(target_directory, f"{dimension}x{dimension}.txt")
    
    # Saving to file
    with open(file_name, 'w') as file:
        file.write(matrix_str)
    
    print(f"Matrix saved to {file_name}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <dimension>")
        sys.exit(1)
    
    try:
        dimension = int(sys.argv[1])
        if dimension <= 0:
            raise ValueError
        generate_matrix_and_save(dimension)
    except ValueError:
        print("Please provide a valid positive integer for the matrix dimension.")
