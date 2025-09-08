import numpy as np

def compare_files(file1_path, file2_path):
    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
        file1_lines = file1.readlines()
        file2_lines = file2.readlines()

    max_lines = max(len(file1_lines), len(file2_lines))
    differences = []

    for i in range(max_lines):
        line1 = file1_lines[i].strip() if i < len(file1_lines) else "<no line>"
        line2 = file2_lines[i].strip() if i < len(file2_lines) else "<no line>"

        if line1 != line2:
            differences.append((i + 1, line1, line2))

    if differences:
        print("Differences found:")
        for line_num, l1, l2 in differences:
            print(f"Line {line_num}:\n  File1: {l1}\n  File2: {l2}")
        return False

    return True

def compare_values(file1_path, file2_path):

    data = np.genfromtxt(file1_path, dtype=str)
    test = (data[:]).astype(float)

    data = np.genfromtxt(file2_path, dtype=str)
    ref = (data[:]).astype(float)

    return np.all(np.isclose(test, ref))

pass_test= np.full(3, True)
pass_test[0] = compare_files('gas_species_names.txt', 'references/gas_species_names_ref.txt')
pass_test[1] = compare_files('aero_species_props.txt', 'references/aero_species_props_ref.txt')
pass_test[2] = compare_values('output.txt', 'references/output_ref.txt')

assert(np.all(pass_test))

