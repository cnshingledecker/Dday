import csv

with open('num_matrices.csv', 'r') as file2:
    reader = csv.reader(file2, delimiter=",")
    for line in reader:
        print(line)