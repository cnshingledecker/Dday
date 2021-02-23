import csv

num1 = 0
num2 = 0
num3 = 0
with open('sum_iter_vars.csv', 'r') as csv_file:
    reader = csv.reader(csv_file)
    for line in reader:
        num1 = line[0]
        num2 = line[1]
        num3 = line[2]

print(num1)
print(num2)
print(num3)