import math

x = 31
y = 15
z = 67

print("X div y is " + str(x/y))
print("X mod y is " + str(x//y))

for i in range(0, 61):
    print("The ceiling of " + str(i) + " div 15 is " + str(math.ceil(i/15)))