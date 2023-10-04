import random

num_of_bodies = 1000
new_file = open(str(num_of_bodies) + '_bodies.txt', 'w')
new_file.write("0 100")
for i in range(num_of_bodies):
    for j in range(7):
        num = random.randint(0,100)
        num = round(num/10, 1)
        if j == 6:
            new_file.write(" 1")
        else:
            new_file.write(" " + str(num))

