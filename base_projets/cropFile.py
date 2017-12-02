file = open("exampleCroped.txt", "w")

with open("example.txt") as f:
    content = f.readlines()

i = 0
nb_lines = 0
step = 5
for x in content:
    if i%step == 0:
        file.write(x)
        nb_lines += 1
    i += 1

print("Number of lines croped from", i, "to", nb_lines)
