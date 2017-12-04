file = open("exampleCroped.data", "w")

with open("data/bear.data") as f:
    content = f.readlines()

i = 1
nb_lines = 0
step = 2
for x in content:
    if i%step == 0:
        file.write(x)
        nb_lines += 1
    i += 1

print("Number of lines croped from", i, "to", nb_lines)
