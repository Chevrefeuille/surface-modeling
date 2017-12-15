from random import *

fichier_original = open("data/sphere.data", "r")
fichier_bruite = open("data/sphere_bruitee.data", "w")

fichier_bruite.write(fichier_original.readline())
for line in fichier_original:
    ligne = line.split()
    X = normalvariate(0,0.04) + float(ligne[0])
    Y = normalvariate(0,0.04) + float(ligne[1])
    Z = normalvariate(0,0.04) + float(ligne[2])
    fichier_bruite.write(str(X) + " " + str(Y) + " " + str(Z) + "\n")

fichier_original.close()
fichier_bruite.close()
