# Projet de Modélisation Surfacique

Sujet choisi : reconstruction de maillage.

Implémentation de la méthode décrite dans l'article de H. Hoppe, T. DeRose, T. Duchamp, J. McDonald, W. Stuetzle : Surface reconstruction from unorganized points. SIGGRAPH'92 (1992)


## Installation et compilation
```{bash}
> mkdir build
> cd build
> cMake ..
> make
```    

## Exécutables principaux

- pour reconstruire un maillage à partir d'un ensemble de points `sphere.data` dans le dossier `data` en utilisant 4 voisins:
```{bash}
> ./reconstruction ../data/sphere.data 4
```
Les paramètres de cet exécutable sont :
    - file_path : le chemin vers le nuage de points à reconstruire
    - K : nombres de voisins utilisés dans l'algorithme
    - rho : la densité du maillage (pour les maillages ouverts)

- pour visualiser le nuage de points contenu dans le fichier `sphere.data` dans le dossier `data`:
```{bash}
> ./visu ../data/sphere.data
```
Le paramètre de cet exécutable est :
    - file_path : le chemin vers le nuage de points à visualiser

## Exécutables secondaires

- pour générer un fichier `sphere.data` dans le dossier `data` à partir d'un fichier `.off` contenu dans le dossier `models`:
```{bash}
> python3 data/sphere.data
```
- pour générer une sphère bruitée:
```{bash}
> python3 bruit.py
```
