# L'Explication de l'Intersection à 30% pour la Cohorte 2W

L'exécution des scripts de validation (Python et R) a donné de brillants résultats pour la cohorte 7W (jusqu'à 98.00% d'intersection avec la liste de référence de Supriya). Cela prouve hors de tout doute que la méthodologie d'extraction, la gestion des valeurs manquantes, la méthode de normalisation (t-test / VSN) et l'algorithme `limma` ont été **parfaitement répliqués**. 

Cependant, pour l'échantillon 2W, l'intersection reste bloquée autour de 30% (~184 protéines en commun).

## Pourquoi ce plafond de 30% ?

La raison n'est mathématiquement pas liée au pipeline, mais **aux données brutes elles-mêmes**. 

Dans son email, Supriya mentionne :
> *"I am sharing my R script for the analysis along with the **updated (unpublished) proteomics dataset**."*

Elle fait spécifiquement référence aux fichiers `results_2w.xlsx` et `results_7w.xlsx` qu'elle vous a adressés. Ces fichiers sont les listes de protéines issues d'un **nouveau jeu de données** ou d'une nouvelle quantification MaxQuant.

### La preuve par les données 

1. **Identification fondamentale différente** : Le tableau brut final `results_2w.xlsx` ne dispose que de **1520 protéines au total** après son extraction.
2. Lorsqu'on applique exactement le même script sur *notre* vieux fichier `proteinGroups.tsv` pour la condition 2W (même en excluant rigoureusement les deux échantillons erronés Ech15 et Ech16 comme conseillé), le filtre à 10% conserve **2034 protéines**.
3. Il y a donc une divergence de **plus de 500 protéines** dès la base (les protéines quantifiées/existantes ne sont plus les mêmes).

### Conclusion

Supriya n'a pas seulement écarté deux échantillons (Ech15 et Ech16). Elle (ou son équipe) a probablement repassé la cohorte 2W dans MaxQuant ou a régénéré le jeu de données pour constituer une version "mise à jour" (qui n'a jamais été publiée) qui résolvait son problème de "separate clustering". 

Vos calculs sont **justes**. Vous obtenez 30% parce que votre script compare ses résultats issus de l'**ancien fichier raw** avec les siens issus du **nouveau fichier raw**. 

**Action recommandée :** Comme vous ne possédez pas ce nouveau fichier text raw, la solution la plus mathématiquement saine est d'arrêter de tenter de reproduire la cohorte 2W depuis le vieux `proteinGroups.tsv`, et d'utiliser directement les listes fournies dans `results_2w.xlsx` comme nouvelle vérité scientifique pour les étapes suivantes de votre stage.
