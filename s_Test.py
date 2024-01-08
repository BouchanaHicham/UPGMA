import pandas as pd
from pandas import DataFrame
from Bio import Phylo
from io import StringIO
import numpy as np

v = []  # Initialize an empty list to store intermediate values during UPGMA
mat_distance = [
    [0, 5, 2, 8],
    [5, 0, 5, 7],
    [2, 5, 0, 8],
    [8, 7, 8, 0]
]
liste = ["A", "B", "C", "D"]
tree = []  # Initialize an empty list to store tree information
cluster_labels = []  # Initialize an empty list to store cluster labels


def tableau(t):
    # Function to convert a lower triangular matrix to a symmetric matrix
    for i in range(len(t)):
        for j in range(i):
            t[i][j] = t[j][i]
    return t


def arbre_à_newick(tree):
    # Convert the tree information to Newick format and visualize the tree
    List1 = []
    Tree1 = []
    splits = []
    for item in tree:
        splits = item.split('+')  # Split clusters if '+' is found
        temp = ''
        for split in splits:
            if temp != '':
                temp += ', '
            if len(split) == 2:
                temp += '(' + split[0] + ', ' + split[1] + ')'
            elif len(split) == 3:
                temp += '((' + split[0] + ', ' + split[1] + ')' + ', ' + split[2] + ')'
            else:
                temp += split
        List1.append(temp)
    tmp = ''
    tmp += '((' + List1[0] + '),' + List1[1] + ')'
    Tree1.append(tmp)
    print()
    print("------------|""Arbre guide :""|---------------")
    print()
    print(Tree1)

    print()
    print("------------|""Plot :""|---------------")
    print()
    handle = StringIO(Tree1[0])
    tree2 = Phylo.read(handle, "newick")
    Phylo.draw(tree2)
    Phylo.draw_ascii(tree2)


def UPGMA(distance):
    global v
    global cluster_labels

    print()
    print("------------|""Original matrix:""|---------------")
    print()
    print(distance)
    cluster_labels = distance.columns
    i = 1
    print()
    while True:
        min_val = distance.replace(0, np.NaN).min(skipna=True).min()
        v.append(min_val / 2)
        print(v)

        if len(distance.columns) == 2:
            for col in distance.columns:
                found = False
                for tree_element in tree:
                    if col == tree_element.replace('+', ''):
                        found = True
                if not found:
                    tree.append(col)
            break

        temp = distance.replace(0, np.NaN).idxmin(skipna=True)

        for row, col in temp.items():
            if distance.at[row, col] == min_val:
                min_row = row
                min_col = col
                break

        print(tree, min_row, min_col)

        if min_row not in tree:
            if min_col not in tree:
                tree.append(min_row + min_col)
                new_col_row_name = min_row + min_col
            else:
                x = tree[tree.index(min_col)] + min_row
                tree[tree.index(min_col)] = x
                new_col_row_name = min_col + min_row
        elif min_col not in tree:
            tree[tree.index(min_row)] += min_col
            new_col_row_name = min_row + min_col
        else:
            if tree.index(min_row) < tree.index(min_col):
                tree[tree.index(min_row)] += "+" + min_col
                del tree[tree.index(min_col)]
                new_col_row_name = min_row + min_col
            else:
                tree[tree.index(min_col)] += "+" + min_row
                del tree[tree.index(min_row)]
                new_col_row_name = min_col + min_row

        new_matrix = distance.drop([min_row, min_col], axis=0).drop([min_row, min_col], axis=1)
        new_matrix = pd.DataFrame(new_matrix)

        # Check if the new column and row should be added at the end or at the beginning
        if min_row in new_matrix.index and min_col in new_matrix.index:
            new_matrix[new_col_row_name] = 0
            inde = [new_col_row_name]
            temp = pd.DataFrame(data=0, columns=new_matrix.columns, index=inde)
            new_matrix = new_matrix.append(temp)

        else:
            new_matrix.insert(0, new_col_row_name, 0.0)
            inde = [new_col_row_name]
            temp = pd.DataFrame(data=0, columns=new_matrix.columns, index=inde)
            new_matrix = pd.concat([temp, new_matrix], axis=0)

        for column in distance:
            if column != min_row and column != min_col:
                avg = distance.at[min_row, column] + distance.at[min_col, column]
                avg = avg / 2.0
                new_matrix.at[new_col_row_name, column] = avg
                new_matrix.at[column, new_col_row_name] = avg

        print()
        print("------------|""Number :", i, "|---------------")
        print()
        distance = new_matrix
        i += 1

    arbre_à_newick(tree)
    return distance


tableau(mat_distance)
distance = DataFrame(mat_distance, columns=liste, index=liste)
distance = UPGMA(distance)