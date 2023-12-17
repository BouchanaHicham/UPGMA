import numpy as np

class Node:
    def __init__(self, label, value):
        self.label = label
        self.value = value
        self.children = []

    def add_child(self, child_label, child_value):
        new_child = Node(child_label, child_value)
        self.children.append(new_child)
        return new_child

    def print_tree(self, indent=0):
        print("  " * indent + f"{self.label}: {self.value}")
        for child in self.children:
            child.print_tree(indent + 1)

    def get_child(self, target_label):
        for child in self.children:
            if child.label == target_label:
                return child
        return None
# ---------------------------------------------------------- [Min Distance] ----------------------------------------------------------
def find_min_distance(matrix):
    min_distance = float('inf') # positive infinity. btw im tired [I Hicham just built diff :p]
    min_i, min_j = -1, -1

    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix[i])):
            if matrix[i][j] < min_distance:
                min_distance = matrix[i][j]
                min_i, min_j = i, j

    return min_i, min_j, min_distance
# ---------------------------------------------------------- [Updating UPGMA Matrix] ----------------------------------------------------------

def update_distance_matrix(matrix, cluster_i, cluster_j, labels):
    new_cluster_label = f"({labels[cluster_i]}/{labels[cluster_j]})"

    # Update the matrix by merging clusters i and j into a new cluster
    matrix[cluster_i] = (matrix[cluster_i] + matrix[cluster_j]) / 2 # update the row at index cluster_i.
    matrix[:, cluster_i] = matrix[cluster_i]  # update the column at index [cluster_i] with the values from the updated row at index [cluster_i]

    # Set the diagonal element of the new cluster to 0
    matrix[cluster_i][cluster_i] = 0

    # Remove the row and column corresponding to cluster j
    matrix = np.delete(matrix, cluster_j, axis=0) #  removes the row at index cluster_j. axis=0 argument specifies the operation along the rows.
    matrix = np.delete(matrix, cluster_j, axis=1) #  removes the column at index cluster_j. axis=1 argument specifies the operation along the columns.

    # Update the labels list
    labels[cluster_i] = new_cluster_label
    del labels[cluster_j]

    return matrix, labels
# ---------------------------------------------------------- [UPGMA] ----------------------------------------------------------
def upgma(distance_matrix):
    clusters = [[i] for i in range(len(distance_matrix))]
    # clusters = [[0], [1], [2], [3]] in n = 4 for example, ill use this as indexes so i can later use with my labels ( alphabets )
    labels = [chr(ord('A') + i) for i in range(len(distance_matrix))]  # A, B, C, ... 
    # ord() => Alpha to Unicode character A=65 
    # chr() => Unicode to Alpha 65 = A :)

    print("Clusters: ",clusters) #just indexes to keep track
    print("labels:" ,labels)
    print(" --------- Initial_Distance_Matrix --------- ")
    print(distance_matrix)
    print(" -------------------------------- ")
    
    Head = Node("Head", 0.0) # Root
    current = Head 
    #prev = Head
    while len(clusters) > 1:
        i, j, min_distance = find_min_distance(distance_matrix)

        # Print merging information only if both i and j are within valid indices
        # So when we don't encouter an issue at the last matrix ( index out of range ) 
        if 0 <= i < len(labels) and 0 <= j < len(labels):
            print(f"Merged clusters {labels[i]} and {labels[j]} with distance {min_distance}. New cluster: {labels[i] + '/' +labels[j]}")
            prev = current

            if(prev.label != "Head"): # if this is not our first iteration, we get the value of the parent, and we use it to calc our merged labels (new - old)
                print("prev value: ",prev.value)
                current.add_child(labels[i], (min_distance/2) - prev.value)
                current.add_child(labels[j], min_distance/2)
                current.add_child("---", min_distance/2)
            else: # if this is our first iteration we calc normally
                current.add_child(labels[i], min_distance/2)
                current.add_child(labels[j], min_distance/2)
                current.add_child("---", min_distance/2)

            
            current = current.get_child("---") # we push through ( continue moving )
            
            # Print the entire tree each time
            # Head.print_tree()

            print(" - - - - - ")
            '''
            print("Left_Node",Left_Node.label,"value",Left_Node.value)
            print("Right_Node",Right_Node.label,"value",Right_Node.value)
            print("Merged_Node",Merged_Node.label,"value",Merged_Node.value)
            '''
        else:
            print(f"One of the clusters is out of index. Merging not performed.")

        # Update the cluster list and distance matrix
        clusters[i].extend(clusters[j]) # Clusters:  [[0], [1], [2], [3]] if i=[0] and j=[2] we extend the i with j => [[0,2],[1],[2],[3]]
        del clusters[j] # we remove the "Zayed =)" =>  [[0,2],[1],[3]]

        distance_matrix, labels = update_distance_matrix(distance_matrix, i, j, labels)

        # Ensure the indices i and j are still valid after updating the matrix
        i = min(i, j)
        j = max(i, j)

        

        print("Cluster List:", clusters)
        print("New Labels",labels)
        print("Distance Matrix:")
        print(distance_matrix)
        print("--------------------------------")
        print()
    Head.print_tree()
    return labels[0]
# ---------------------------------------------------------- [Execution] ----------------------------------------------------------
# Course Distance Matrix
distance_matrix = np.array([[0.0, 5.0, 2.0, 8.0],
                            [5.0, 0.0, 5.0, 7.0],
                            [2.0, 5.0, 0.0, 8.0],
                            [8.0, 7.0, 8.0, 0.0]])


'''
# Test Distanec Matrix
distance_matrix = np.array([[0.0, 17.0, 21.0, 27.0, 30.0],
                            [17.0, 0.0, 12.0, 18.0, 23],
                            [21.0, 12.0, 0.0, 14.0, 17.0],
                            [27.0, 18.0, 14.0, 0.0, 9.0],
                            [30.0, 23.0, 17.0, 9.0, 0.0]])
'''
upgma_result = upgma(distance_matrix)
