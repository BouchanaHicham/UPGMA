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
    
    def get_parent(self, target_label, parent=None):
        if self.label == target_label:
            return parent
        for child in self.children:
            result = child.get_parent(target_label, self)
            if result:
                return result
        return None

def find_min_distance(matrix):
    min_distance = float('inf')
    min_i, min_j = -1, -1

    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix[i])):
            if matrix[i][j] < min_distance:
                min_distance = matrix[i][j]
                min_i, min_j = i, j

    return min_i, min_j, min_distance

def update_distance_matrix(matrix, cluster_i, cluster_j, labels):
    new_cluster_label = f"({labels[cluster_i]}/{labels[cluster_j]})"

    # Update the matrix by merging clusters i and j into a new cluster
    matrix[cluster_i] = (matrix[cluster_i] + matrix[cluster_j]) / 2
    matrix[:, cluster_i] = matrix[cluster_i]

    # Set the diagonal element of the new cluster to 0
    matrix[cluster_i][cluster_i] = 0

    # Remove the row and column corresponding to cluster j
    matrix = np.delete(matrix, cluster_j, axis=0)
    matrix = np.delete(matrix, cluster_j, axis=1)

    # Update the labels list
    labels[cluster_i] = new_cluster_label
    del labels[cluster_j]

    return matrix, labels

def upgma(distance_matrix):
    clusters = [[i] for i in range(len(distance_matrix))]
    labels = [chr(ord('A') + i) for i in range(len(distance_matrix))]  # A, B, C, ...
    print("Clusters: ",clusters)
    print("labels:" ,labels)
    print(" --------- Initial_Distance_Matrix --------- ")
    print(distance_matrix)
    print(" -------------------------------- ")
    Head = Node("Head", 0.0)
    current = Head
    prev = Head
    while len(clusters) > 1:
        i, j, min_distance = find_min_distance(distance_matrix)

        # Print merging information only if both i and j are within valid indices
        if 0 <= i < len(labels) and 0 <= j < len(labels):
            print(f"Merged clusters {labels[i]} and {labels[j]} with distance {min_distance}. New cluster: {labels[i]}")
            prev = current

            if(prev.label != "Head"):
                print("prev value: ",prev.value)
                current.add_child(labels[i], (min_distance/2) - prev.value)
                current.add_child(labels[j], min_distance/2)
                current.add_child("---", min_distance/2)
            else:
                current.add_child(labels[i], min_distance/2)
                current.add_child(labels[j], min_distance/2)
                current.add_child("---", min_distance/2)

            
            current = current.get_child("---")
            
            # Print the entire tree
            Head.print_tree()
            
            print(" - - - - - ")
            '''
            print("Left_Node",Left_Node.label,"value",Left_Node.value)
            print("Right_Node",Right_Node.label,"value",Right_Node.value)
            print("Merged_Node",Merged_Node.label,"value",Merged_Node.value)
            '''
        else:
            print(f"One of the clusters is out of index. Merging not performed.")

        # Update the cluster list and distance matrix
        clusters[i].extend(clusters[j])
        del clusters[j]

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

    return labels[0]

# Example usage
distance_matrix = np.array([[0.0, 5.0, 2.0, 8.0],
                            [5.0, 0.0, 5.0, 7.0],
                            [2.0, 5.0, 0.0, 8.0],
                            [8.0, 7.0, 8.0, 0.0]])

upgma_result = upgma(distance_matrix)
