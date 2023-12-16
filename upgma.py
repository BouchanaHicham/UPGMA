import numpy as np

def find_min_distance(matrix):
    min_distance = float('inf')
    min_i, min_j = -1, -1

    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix[i])):
            if matrix[i][j] < min_distance:
                min_distance = matrix[i][j]
                min_i, min_j = i, j

    return min_i, min_j, min_distance

def update_distance_matrix(matrix, cluster_i, cluster_j):
    new_cluster = (cluster_i + cluster_j) / 2

    # Update the matrix by merging clusters i and j into a new cluster
    matrix[cluster_i] = (matrix[cluster_i] + matrix[cluster_j]) / 2
    matrix[:, cluster_i] = matrix[cluster_i]

    # Remove the row and column corresponding to cluster j
    matrix = np.delete(matrix, cluster_j, axis=0)
    matrix = np.delete(matrix, cluster_j, axis=1)

    return matrix, new_cluster

def upgma(distance_matrix):
    clusters = [[i] for i in range(len(distance_matrix))]

    while len(clusters) > 1:
        i, j, min_distance = find_min_distance(distance_matrix)

        # Update the cluster list and distance matrix
        clusters[i].extend(clusters[j])
        del clusters[j]

        distance_matrix, new_cluster = update_distance_matrix(distance_matrix, i, j)

        print(f"Merged clusters {i} and {j} with distance {min_distance}. New cluster: {new_cluster}")
        print("Cluster List:", clusters)
        print("Distance Matrix:")
        print(distance_matrix)
        print()

    return clusters[0]

# Example usage
distance_matrix = np.array([[0, 5, 2, 8],
                            [5, 0, 5, 7],
                            [2, 5, 0, 8],
                            [8, 7, 8, 0]])

upgma_result = upgma(distance_matrix)
