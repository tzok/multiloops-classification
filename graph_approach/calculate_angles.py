import matplotlib.pyplot as plt
import random
import numpy as np
import json
import pandas as pd
from rich import print
from rich.pretty import pprint
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os




def angle_between_vectors_on_plane(A, B, normal_vector):
    """
    Oblicza kąt pomiędzy rzutami dwóch wektorów na płaszczyznę.

    Parametry:
    A, B -- wektory w przestrzeni 3D 
    normal_vector -- wektor normalny do płaszczyzny 

    Zwraca:
    kąt w radianach pomiędzy rzutami wektorów na płaszczyznę
    """
    n = normal_vector / np.linalg.norm(normal_vector)

    A_proj = A - np.dot(A, n) * n
    B_proj = B - np.dot(B, n) * n

    # Sprawdzenie czy któryś z rzutów nie jest wektorem zerowym
    if np.linalg.norm(A_proj) < 1e-10 or np.linalg.norm(B_proj) < 1e-10:
        return None 

    cos_theta = np.dot(A_proj, B_proj) / (np.linalg.norm(A_proj) * np.linalg.norm(B_proj))

    cos_theta = max(min(cos_theta, 1.0), -1.0)

    angle_rad = np.arccos(cos_theta)

    cross_product = np.cross(A_proj, B_proj)
    if np.dot(cross_product, normal_vector) < 0:
        angle_rad = -angle_rad 

    return np.degrees(angle_rad)

def distance_3d_numpy(point1, point2):
    """
    Calculate euclidian distance between two points in 3D space.
    """


    return np.linalg.norm(point2 - point1)

def height_from_triangle_point(p1, p2, p3, from_point='p1'):
    """
    Oblicza wysokość (odległość) z podanego punktu (p1, p2 lub p3) na prostą wyznaczoną przez dwa pozostałe.

    Parametry:
    - p1, p2, p3: współrzędne punktów w 3D
    - from_point: 'p1', 'p2' lub 'p3' — punkt, z którego wychodzi wysokość

    Zwraca:
    - wysokość jako float
    """
    # Zamień na NumPy
    p1, p2, p3 = np.array(p1), np.array(p2), np.array(p3)

    if from_point == 'p1':
        point = p1
        line_start, line_end = p2, p3
    elif from_point == 'p2':
        point = p2
        line_start, line_end = p1, p3
    elif from_point == 'p3':
        point = p3
        line_start, line_end = p1, p2
    else:
        raise ValueError("from_point musi być 'p1', 'p2' lub 'p3'")

    line_vec = line_end - line_start
    point_vec = point - line_start

    cross_prod = np.cross(line_vec, point_vec)
    height = np.linalg.norm(cross_prod) / np.linalg.norm(line_vec)
    return height

def main():
    # Load JSON data from file
    labeled_df = pd.read_csv('..//data/3way_junctions_labeled.csv')
    # Split the 'Junction_ID' column into 'pdb' and 'id'
    labeled_df[['pdb', 'id']] = labeled_df['Junction_ID'].str.split('_', expand=True)
    
    folder_path = './forgi_graph_files'
    json_files = [f for f in os.listdir(folder_path) if f.endswith('.json')]
    features = []
    families = []
    structure_names = []
    for json_file in json_files:
        file_path = os.path.join(folder_path, json_file)
        pdb_id = json_file.split('.')[0]

        with open(file_path, 'r') as file:
            cg_model = json.load(file)
        if len(cg_model.keys()) == 0 or '3' not in cg_model or len(cg_model['3']) == 0:
            print(f"File {json_file} has no 3-way junctions.")
            continue
        for junction in cg_model['3']:
            try:
                temp_features = []
                structure_names.append(json_file)
                ranges = [junction[0]['m_residues']] + [(junction[1]['m_residues'])]+ [(junction[2]['m_residues'])]
                temp_family = ''
                for _, x in labeled_df[labeled_df['pdb'] == pdb_id].iterrows():
                    list_of_matches = []
                    list_of_matches.append(any(int(start) <= int(x['id']) <= int(end) for start, end in ranges))
   
                    list_of_matches.append(any(int(start) <= int(x['id']) <= int(end) for start, end in ranges))
                    if any(list_of_matches):
                        temp_family = x['Family']
                if temp_family != '':
                    families.append(temp_family)
                else:
                    families.append(random.choice(['A', 'B', 'C']))
                    # families.append('Unknown')
    


                p1 = np.array(junction[0]['next_coords'])
                p2 = np.array(junction[1]['next_coords'])
                p3 = np.array(junction[2]['next_coords'])

                s_p1 = np.array(junction[0]['stem']['next']['coords'])
                s_p2 = np.array(junction[1]['stem']['next']['coords'])
                s_p3 = np.array(junction[2]['stem']['next']['coords'])

                # Calculate the normal vector to the plane defined by p1, p2, and p3
                normal = np.cross(p2 - p1, p3 - p1)
                normal /= np.linalg.norm(normal)

                centroid = np.mean(np.array([p1,p2,p3]), axis=0)
                a_c = centroid - p1
                a = s_p1 - p1
                b_c = centroid - p2
                b = s_p2 - p2
                c_c = centroid - p3
                c = s_p3 - p3

                temp_features.append(angle_between_vectors_on_plane(a,a_c, normal))
                temp_features.append(angle_between_vectors_on_plane(b,b_c, normal))
                temp_features.append(angle_between_vectors_on_plane(c,c_c, normal))

                temp_features.append(height_from_triangle_point(p1, p2, p3, from_point='p1'))
                temp_features.append(height_from_triangle_point(p1, p2, p3, from_point='p2'))
                temp_features.append(height_from_triangle_point(p1, p2, p3, from_point='p3'))
                features.append(temp_features)
            except Exception as e:
                print(f"Error processing junction in file {json_file}: {e}")
                families.pop(0) 
                continue
     

    # Perform PCA analysis on features
    # Convert features to a NumPy array
    features_array = np.array(features)
    print(features_array.shape)
    print(len(families))
    # Standardize the features
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features_array)

    # Apply PCA
    pca = PCA(n_components=2)  # Reduce to 2 principal components
    features_pca = pca.fit_transform(features_scaled)

    # Print explained variance ratio
    print("Explained variance ratio:", pca.explained_variance_ratio_)

    # Combine features and families into a single DataFrame
    dataset = pd.DataFrame(features, columns=['angle1', 'angle2', 'angle3', 'height1', 'height2', 'height3'])
    dataset['family'] = families
    saved_path = '../data/3way_junctions_features_graphs.csv'
    dataset.to_csv(saved_path, index=False)

    # Plot PCA-transformed features
    plt.figure(figsize=(8, 6))
    plt.scatter(features_pca[:, 0], features_pca[:, 1], s=100, alpha=0.7, edgecolor='k')
    unique_families = list(set(families))
    cmap = plt.get_cmap('tab10')
    colors = cmap(np.linspace(0, 1, len(unique_families)))
    family_color_map = {family: colors[i] for i, family in enumerate(unique_families)}

    for i, family in enumerate(families):
        plt.scatter(features_pca[i, 0], features_pca[i, 1], s=100, alpha=1, edgecolor='k', color=family_color_map[family], label=family if family not in plt.gca().get_legend_handles_labels()[1] else "")


    plt.legend(title="Families", loc="best")
    plt.title('PCA of Features for 3-way junctions')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()