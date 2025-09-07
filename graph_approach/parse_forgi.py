import numpy as np
from rich import print
from rich.pretty import pprint
import json
import argparse


# class CGmodel:
#     def __init__(self, coords_key, coords_values, multi_list):
#         self.coords_key = coords_key
#         self.coords_values = coords_values
#         self.multi_list = multi_list



def parse_file(cg_file_path):
    coords_key = {}
    coords_values = {}
    multi_list = []
    m_length = {}
    s_residues = {}
    m_nucleotides = {}


    try:
        with open(cg_file_path, encoding='utf-8') as file:
            for f in file: 
                line = [el.replace('\n','') for el in f.split()]
                if not line or len(line) == 0:
                    continue
                if line[0] == 'define' and 'm' in line[1]:
                    if len(line)>=4:
                        m_length[line[1]] = int(line[3]) - int(line[2]) + 1
                        m_nucleotides[line[1]] = (int(line[2]),int(line[3]))
                    else:
                        m_length[line[1]] = 0
                        m_nucleotides[line[1]] = (0,0)
                if line[0] == 'define' and 's' in line[1]:
                    s_residues[line[1]] = [(int(line[2]), int(line[3])), (int(line[4]),int(line[5]))]

                if line[0] == 'coord' and ("s" in line[1] or "m" in line[1]):
                    key1 = (float(line[2]),float(line[3]),float(line[4]))
                    key2 = (float(line[5]),float(line[6]),float(line[7]))
                    coords_key.setdefault(key1, []).append(line[1])
                    coords_key.setdefault(key2, []).append(line[1])
                    coords_values[line[1]] =  [key1, key2]
                    if 'm' in line[1]:
                        multi_list.append(line[1])
     
    except (OSError, IOError) as e:
        print(f"Error reading file: {cg_file_path} - {e}")
        return None

    return coords_key, coords_values, multi_list, m_length, s_residues, m_nucleotides

def extract_junction(coords_key, coords_values, multi_list, m_length, s_residues, m_nucleotides):


    multi_set = set(multi_list)
    ml = {}
    while multi_set:
        curr = multi_set.pop()  # pobierz dowolny punkt startowy
        mle_temp = []
        visited = set()
        end = curr
        is_junction = True
        while is_junction:
            if curr in visited:
                break  # zapobiega zapętleniu
            visited.add(curr)

            mle = {"name": curr, "length": m_length[curr], "stem":{'previous':{},'next':{}}, 'm_residues': m_nucleotides[curr]}
            start_coords = coords_values[curr][0]
            for el in coords_key.get(start_coords, []):
                if "s" in el:
                    s_coords = coords_values[el]
                    s_coords_outside = s_coords[0] if s_coords[1] == start_coords else s_coords[1]
                    mle['stem']['previous'] = {'name':el, 'coords':s_coords_outside, 'residues':s_residues[el] if el in s_residues else None}
                if "m" in el and el != curr:
                    mle["previous"] = el
                    mle['previous_coords'] = start_coords
            
            end_coords = coords_values[curr][1]
                    
            for el in coords_key.get(end_coords, []):
                if "s" in el:
                    s_coords = coords_values[el]
                    s_coords_outside = s_coords[0] if s_coords[1] == end_coords else s_coords[1]
                    mle['stem']['next'] = {'name':el, 'coords':s_coords_outside, 'residues':s_residues[el] if el in s_residues else None}
                if "m" in el and el != curr:
                    mle["next"] = el
                    mle['next_coords'] = end_coords

            mle_temp.append(mle)

            # wyjdź, jeśli nie ma następnego elementu
            if "next" not in mle or mle["next"] == end:
                break

            curr = mle["next"]
            multi_set.discard(curr)

        n_junction = len(mle_temp)
        ml.setdefault(n_junction, []).append(mle_temp)


    return ml       
       # junctions['edges']



def main():
    # Example usage
    parser = argparse.ArgumentParser(description="Process FORGI file and save output as JSON.")
    parser.add_argument("file_path", type=str, help="Path to the FORGI file.")
    parser.add_argument("output_file", type=str, help="Path to save the output JSON file.")
    args = parser.parse_args()

    file_path = args.file_path
    output_file = args.output_file

    coords_key, coords_values, multi_list, m_length, s_residues, m_nucleotides = parse_file(file_path)
    ml = extract_junction(coords_key, coords_values, multi_list, m_length, s_residues, m_nucleotides)

    # Save the dictionary to a JSON file
    try:
        with open(output_file, 'w', encoding='utf-8') as json_file:
            json.dump(ml, json_file, ensure_ascii=False, indent=4)
        print(f"Dictionary saved to {output_file}")
    except (OSError, IOError) as e:
        print(f"Error saving JSON file: {output_file} - {e}")


if __name__ == "__main__":
    main()