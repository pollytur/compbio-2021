from processing import file_extractor
import matplotlib.pyplot as plt
from processing import file_extractor
from processing import filter_out
from occurence import dictionnary_writer
from tqdm import tqdm
import numpy as np

file_path_initial = "dic_initial.txt"
file_path_variant = "dic_variant.txt"
output_file_path_initial = "dic_initial_filtered.txt"
output_file_path_variant = "dic_variant_filtered.txt"

def comparison_files(file_path_1, file_path_2):
    initial_dic = file_extractor(file_path_1)
    variant_dic = file_extractor(file_path_2)
    res = []
    for key in variant_dic.keys():
        if key not in initial_dic.keys():
            res.append(key)
    return res

def comparison_dictionnaries(dic_initial, dic_variant):
    res = []
    for key in dic_variant.keys():
        if key not in dic_initial.keys():
            res.append(key)
    return res

def hamming(string_1, string_2):
    assert len (string_1) == len (string_2)
    return sum(chaine1 != chaine2 for chaine1, chaine2 in zip(string_1, string_2))

def opti_hamming(string_1, string_2):
    assert len (string_1) == len (string_2)
    err_counter = 0
    res = [0, string_2]
    for i in range(len(string_1)):
        if string_1[i] != string_2[i]:
            err_counter += 1
            res[0] = i
        if err_counter >= 2:
            res = [-1, string_2]
            break
    return res


def hamming_filter(difference_dic, initial_dic):
    res = []
    for key in difference_dic:
        for initial_key in initial_dic.keys():
            if hamming(key, initial_key) == 1:
                res.append(key)
    return res

def fast_hamming(difference_dic, initial_dic):
    one_hamming = defaultdict(tuple)
    letters = ['A', 'T', 'G', 'C']
    for key in difference_dic:
        more_than_one = 0
        more_than_one_per_elem = np.zeros(10)
        for i in range(0, len(key)):
            for letter in letters:
                temp_key = list(key)
                temp_key[i] = letter
                temp_key = ''.join(temp_key)
                if temp_key in initial_dic:
                    # that way we keep the last match only
                    one_hamming[key] = (temp_key, i, key[i], letter)
                    more_than_one += 1
                    more_than_one_per_elem[i] += 1
        if more_than_one >= 0:
            print(f'`FOR {key} it was {more_than_one} changes')
            print('shit = ', more_than_one_per_elem)
    return one_hamming


if __name__ == "__main__":
    difference_dic = comparison_files(output_file_path_initial, output_file_path_variant)
    dic = file_extractor(file_path_initial)
    res = hamming_filter(difference_dic, dic)
    print(res)

    

    

