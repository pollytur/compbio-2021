from processing import file_extractor
from processing import filter_out
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from occurence import dictionnary_writer

file_path_1 = "outputs/dic_initial_filtered.txt"
file_path_2 = "outputs/dic_variant_filtered.txt"
file_path_initial = "outputs/dic_initial.txt"
file_path_variant = "outputs/dic_variant.txt"


def comparison(file_path_1, file_path_2):
    initial_dic = file_extractor(file_path_1)
    variant_dic = file_extractor(file_path_2)
    print(len(initial_dic.keys()))
    print(len(variant_dic.keys()))
    diffs = {k: variant_dic[k] for k in variant_dic if not k in initial_dic}
    print("diffs", len(diffs.keys()))
    return diffs


def draw_hist(dic, xlim):
    print("dic = ", dic)
    plt.plot(range(0, xlim), dic)
    plt.show()


def variance_over_threshold():
    filtered_dic_initial = file_extractor(file_path_initial)
    filtered_dic_variant = file_extractor(file_path_variant)
    diff_old = defaultdict(lambda: 0)
    N = 30
    diffs = np.zeros(N)
    for threshold in range(0, N, 1):
        print(threshold)
        filtered_dic_initial = filter_out(filtered_dic_initial, threshold)
        filtered_dic_variant = filter_out(filtered_dic_variant, threshold)
        diff_new = {k: filtered_dic_variant[k] for k in filtered_dic_variant if not k in filtered_dic_initial}
        step_diff = {k: diff_old[k] for k in diff_old if not k in diff_new}
        diffs[threshold] = len(diff_new.keys())
        print("Number of old values not presented in the new one: ", len(step_diff.keys()))
        diff_old = diff_new
        if threshold == 12:
            dictionnary_writer(diff_old, "./outputs/diff_old.txt")
            dictionnary_writer(diff_new, "./outputs/diff_new.txt")
    draw_hist(diffs, N)


def hamming(string_1, string_2):
    assert len(string_1) == len(string_2)
    for i in range(0, len(string_1)):
        if string_1[i] != string_2[i]:
            return i
    return -1


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
        for i in range(0, len(key)):
            for letter in letters:
                temp_key = list(key)
                temp_key[i] = letter
                temp_key = ''.join(temp_key)
                if temp_key in initial_dic:
                    one_hamming[key] = (temp_key, i, key[i], letter)
                    more_than_one += 1
        if more_than_one > 1:
            print(f'`FOR {key} it was {more_than_one} changes')
    return one_hamming


if __name__ == "__main__":
    different_keys = comparison(file_path_1, file_path_2)
    dictionnary_writer(different_keys, "./outputs/diff_thld1.txt")
    print(len(different_keys))
    rr = fast_hamming(different_keys, file_extractor(file_path_1))
    print('Length of my new dic', len(rr))
    print(rr)
    # variance_over_threshold()
