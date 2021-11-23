from processing import file_extractor

file_path_1 = "outputs/dic_initial_filtered.txt"
file_path_2 = "outputs/dic_variant_filtered.txt"

def comparison(file_path_1, file_path_2):
    initial_dic = file_extractor(file_path_1)
    variant_dic = file_extractor(file_path_2)
    print(len(initial_dic.keys()))
    print(len(variant_dic.keys()))
    return set(initial_dic) - set(variant_dic)

if __name__ == "__main__":
    print(len(comparison(file_path_1, file_path_2)))
