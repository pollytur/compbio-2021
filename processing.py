from occurence import *

file_path_initial = "outputs/dic_initial.txt"
file_path_variant = "outputs/dic_variant.txt"
output_file_path_initial = "outputs/dic_initial_filtered.txt"
output_file_path_variant = "outputs/dic_variant_filtered.txt"


def file_extractor(file_path):
    """
    returns a dictionnary which is the same than the one created in occurence.py
    """
    d = defaultdict(int)
    with open(file_path, "r") as f:
        for line in f.readlines():
            S = line.split()
            d[S[0]] = int(S[1])
    return d

def filter_out(dic, threshold):
    """
    returns the dictionnary generated by file_extractor but filtered
    """
    return dict(filter(lambda occurence: int(occurence[1]) > threshold, dic.items()))


threshold = 20

if __name__ == "__main__":
    #as previously we write our filtered result in a file to eventually
    #process it later
    # Extracting the dictionary files and filtering out all sequences which appear only once
    D = file_extractor(file_path_initial)
    print("Filtering out errors for dictionary of 'original' sequence...")
    filtered_dic = filter_out(D, threshold)
    dictionary_writer(filtered_dic, output_file_path_initial)

    D = file_extractor(file_path_variant)
    print("Filtering out errors for dictionary of variant sequence...")
    filtered_dic = filter_out(D, threshold)
    dictionary_writer(filtered_dic, output_file_path_variant)
    print("Done.")