from processing import file_extractor
import matplotlib.pyplot as plt
import math as m

file_path_1 = "./outputs/dic_initial.txt"
file_path_2 = "./outputs/dic_variant.txt"
# file_path_2 = "./outputs/filtered.txt"

def draw_hist(dic, xlim, name):
    """
    plots a hist of the dic with some lim on the x axis
    """
    # bi = m.ceil((dic.values.max() - dic.values.min())/5)
    plt.hist(dic.values(), range(1,xlim))
    plt.title(name)
    plt.show()


if __name__ == "__main__":
    initial_dic = file_extractor(file_path_1)
    variant_dic = file_extractor(file_path_2)
    # xlim = max(variant_dic.values())
    # xlim = 500
    plt.figure()
    plt.hist(initial_dic.values(), range(1,max(initial_dic.values())))
    plt.title("Initial sequence")
    
    plt.figure()
    plt.hist(variant_dic.values(), range(1,max(variant_dic.values())))
    plt.title("Variant sequence")

    plt.show()
    