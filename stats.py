from processing import file_extractor
import matplotlib.pyplot as plt

file_path_1 = "outputs/dic.txt"
file_path_2 = "outputs/filtered.txt"

def draw_hist(dic, xlim):
    """
    plots a hist of the dic with some lim on the x axis
    """
    plt.hist(dic.values(), range(1, xlim))
    plt.show()


if __name__ == "__main__":
    initial_dic = file_extractor(file_path_1)
    filtered_dic = file_extractor(file_path_2)
    xlim = max(initial_dic.values())
    #draw_hist(filtered_dic, xlim)
    draw_hist(initial_dic, xlim)
    