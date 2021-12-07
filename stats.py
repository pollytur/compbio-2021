from processing import file_extractor
import matplotlib.pyplot as plt

file_path_1 = "outputs/dic_initial.txt"
file_path_2 = "outputs/dic_initial_filtered.txt"

def draw_hist(dic, xlim):
    """
    plots a hist of the dic with some lim on the x axis
    """
    plt.hist(dic.values(), range(16587, 1676784))
    plt.show()

if __name__ == "__main__":
    initial_dic = file_extractor(file_path_1)
    #filtered_dic = file_extractor(file_path_2)
    xlim = 1676784  #max(initial_dic.values())
    #draw_hist(filtered_dic, xlim)
    draw_hist(initial_dic, xlim)
    