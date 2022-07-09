# Generates strain maps in the material based on the diffraction data acquired using
# Precession Electron Diffraction (PED)
# May also be used for 4D-STEM data
# The first part of the algorithm filters the PED data
# The second part calculates the distance in diffraction patterns
# The third part generates strain maps

import time
from os import remove, path

import numpy as np
import plotly.express as px
from csv import writer
from ctypes import c_double
from hyperspy.api import load
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from multiprocessing import Array, Pool
import tqdm
from numpy import sqrt, array, ndenumerate, arange, percentile, mean
from numpy.ctypeslib import as_array
from pandas import DataFrame
from PIL import Image, ImageTk
from pixstem.api import PixelatedSTEM
import tkinter as tk
from tkinter import filedialog
from skimage.feature import match_template
from scipy.ndimage import gaussian_filter
from skimage.restoration import estimate_sigma


# global variables
file = None  # input file
distances = None  # distances for all peak points in image
single_values = None  # input values
curr_func = None  # current function
shared_array = None
point1 = None  # user-selected point 1
point2 = None  # user-selected point 2

# changes curr_func to user selected function (used in buttons in UI)
# displays instructions in UI
def set_curr_func(func_name):
    global curr_func, file, single_values
    curr_func = str(func_name)
    entry.delete(0, tk.END)
    if curr_func == "to_csv":
        if file is None:
            label1['text'] = "Please load a file before saving data.\n"
        elif single_values is None:
            label1['text'] = "Please analyze the file before saving data.\n"
        else:
            entry.bind("<Return>", get_entry)
            label1['text'] = "Please enter the path of the file you want to save to in the\ntext box " \
                             "provided then press Enter.\n "
    elif curr_func == "analysis":
        if file is None:
            label1['text'] = "Please load a file before starting analysis.\n"
        else:
            entry.bind("<Return>", get_entry)
            label1['text'] = "Please enter the number of rows and columns you would like to " \
                             "analyze,\nas integers, separated by spaces. Press Enter when ready.\n"
    elif curr_func == 'heat map':
        if file is None:
            label1['text'] = "Please load a file before creating a heat map.\n"
        elif distances is None:
            label1['text'] = "Please analyze the file before creating a heat map.\n"
        else:
            label1['text'] = 'Please enter a strain-free distance value for comparison. \nFor reference, ' \
                             f'average distance of the dataset is {round(np.average(single_values), 3)}\n' \
                             'If the distribution is bimodal or multimodal, reference the histogram \nfor the desired '\
                             'distance.'
            entry.bind("<Return>", get_entry)


# executes current function upon keyboard event (used in set_curr_func)
def get_entry(event):
    global curr_func
    if curr_func == "load_file":
        entry.unbind("<Return>")
        load_file(entry.get())
    elif curr_func == "analysis":
        entry.unbind("<Return>")
        start_analysis(entry.get())
    elif curr_func == "to_csv":
        entry.unbind("<Return>")
        to_csv(entry.get())
    elif curr_func == 'heat map':
        entry.unbind("<Return>")
        heat_map(entry.get())


# opens file explorer for user to select file
# assigns selected file to global file variable
def load_file(filename=None):
    global file
    label1['text'] = "Loading file...\n"
    input_file = filedialog.askopenfilename()
    root.update()
    try:
        file = load(input_file)
        label1['text'] = label1['text'] + "File loaded.\n"
    except:
        label1['text'] = label1['text'] + "Error loading. Please check path and try again.\n"
    entry.delete(0, tk.END)
    # entry.unbind("<Return>")


# calculates distance between two points
def distance(x1, y1, x2, y2):
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2))


# finds coordinate of center spot in image
def find_center(im, peak, xy2):
    center = (xy2[0], xy2[1])
    minimum = 144
    for (i, j) in ndenumerate(peak):
        for (x, y) in j:
            length = len(im)
            dist = distance(xy2[0], xy2[1], y, x)
            # d = distance(length/2, length/2, b, a)
            if int(xy2[0]*0.9) < int(x) < int(xy2[0]*1.1) and int(xy2[1]*0.9) < int(y) < int(xy2[1]*1.1)\
                    and dist < minimum:
                minimum = dist
                center = (y, x)
    return center


# applies filter to image
# performs template matching and finds peak points in diffraction spots
# calculates distances for all peak points in image
def multiprocessing_func(values):

    original = values[1]
    ####################################################################################################################
    # # FILTERS
    # Gaussian is the default filter
    # Uncomment nlm (non-local means) or wien (wiener) to use those filters

    sigma_est = mean(estimate_sigma(original, ))
    gaussian = gaussian_filter(original, 1.15 * sigma_est)

    # # patch_size for 580 - 1, patch_size for 144 = 3
    # nlm = denoise_nl_means(original, h=1.15*sigma_est, fast_mode=True, patch_size=1, patch_distance=6, )
    # wien = wiener(original, 5, 3)

    filtered = array(gaussian)
    ####################################################################################################################
    # # PIXSTEM

    # s = s.rotate_diffraction(0,show_progressbar=False)
    # st = s.template_match_disk(disk_r=5, lazy_result=False, show_progressbar=False)
    # peak_array = st.find_peaks(lazy_result=False, show_progressbar=False)
    # peak_array_com = s.peak_position_refinement_com(peak_array, lazy_result=False, show_progressbar=False)
    # s_rem = s.subtract_diffraction_background(lazy_result=False, show_progressbar=False)
    # peak_array_rem_com = s_rem.peak_position_refinement_com(peak_array_com, lazy_result=False, show_progressbar=False)
    ####################################################################################################################
    # MY METHOD

    # defines template and templates matches
    # spot for 580 - [265:320, 265:320]
    # spot for 144 - [65:80, 65:80]
    template = filtered[265:320, 265:320]
    result = match_template(filtered, template, pad_input=True)
    # only takes points greater than the threshold r-value
    tempList = []
    for i in range(len(result)):
        for j in range(len(result[i])):
            if result[i][j] > 0.87:  # correlation value
                tempList.append((i, j))

    # removes duplicate spots that are too close to each other
    peaks_list = []
    while len(tempList) > 0:
        j = 0
        temp = []
        point = tempList[0]
        while j < len(tempList):
            if distance(point[0], point[1], tempList[j][0], tempList[j][1]) < 25:  # min center dist can be changed
                temp.append(tempList[j])
                tempList.pop(j)
            else:
                j = j + 1
        max = 0
        if temp:
            pnt = temp[0]
        for j in range(len(temp)):
            if result[pnt[0]][pnt[1]] < result[temp[j][0]][temp[j][1]]:
                max = result[temp[j][0]][temp[j][1]]
                pnt = temp[j]
        peaks_list.append(pnt)
    peak_array_rem_com = [[], peaks_list]
    ####################################################################################################################
    center = find_center(filtered, peak_array_rem_com, [values[0][4], values[0][5]])
    # finds the specific spot and adding that distance to the array
    pos_distance = 0
    closest_point = center
    idx = 0
    length = len(filtered)

    # calculates distances for all peak points in image
    for (i, j) in ndenumerate(peak_array_rem_com):
        minimum = 999999
        for (x, y) in j:
            if 2 < y < length - 2 and 2 < x < length - 2:
                di = distance(center[0], center[1], y, x)
                idx += 1
            dis = distance(values[0][2], values[0][3], y, x)
            if dis < minimum and dis < length / 10:
                minimum = dis
                closest_point = (y, x)
    pos_distance = distance(closest_point[0], closest_point[1], center[0], center[1])
    return values[0][0], values[0][1], closest_point, pos_distance, center


# creates pop-up UI where user can select 2 points on the image
# calls analysis function with selected 2 points
def start_analysis(values=None):
    global file, curr_func, point1, point2

    def reset_coords():
        global point1, point2
        point1 = None
        point2 = None
        analysis_log['text'] = "Strain mapping: please click on the two points you would like to " \
                               "analyze from the diffraction pattern above.\n"

    def confirm_coords():
        global point1, point2
        if point1 is not None and point2 is not None:
            print("Selected points are ", point1, " and ", point2)
            analysis_log['text'] = analysis_log['text'] + "Starting analysis...\n"
            analysis(point1, values, point2)
            remove("temp.png")
            c2.unbind('<Button-1>')
            r.destroy()
            label1['text'] = label1['text'] + "Analysis complete.\n"
        else:
            reset_coords()

    def mouse_coords(event):
        s = PixelatedSTEM(file.inav[25, 25])
        length = len(array(s))
        global point1, point2
        if point1 is None:
            point1 = (int(event.x * length / 400), int(event.y * length / 400))  # get the mouse position from event
            analysis_log['text'] = analysis_log['text'] + str(point1[0]) + " " + str(point1[1]) + "\n"
            print("point1 is ", point1)
        elif point2 is None:
            point2 = (int(event.x * length / 400), int(event.y * length / 400))  # get the mouse position from event
            analysis_log['text'] = analysis_log['text'] + str(point2[0]) + " " + str(point2[1]) + "\n"
            print("point2 is ", point2)
        r.update()

    s = PixelatedSTEM(file.inav[25, 25])
    s.save("temp.png")
    img = Image.open("temp.png")
    img = img.resize((400, 400))
    img.save('temp.png')

    r = tk.Toplevel(root)

    c = tk.Canvas(r, height=640, width=840)
    c.pack()
    f = tk.Frame(r, bg='#FFFFFF')
    f.place(relwidth=1, relheight=1)
    analysis_log = tk.Message(f, bg='#FFFFFF', font=('Calibri', 15), anchor='nw', justify='left',
                              highlightthickness=0, bd=0, width=756)
    analysis_log.place(relx=0.05, rely=0.65, relwidth=0.9, relheight=0.2)

    reset_button = tk.Button(f, text='Reset', bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0,
                             bd=0, activebackground='#D4D4D4', activeforeground='#252525',
                             command=lambda: reset_coords(), pady=0.02, fg='#373737', borderwidth='2',
                             relief="groove")
    reset_button.place(relx=0.15, rely=0.90, relwidth=0.30, relheight=0.07)

    confirm_button = tk.Button(f, text='Confirm', bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0,
                             bd=0, activebackground='#D4D4D4', activeforeground='#252525',
                             command=lambda: confirm_coords(), pady=0.02, fg='#373737', borderwidth='2',
                             relief="groove")
    confirm_button.place(relx=0.55, rely=0.90, relwidth=0.30, relheight=0.07)

    c2 = tk.Canvas(r, width=400, height=400)
    c2.place(relx=0.5, anchor='n')
    img = ImageTk.PhotoImage(Image.open("temp.png"))
    c2.create_image(0, 0, anchor='nw', image=img)
    c2.bind('<Button-1>', mouse_coords)
    analysis_log['text'] = "Strain mapping: please click on the two points you would like to "\
                           "analyze from the diffraction pattern above.\n"

    r.mainloop()
    if path.exists("temp.png"):
        remove("temp.png")


# calls multiprocessing function to calculate all distances
# saves distances to new file
def analysis(pointxy1, values, pointxy2):
    global file, single_values, distances
    t = values.split(" ")
    ROW = int(t[0])
    COL = int(t[1])

    list = [[]]
    i = 0
    img_array = array(PixelatedSTEM(file))
    for r in range(ROW):
        for c in range(COL):
            list.append([])
            list[i].append([r, c, pointxy1[0], pointxy1[1], pointxy2[0], pointxy2[1]])
            list[i].append(img_array[r, c])
            i += 1
    del list[-1]
    del img_array
    shared_array_base = Array(c_double, ROW * COL)
    single_values = as_array(shared_array_base.get_obj())
    single_values = single_values.reshape(COL, ROW)

    shared_array = Array(c_double, ROW * COL * 50)
    distances = as_array(shared_array.get_obj())
    distances = distances.reshape(COL, ROW, 50)

    start_time = time.time()

    # Running the multiprocessing_func with a pool of processes to increase calculation time.
    # tqdm is a progress bar function that prints to the console.
    results = []
    pool = Pool(processes=None)
    for output in tqdm.tqdm(pool.imap_unordered(multiprocessing_func, list), total=len(list)):
        results.append(output)
        pass
    pool.close()

    for i in range(len(results)):
        single_values[results[i][1]][results[i][0]] = results[i][3]

    print('Calculation time:', round(time.time() - start_time, 2))

    entry.delete(0, tk.END)
    f = open("Distances", "w")
    w = writer(f)
    for i in distances:
        w.writerow(i)
    f.close()
    label1['text'] = label1['text'] + "File saved.\n"
    entry.delete(0, tk.END)
    # entry.unbind("<Return>")


# creates and saves data to a csv file
def to_csv(filename=None):
    global single_values
    f = open(filename, "w")
    w = writer(f)
    for i in single_values:
        w.writerow(i)
    f.close()
    label1['text'] = label1['text'] + "File saved.\n"
    entry.delete(0, tk.END)
    # entry.unbind("<Return>")


# creates bar chart pop-up UI
def bar_chart(INTERVAL=0.1):
    global distances
    if file is None:
        label1['text'] = "Please load a file before creating a bar chart.\n"
    elif distances is None:
        label1['text'] = "Please analyze the file before creating a bar chart.\n"
    else:
        label1['text'] = "Creating bar chart.\n(This might take several minutes depending on the size " \
                         "of data.)\n"
        root.update()
        dist = single_values.flatten()

        fig, a = plt.subplots(figsize=(6, 5.5))
        plt.xlabel('Distance from center peak', fontsize=10)
        plt.ylabel('Counts', fontsize=10)
        plt.title('Distance Counts', fontsize=10)
        # plt.bar(y_pos, counts, align='center', alpha=0.95) # creates the bar plot
        plt.hist(dist, bins=500)

        bar_chart_window = tk.Toplevel(root)
        bar_chart_window.geometry('600x600')
        chart_type = FigureCanvasTkAgg(plt.gcf(), bar_chart_window)
        chart_type.draw()
        chart_type.get_tk_widget().place(relx=0.0, rely=0.0, relwidth=1)


# calculates and returns thresholds for lower and upper outliers in a given array
def outlier(data):
    data = data.flatten()
    q1 = percentile(data, 25)
    q3 = percentile(data, 75)
    iqr = q3 - q1
    minimum = q1 - (1.5 * iqr)
    maximum = q3 + (1.5 * iqr)
    return minimum, maximum


# opens a heat map as a tab in browser
# creates pop-up UI that can bring up user-requested diffraction patterns
def heat_map(input_distance):
    import hyperspy.api as hs
    global single_values
    data = single_values.copy()
    d0 = float(input_distance)

    strain_values = []
    for i in range(len(data)):
        strain_values.append([])
        for j in range(len(data[i])):
            if data[i][j] == 0:
                strain_values[i].append(float('NaN'))
            else:
                strain_values[i].append((d0 / int(data[i][j])) - 1)

    df = DataFrame(strain_values, columns=arange(len(strain_values[0])), index=arange(len(strain_values)))
    print(df)
    fig = px.imshow(df, color_continuous_midpoint=0, zmin=-0.05, zmax=0.05, color_continuous_scale='turbo')
    fig.show()

    # creates diffraction pattern pop-up windows
    # def image_gallery(event):
    #     global file
    #     values = e.get().split(" ")
    #     x0 = int(values[0])
    #     y0 = int(values[1])
    #     x1 = int(values[2])
    #     y1 = int(values[3])
    #     indexx = x1
    #
    #     for x in range(x1 - x0 + 1):
    #         indexy = y1
    #         for y in range(y1 - y0 + 1):
    #             s = PixelatedSTEM(hs.signals.Signal2D(file.inav[indexx, indexy]))
    #             st = s.template_match_ring(r_inner=1, r_outer=6, lazy_result=True, show_progressbar=False)
    #             peak_array = st.find_peaks(method='dog', min_sigma=0.8, max_sigma=15, sigma_ratio=1.9,
    #                                        threshold=0.42, overlap=0.5, lazy_result=False, show_progressbar=True)
    #             s.add_peak_array_as_markers(peak_array)
    #             # plt.plot(s)
    #             s.plot()
    #             ax = s._plot.signal_plot.ax
    #             ax.set_xlabel("pixel(" + str(indexx) + "_" + str(indexy) + ")")
    #             # plt.title("pixel(" + str(indexx) + "_" + str(indexy) + ")")
    #             # plt.show()
    #             indexy -= 1
    #             y += 1
    #         indexx -= 1
    #         x += 1
    #     plt.show()

    # creates pop-up UI that calls image_gallery upon user input
    # bar_chart_window = tk.Toplevel(root)
    # bar_chart_window.geometry('500x400')
    # m = tk.Message(bar_chart_window, font=('Calibri', 15), highlightthickness=0, bd=0, justify='left')
    # m['text'] = "A new window should open displaying the heatmap created. If you would like to view specific " \
    #             "diffraction patterns, enter the starting x and the y value and the ending x and y value " \
    #             "separated by a space. Press Enter to display these diffraction patterns."
    # m.place(relx=0.05, rely=0.05, relwidth=0.9, relheight=0.5)
    # e = tk.Entry(bar_chart_window, bg='#F4F4F4', font=('Calibri', 15), justify='left', highlightthickness=0,
    #              bd=0, fg='#373737', borderwidth=2, relief="groove")
    # e.place(relx=0.1, rely=0.7, relwidth=0.8, relheight=0.1)
    # e.bind("<Return>", image_gallery)

    # bar_chart_window.mainloop()


# main menu UI
if __name__ == "__main__":

    HEIGHT = 700
    WIDTH = 800

    root = tk.Tk()

    canvas = tk.Canvas(root, height=HEIGHT, width=WIDTH)
    canvas.pack()
    frame = tk.Frame(root, bg='#FFFFFF')
    frame.place(relwidth=1, relheight=1)

    # TAMU MSEN logo
    image = Image.open('msen.png')
    image = image.resize((200, 40))
    img = ImageTk.PhotoImage(image)
    label1 = tk.Label(frame, image=img, bg='#FFFFFF')
    label1.place(relx=0.05, rely=0.05, anchor='w')

    # Menu Label
    label = tk.Label(frame, text='PED Strain Mapping', bg='#FFFFFF', font=('Times New Roman', 40), fg='#373737')
    label.place(relx=0.20, rely=0.1, relwidth=0.6, relheight=0.1)

    # Text Output box
    label1 = tk.Message(frame, bg='#F3F3F3', font=('Calibri', 15), anchor='nw', justify='left', highlightthickness=0,
                        bd=0, width=1500, fg='#373737', borderwidth=2, relief="groove")
    label1['text'] = "This program was originally designed by Aniket Patel and Aaron Barbosa \nand modified by " \
                     "Marcus Hansen and Ainiu Wang."
    label1.place(relx=0.1, rely=0.54, relwidth=0.8, relheight=0.32)

    # Entry box
    entry = tk.Entry(frame, bg='#F3F3F3', font=('Calibri', 15), justify='left', highlightthickness=0,
                     bd=0, width=1500, fg='#373737', borderwidth=2, relief="groove")
    entry.place(relx=0.1, rely=0.88, relwidth=0.8, relheight=0.05)

    # Buttons
    button = tk.Button(frame, text='Load File', bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0, bd=0,
                       activebackground='#D4D4D4', activeforeground='#252525',
                       command=lambda: load_file(), pady=0.02, fg='#373737', borderwidth='2',
                       relief="groove")
    button.place(relx=0.29, rely=0.22, relwidth=0.42, relheight=0.05)

    button1 = tk.Button(frame, text='Start Analysis', bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0, bd=0,
                        activebackground='#D4D4D4', activeforeground='#252525',
                        command=lambda: set_curr_func("analysis"), pady=0.02, fg='#373737', borderwidth='2',
                        relief="groove")
    button1.place(relx=0.29, rely=0.28, relwidth=0.42, relheight=0.05)

    button2 = tk.Button(frame, text='Create Bar Chart', bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0, bd=0,
                        activebackground='#D4D4D4', activeforeground='#252525',
                        command=lambda: bar_chart(), pady=0.02, fg='#373737', borderwidth='2',
                        relief="groove")
    button2.place(relx=0.29, rely=0.34, relwidth=0.42, relheight=0.05)

    button3 = tk.Button(frame, text='Create Strain Heat Map', bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0, bd=0,
                        activebackground='#D4D4D4', activeforeground='#252525',
                        command=lambda: set_curr_func('heat map'), pady=0.02, fg='#373737', borderwidth='2',
                        relief="groove")
    button3.place(relx=0.29, rely=0.46, relwidth=0.42, relheight=0.05)

    button4 = tk.Button(frame, text='Export Distance Data to .csv', bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0,
                        bd=0, activebackground='#D4D4D4', activeforeground='#252525',
                        command=lambda: set_curr_func("to_csv"), pady=0.02, fg='#373737', borderwidth='2',
                        relief="groove")
    button4.place(relx=0.29, rely=0.40, relwidth=0.42, relheight=0.05)

    root.mainloop()
    if path.exists("temp.png"):
        remove("temp.png")
