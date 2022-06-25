# Generates strain maps in the material based on the diffraction data acquired using
# Precession Electron Diffraction (PED)
# May also be used for 4D-STEM data
# The first part of the algorithm filters the PED data
# The second part calculates the distance in diffraction patterns
# The third part generates strain maps

from os import remove, path
import plotly.express as px
from csv import writer
from ctypes import c_double
from hyperspy.api import load
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import cm
from matplotlib import colors
import matplotlib.pyplot as plt
from multiprocessing import Array
from numpy import sqrt, array, ndenumerate, arange, percentile, linspace, zeros, mean
from numpy.ctypeslib import as_array
from pandas import DataFrame
from PIL import Image, ImageTk
from pixstem.api import PixelatedSTEM
from seaborn import heatmap
import tkinter as tk
from tkinter import filedialog
from skimage.feature import match_template
from scipy.ndimage import gaussian_filter
from skimage.restoration import estimate_sigma, denoise_nl_means
from scipy.signal import wiener

# global variables
file = None  # input file
distances = None  # distances for all peak points in image
single_values = None  # input values
curr_func = None  # current function


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
def find_center(im, peak):
    center = (352, 382)
    minimum = 144
    for (i, j) in ndenumerate(peak):
        for (x, y) in j:
            length = len(im)
            dist = distance(350, 380, y, x)
            # d = distance(length/2, length/2, b, a)
            if 340 < int(x) < 390 and 340 < int(y) < 370 and dist < minimum:
                minimum = dist
                center = (y, x)
    return center


# applies filter to image
# performs template matching and finds peak points in diffraction spots
# calculates distances for all peak points in image
def multiprocessing_func(values):
    global single_values, distances
    s = PixelatedSTEM(file.inav[values[0], values[1]])

    original = array(s)
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
        pnt = temp[0]
        for j in range(len(temp)):
            if result[pnt[0]][pnt[1]] < result[temp[j][0]][temp[j][1]]:
                max = result[temp[j][0]][temp[j][1]]
                pnt = temp[j]
        peaks_list.append(pnt)
    peak_array_rem_com = [[], peaks_list]
    ####################################################################################################################

    center = find_center(filtered, peak_array_rem_com)
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
                distances[values[1]][values[0]][idx] = round(di, 3)
                idx += 1
            dis = distance(values[2], values[3], y, x)
            if dis < minimum and dis < length / 10:
                minimum = dis
                closest_point = (y, x)
    pos_distance = distance(closest_point[0], closest_point[1], center[0], center[1])
    single_values[values[1]][values[0]] = round(pos_distance, 4)
    print(values[0], values[1], closest_point, pos_distance, center)


# creates pop-up UI where user can select a point on the image
# calls analysis function with selected point
def start_analysis(values=None):
    global file, curr_func
    pointxy = None
    if pointxy is None:
        def mouse_coords(event):
            s = PixelatedSTEM(file.inav[25, 25])
            length = len(array(s))
            pointxy = (int(event.x * length / 400), int(event.y * length / 400))  # get the mouse position from event
            analysis_log['text'] = analysis_log['text'] + str(pointxy[0]) + " " + str(pointxy[1]) + "\n"
            analysis_log['text'] = analysis_log['text'] + "Starting analysis...\n"
            r.update()
            analysis(pointxy, values)
            remove("temp.png")
            c2.unbind('<Button-1>')
            r.destroy()
            label1['text'] = label1['text'] + "Analysis complete.\n"

        s = PixelatedSTEM(file.inav[25, 25])
        s.save("temp.png")
        img = Image.open("temp.png")
        img = img.resize((400, 400), Image.ANTIALIAS)
        img.save('temp.png')

        r = tk.Toplevel(root)

        c = tk.Canvas(r, height=600, width=840)
        c.pack()
        f = tk.Frame(r, bg='#FFFFFF')
        f.place(relwidth=1, relheight=1)
        analysis_log = tk.Message(f, bg='#FFFFFF', font=('Calibri', 15), anchor='nw', justify='left',
                                  highlightthickness=0, bd=0, width=756)
        analysis_log.place(relx=0.05, rely=0.7, relwidth=0.9, relheight=0.2)
        c2 = tk.Canvas(r, width=400, height=400)
        c2.place(relx=0.5, anchor='n')
        img = ImageTk.PhotoImage(Image.open("temp.png"))
        c2.create_image(0, 0, anchor='nw', image=img)
        c2.bind('<Button-1>', mouse_coords)
        analysis_log['text'] = analysis_log['text'] + "Strain mapping: \nPlease click on the point you would like to "\
                                                      "analyze from the diffraction pattern above.\n"
        r.mainloop()
        if path.exists("temp.png"):
            remove("temp.png")


# calls multiprocessing function to calculate all distances
# saves distances to new file
def analysis(pointxy, values):
    global file, single_values, distances
    t = values.split(" ")
    ROW = int(t[0])
    COL = int(t[1])

    list = []
    for r in range(ROW):
        for c in range(COL):
            list.append((r, c, pointxy[0], pointxy[1]))

    shared_array_base = Array(c_double, ROW * COL)
    single_values = as_array(shared_array_base.get_obj())
    single_values = single_values.reshape(COL, ROW)

    shared_array = Array(c_double, ROW * COL * 50)
    distances = as_array(shared_array.get_obj())
    distances = distances.reshape(COL, ROW, 50)

    # calling multiprocessing function with a single thread (to be updated)
    for i in range(len(list)):
        multiprocessing_func(list[i])

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


# creates heat map (used by bar_chart function)
def heat_map_maker(minimum, maximum, parity=0):
    global single_values, distances

    if parity == 0:
        data = single_values.copy()
        df = DataFrame(data, columns=arange(len(data[0])), index=arange(len(data)))
        _, a = plt.subplots(figsize=(6, 5.5))
        chart1 = heatmap(df, cmap=cm.get_cmap("rainbow"), ax=a, vmin=minimum, vmax=maximum, square=True)
        return chart1.get_figure()
    else:
        data = zeros((len(single_values), len(single_values[0])), dtype=float)
        for i in range(len(distances)):
            for j in range(len(distances[i])):
                sum = 0
                num = 0
                for k in distances[i][j]:
                    if minimum < k < maximum:
                        sum += k
                        num += 1
                if num > 0:
                    data[i][j] = round(sum / num, 1)

        df = DataFrame(data, columns=arange(len(data[0])), index=arange(len(data)))
        _, a = plt.subplots(figsize=(6, 5.5))
        gray = cm.get_cmap('gray', 512)
        newcolors = gray(linspace(0.15, 0.85, 2048))
        white = array([255 / 256, 255 / 256, 255 / 256, 1])
        newcolors[:1, :] = white
        newcolors[2047:, :] = white
        newcmp = colors.ListedColormap(newcolors)
        chart = heatmap(df, cmap=newcmp, vmin=minimum, vmax=maximum, square=True)
        return chart.get_figure()


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
        plt.xlabel('Distance from center peek', fontsize=10)
        plt.ylabel('Counts', fontsize=10)
        plt.title('Distance Counts', fontsize=10)
        # plt.bar(y_pos, counts, align='center', alpha=0.95) # creates the bar plot
        plt.hist(dist, bins=500)

        # calls heat_map_maker function upon user input
        def scope_heat_map(event):
            values = e.get().split(" ")
            minimum = float(values[0])
            maximum = float(values[1])
            f = heat_map_maker(minimum, maximum, 1)
            chart_type = FigureCanvasTkAgg(f, bar_chart_window)
            chart_type.draw()
            chart_type.get_tk_widget().place(relx=0.51, rely=0.2)

        bar_chart_window = tk.Toplevel(root)
        bar_chart_window.geometry('1200x760')
        chart_type = FigureCanvasTkAgg(plt.gcf(), bar_chart_window)
        chart_type.draw()
        chart_type.get_tk_widget().place(relx=0.0, rely=0.2, relwidth=0.5)
        m = tk.Message(bar_chart_window, font=('Calibri', 15), highlightthickness=0, bd=0, width=1000, justify='center')
        m['text'] = "Enter the minimum value and the maximum value (exclusive) separated by a space.\nPress Enter to " \
                    "create the heatmap with these specifications."
        m.place(relx=0.25, rely=0.05)
        e = tk.Entry(bar_chart_window, font=('Calibri', 15))
        e.place(relx=0.44, rely=0.14)
        e.bind("<Return>", scope_heat_map)


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
def heat_map():
    import hyperspy.api as hs
    global single_values
    if file is None:
        label1['text'] = "Please load a file before creating a heat map.\n"
    elif distances is None:
        label1['text'] = "Please analyze the file before creating a heat map.\n"
    else:
        data = single_values.copy()
        df = DataFrame(data, columns=arange(len(data[0])), index=arange(len(data)))
        print(df)
        print(single_values)
        fig = px.imshow(df, color_continuous_scale=["blue", "green", "red"])
        fig.show()

        # creates diffraction pattern pop-up windows
        def image_gallery(event):
            global file
            values = e.get().split(" ")
            x0 = int(values[0])
            y0 = int(values[1])
            x1 = int(values[2])
            y1 = int(values[3])
            indexx = x1

            for x in range(x1 - x0 + 1):
                indexy = y1
                for y in range(y1 - y0 + 1):
                    s = PixelatedSTEM(hs.signals.Signal2D(file.inav[indexx, indexy]))
                    st = s.template_match_ring(r_inner=1, r_outer=6, lazy_result=True, show_progressbar=False)
                    peak_array = st.find_peaks(method='dog', min_sigma=0.8, max_sigma=15, sigma_ratio=1.9,
                                               threshold=0.42, overlap=0.5, lazy_result=False, show_progressbar=True)
                    s.add_peak_array_as_markers(peak_array)
                    # plt.plot(s)
                    s.plot()
                    ax = s._plot.signal_plot.ax
                    ax.set_xlabel("pixel(" + str(indexx) + "_" + str(indexy) + ")")
                    # plt.title("pixel(" + str(indexx) + "_" + str(indexy) + ")")
                    # plt.show()
                    indexy -= 1
                    y += 1
                indexx -= 1
                x += 1
            plt.show()

        # creates pop-up UI that calls image_gallery upon user input
        bar_chart_window = tk.Toplevel(root)
        bar_chart_window.geometry('500x400')
        m = tk.Message(bar_chart_window, font=('Calibri', 15), highlightthickness=0, bd=0, justify='left')
        m['text'] = "A new window should open displaying the heatmap created. If you would like to view specific " \
                    "diffraction patterns, enter the starting x and the y value and the ending x and y value " \
                    "separated by a space. Press Enter to display these diffraction patterns."
        m.place(relx=0.05, rely=0.05, relwidth=0.9, relheight=0.5)
        e = tk.Entry(bar_chart_window, bg='#F4F4F4', font=('Calibri', 15), justify='left', highlightthickness=0,
                     bd=0, fg='#373737', borderwidth=2, relief="groove")
        e.place(relx=0.1, rely=0.7, relwidth=0.8, relheight=0.1)
        e.bind("<Return>", image_gallery)

        bar_chart_window.mainloop()


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
    label = tk.Label(frame, text='Menu', bg='#FFFFFF', font=('Times New Roman', 40), fg='#373737')
    label.place(relx=0.40, rely=0.1, relwidth=0.2, relheight=0.1)

    # Text Output box
    label1 = tk.Message(frame, bg='#F3F3F3', font=('Calibri', 15), anchor='nw', justify='left', highlightthickness=0,
                        bd=0, width=1500, fg='#373737', borderwidth=2, relief="groove")
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
    button.place(relx=0.33, rely=0.22, relwidth=0.34, relheight=0.05)

    button1 = tk.Button(frame, text='Start Analysis', bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0, bd=0,
                        activebackground='#D4D4D4', activeforeground='#252525',
                        command=lambda: set_curr_func("analysis"), pady=0.02, fg='#373737', borderwidth='2',
                        relief="groove")
    button1.place(relx=0.33, rely=0.28, relwidth=0.34, relheight=0.05)

    button2 = tk.Button(frame, text='Create Bar Chart', bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0, bd=0,
                        activebackground='#D4D4D4', activeforeground='#252525',
                        command=lambda: bar_chart(), pady=0.02, fg='#373737', borderwidth='2',
                        relief="groove")
    button2.place(relx=0.33, rely=0.34, relwidth=0.34, relheight=0.05)

    button3 = tk.Button(frame, text='Create Heat Map', bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0, bd=0,
                        activebackground='#D4D4D4', activeforeground='#252525',
                        command=lambda: heat_map(), pady=0.02, fg='#373737', borderwidth='2',
                        relief="groove")
    button3.place(relx=0.33, rely=0.40, relwidth=0.34, relheight=0.05)

    button4 = tk.Button(frame, text='Transfer Data to .csv', bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0,
                        bd=0, activebackground='#D4D4D4', activeforeground='#252525',
                        command=lambda: set_curr_func("to_csv"), pady=0.02, fg='#373737', borderwidth='2',
                        relief="groove")
    button4.place(relx=0.33, rely=0.46, relwidth=0.34, relheight=0.05)

    root.mainloop()
    if path.exists("temp.png"):
        remove("temp.png")
