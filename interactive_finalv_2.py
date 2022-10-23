import requests
import tkinter as tk
from tkinter import filedialog
from hyperspy.api import load
from PIL import Image, ImageTk, ImageEnhance
from pixstem.api import PixelatedSTEM
import numpy as np
import cv2
import time
import io
import tqdm
from multiprocessing import Pool
from skimage.feature import match_template
from scipy.ndimage import gaussian_filter
from skimage.restoration import estimate_sigma
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
from pandas import DataFrame
from statistics import mode

file = None
analysis_size = None
point1 = None
point2 = None
distance_arr = []
reciprocal_distance_arr = []
strain_arr = []


def point_selection(entry_event):
    global file, analysis_size, point1, point2, reciprocal_distance_arr, strain_arr

    reciprocal_distance_arr = []
    strain_arr = []
    image_len = len(file.data[0][0])

    analysis_size = [int(x) for x in user_entry.get().split(" ")]
    user_entry.delete(0, tk.END)

    def reset_coordinates():
        global point1, point2
        point1 = None
        point2 = None
        analysis_log['text'] = "Strain mapping: please click on the two points you would like to " \
                               "analyze from the diffraction pattern above.\n"

    def confirm_coordinates():
        global point1, point2
        if point1 is not None and point2 is not None:
            print("Selected points are ", point1, " and ", point2)
            analysis_log['text'] = analysis_log['text'] + "Starting analysis...\n"
            multiprocessing_func(point1, point2)
            c2.unbind('<Button-1>')
            r.destroy()
            label_output['text'] = "Analysis complete.\n"
        else:
            reset_coordinates()

    def get_mouse_coordinates(mouse_event):
        global point1, point2
        if point1 is None:
            point1 = (int(mouse_event.x * image_len / 400), int(mouse_event.y * image_len / 400))  # get mouse position
            analysis_log['text'] = analysis_log['text'] + "point1 = " + str(point1[0]) + " " + str(point1[1]) + "\n"
            print("point1 is ", point1)
        elif point2 is None:
            point2 = (int(mouse_event.x * image_len / 400), int(mouse_event.y * image_len / 400))  # get mouse position
            analysis_log['text'] = analysis_log['text'] + "point2 = " + str(point2[0]) + " " + str(point2[1]) + "\n"
            print("point2 is ", point2)
        r.update()

    r = tk.Toplevel(root)
    r.title('')

    c = tk.Canvas(r, height=640, width=840)
    c.pack()

    f = tk.Frame(r, bg='#FFFFFF')
    f.place(relwidth=1, relheight=1)

    analysis_log = tk.Message(f, bg='#FFFFFF', font=('Calibri', 15), anchor='nw', justify='left',
                              highlightthickness=0, bd=0, width=756)
    analysis_log.place(relx=0.05, rely=0.65, relwidth=0.9, relheight=0.2)

    reset_button = tk.Button(f, text='Reset',
                             bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0,
                             bd=0, activebackground='#D4D4D4', activeforeground='#252525',
                             command=lambda: reset_coordinates(),
                             pady=0.02, fg='#373737', borderwidth='2', relief="groove")
    reset_button.place(relx=0.15, rely=0.90, relwidth=0.30, relheight=0.07)

    confirm_button = tk.Button(f, text='Confirm',
                               bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0,
                               bd=0, activebackground='#D4D4D4', activeforeground='#252525',
                               command=lambda: confirm_coordinates(),
                               pady=0.02, fg='#373737', borderwidth='2', relief="groove")
    confirm_button.place(relx=0.55, rely=0.90, relwidth=0.30, relheight=0.07)

    c2 = tk.Canvas(r, width=400, height=400)
    c2.place(relx=0.5, anchor='n')

    preview_data = file.data[25, 25]
    preview_img = Image.fromarray(preview_data)
    preview_img = preview_img.resize((400, 400))
    preview_img = ImageTk.PhotoImage(preview_img)

    c2.create_image(0, 0, anchor='nw', image=preview_img)
    preview_filtered = filter_img([preview_data, 0, 0, 0, 0])
    p_list = find_peaks(preview_filtered)
    for p in p_list[0]:
        p = (p[0] * 400 / 580, p[1] * 400 / 580)  # resize for preview image
        c2.create_oval(p[1] - 4, p[0] - 4, p[1] + 4, p[0] + 4, fill='#ff0000', outline='#ff0000')
    c2.bind('<Button-1>', get_mouse_coordinates)
    analysis_log['text'] = "Strain mapping: please click on the two points you would like to " \
                           "analyze from the diffraction pattern above.\n"

    r.mainloop()


def multiprocessing_func(p1, p2):
    global file, analysis_size, distance_arr, reciprocal_distance_arr, strain_arr

    distance_arr, reciprocal_distance_arr, strain_arr = [], [], []
    img_list = []
    for i in range(analysis_size[0]):
        for j in range(analysis_size[1]):
            img_list.append([file.data[i][j], i, j, p1, p2])

    # Filtering
    filtered_img_list = []
    print("Filtering all diffraction patterns: ")
    start_time = time.time()
    pool = Pool(processes=None)
    for output in tqdm.tqdm(pool.imap_unordered(filter_img, img_list), total=len(img_list)):
        filtered_img_list.append(output)
        pass
    pool.close()
    del img_list
    print('Filtering time:', round(time.time() - start_time, 2))

    # Peak finding
    results_peaks = []
    print("Finding peaks for all diffraction patterns: ")
    start_time = time.time()
    pool = Pool(processes=None)
    for output in tqdm.tqdm(pool.imap_unordered(find_peaks, filtered_img_list), total=len(filtered_img_list)):
        results_peaks.append(output)
        pass
    pool.close()
    print('Peak-finding time:', round(time.time() - start_time, 2))

    # Distance calculation
    results_distance = []
    print("Matching points and calculating distance for all diffraction patterns: ")
    start_time = time.time()
    pool = Pool(processes=None)
    for output in tqdm.tqdm(pool.imap_unordered(calculate_points_distance, results_peaks),
                            total=len(filtered_img_list)):
        results_distance.append(output)
        pass
    pool.close()
    del results_peaks
    sorted_distance = sorted(results_distance, key=lambda x: (x[1] * 10 + x[2]))
    del results_distance

    # Creating distance and reciprocal distance list
    for data in sorted_distance:
        distance_arr.append(data[0])
        if data[0] != float('NaN'):
            reciprocal_distance_arr.append(1 / data[0])
        else:
            reciprocal_distance_arr.append(float('NaN'))
    del sorted_distance
    distance_arr = np.reshape(distance_arr, (-1, analysis_size[1]))
    reciprocal_distance_arr = np.reshape(reciprocal_distance_arr, (-1, analysis_size[1]))
    print("distance arr:", distance_arr)
    print('Calculation time:', round(time.time() - start_time, 2))


def filter_img(input_data):
    input_img, row, col, p1, p2 = input_data[0], input_data[1], input_data[2], input_data[3], input_data[4]


    sigma_est = np.mean(estimate_sigma(input_img, ))
    input_img = gaussian_filter(input_img, 1.15 * sigma_est)

    return [input_img, row, col, p1, p2]


# finds the centers of diffraction spots in filtered image
def find_peaks(input_data):
    input_img, row, col, p1, p2 = input_data[0], input_data[1], input_data[2], input_data[3], input_data[4]

    # kernel = np.ones((10, 10), np.uint8)
    # erosion = cv2.erode(input_img, kernel, iterations=1)
    # input_img = erosion

    # input_img = cv2.imread(input_img)
    # input_img = cv2.resize(input_img, (580, 580))
    template = input_img[260:315, 257:312]  # CHANGE BASED ON IMAGE SIZE
    # slice1 = int(len(input_img)/2.2)
    # slice2 = int(len(input_img)/1.8)
    # template = input_img[slice1:slice2, slice1:slice2]
    result = match_template(input_img, template, pad_input=True)
    # only takes points greater than the threshold r-value
    temp_list = []
    for i in range(len(result)):
        for j in range(len(result[i])):
            if result[i][j] > 0.87:  # correlation value
                temp_list.append((i, j))
    # removes duplicate spots that are too close to each other
    peaks_list = []
    while len(temp_list) > 0:
        j = 0
        temp = []
        point = temp_list[0]
        while j < len(temp_list):
            if calculate_distance(point[0], point[1], temp_list[j][0], temp_list[j][1]) < 25:
                # min center dist can be changed
                temp.append(temp_list[j])
                temp_list.pop(j)
            else:
                j = j + 1
        pnt = temp[0]
        for j in range(len(temp)):
            if result[pnt[0]][pnt[1]] < result[temp[j][0]][temp[j][1]]:
                pnt = temp[j]
        peaks_list.append(pnt)
    return [peaks_list, row, col, p1, p2]


# calculates and returns distance between two points
def calculate_distance(x1, y1, x2, y2):
    return np.sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2))


def calculate_points_distance(input_data):
    peaks, row, col, p1, p2 = input_data[0], input_data[1], input_data[2], input_data[3], input_data[4]

    peak1 = None
    peak2 = None
    for peak in peaks:
        if calculate_distance(peak[0], peak[1], p1[0], p1[1]) < 30:
            peak1 = peak
        if calculate_distance(peak[0], peak[1], p2[0], p2[1]) < 30:
            peak2 = peak

    distance = float('NaN')
    if peak1 is not None and peak2 is not None:
        distance = calculate_distance(peak1[0], peak1[1], peak2[0], peak2[1])

    return [distance, row, col]


def calculate_strain(event):
    global reciprocal_distance_arr, strain_arr, analysis_size

    # avg_reciprocal_dist = np.nanmean(reciprocal_distance_arr)
    user_selected_stat = float(user_entry.get())

    for r_dist in np.array(reciprocal_distance_arr).flatten():
        if r_dist == float('NaN'):
            strain_arr.append(float('NaN'))
        else:
            strain = (r_dist - user_selected_stat) / user_selected_stat
            strain_arr.append(strain)
    strain_arr = np.reshape(strain_arr, (-1, analysis_size[1]))
    print("strain arr: ", strain_arr)
    create_strain_map()


def create_strain_map():
    global file, strain_arr, distance_arr, reciprocal_distance_arr, analysis_size

    df = DataFrame(strain_arr, columns=np.arange(analysis_size[1]), index=np.arange(analysis_size[0]))
    print(df)

    if len(strain_arr) > 0:
        print("Mean strain: ", round(np.nanmean(strain_arr), 8))
        print("Strain std: ", round(np.nanstd(strain_arr), 8))

    fig = px.imshow(df, color_continuous_midpoint=0, zmin=-0.05, zmax=0.05, color_continuous_scale='turbo')
    buf = io.BytesIO()
    pio.write_image(fig, buf)
    heat_img = Image.open(buf)
    heat_img.show()
    fig.show()


########################################################################################################################
# BUTTON FUNCTIONS

# prompts file dialog for user to select file
def load_file():
    global file
    label_output['text'] = "Loading file...\n"
    input_file = filedialog.askopenfilename()
    root.update()
    try:
        file = PixelatedSTEM(load(input_file))
        label_output['text'] = label_output['text'] + "File loaded.\n"
    except ValueError:
        label_output['text'] = label_output['text'] + "Please select a file and try again.\n"
    except OSError:
        label_output['text'] = label_output['text'] + "Error loading. Please check the file path and try again.\n"


def analysis():
    global file
    if file is None:
        label_output['text'] = "Please load a file before analyzing.\n"
    else:
        label_output['text'] = "Please enter the number of rows and columns you would like to " \
                               "\nanalyze, as integers separated by spaces. Press Enter when ready.\n"
        user_entry.bind("<Return>", point_selection)


# creates bar chart pop-up UI
def bar_chart():
    global distance_arr
    if file is None:
        label_output['text'] = "Please load a file before creating a bar chart.\n"
    elif len(distance_arr) == 0:
        label_output['text'] = "Please analyze the file before creating a bar chart.\n"
    else:
        label_output['text'] = "Creating bar chart.\n"
        root.update()

        cleaned_distance_arr = []
        for dist in np.array(distance_arr).flatten():
            if dist != float('NaN'):
                cleaned_distance_arr.append(dist)

        fig, a = plt.subplots(figsize=(6, 5.5))
        plt.xlabel('Distance', fontsize=10)
        plt.ylabel('Counts', fontsize=10)
        plt.title('Distance Counts', fontsize=10)
        plt.hist(np.array(cleaned_distance_arr), bins=500)

        bar_chart_window = tk.Toplevel(root)
        bar_chart_window.geometry('600x600')
        chart_type = FigureCanvasTkAgg(plt.gcf(), bar_chart_window)
        chart_type.draw()
        chart_type.get_tk_widget().place(relx=0.0, rely=0.0, relwidth=1)


# opens a heat map as a tab in browser
# creates pop-up UI that can bring up user-requested diffraction patterns
def begin_strain_map():
    global strain_arr
    if file is None:
        label_output['text'] = "Please load a file before creating a heat map.\n"
    elif len(reciprocal_distance_arr) == 0:
        label_output['text'] = "Please analyze the file before creating a heat map.\n"
    else:
        label_output['text'] = 'Please enter a strain-free distance value for comparison. \n' \
                               f'For reference, here are some statistics from the dataset: \n' \
                               f'Mean reciprocal distance = {round(np.nanmean(reciprocal_distance_arr), 8)} \n' \
                               f'Median reciprocal distance = {round(np.nanmedian(reciprocal_distance_arr), 8)} \n' \
                               f'Mode reciprocal distance = ' \
                               f'{round(mode(reciprocal_distance_arr.flatten().tolist()), 8)} \n' \
                               f'STD reciprocal distance = {round(np.std(reciprocal_distance_arr), 8)} \n' \
                               'If the distribution is bimodal or multimodal, reference the histogram \n' \
                               'for the desired distance.'
        user_entry.bind("<Return>", calculate_strain)


# main UI
if __name__ == "__main__":
    HEIGHT = 700
    WIDTH = 800

    root = tk.Tk()
    root.title('')

    canvas = tk.Canvas(root, height=HEIGHT, width=WIDTH)
    canvas.pack()
    frame = tk.Frame(root, bg='#FFFFFF')
    frame.place(relwidth=1, relheight=1)

    # TAMU MSEN logo
    url = 'https://github.com/TAMU-Xie-Group/PED-Strain-Mapping/blob/main/msen.png?raw=true'
    image = Image.open(requests.get(url, stream=True).raw)
    image = image.resize((200, 40))
    img = ImageTk.PhotoImage(image)
    label_logo = tk.Label(frame, image=img, bg='#FFFFFF')
    label_logo.place(relx=0.05, rely=0.05, anchor='w')

    # Menu Label
    label_title = tk.Label(frame, text='PED Strain Mapping',
                           bg='#FFFFFF', font=('Times New Roman', 40), fg='#373737')
    label_title.place(relx=0.20, rely=0.1, relwidth=0.6, relheight=0.1)

    # Text Output box
    label_output = tk.Message(frame, bg='#F3F3F3', font=('Calibri', 15),
                              anchor='nw', justify='left', highlightthickness=0,
                              bd=0, width=1500, fg='#373737', borderwidth=2, relief="groove")
    label_output['text'] = "This program was designed by Ainiu Wang, using code developed by " \
                           "\nAniket Patel, Aaron Barbosa, and Marcus Hansen."
    label_output.place(relx=0.1, rely=0.54, relwidth=0.8, relheight=0.32)

    # Entry box
    user_entry = tk.Entry(frame, bg='#F3F3F3', font=('Calibri', 15), justify='left', highlightthickness=0,
                          bd=0, width=1500, fg='#373737', borderwidth=2, relief="groove")
    user_entry.place(relx=0.1, rely=0.88, relwidth=0.8, relheight=0.05)

    # Buttons
    button1 = tk.Button(frame, text='Load File',
                        bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0,
                        bd=0, activebackground='#D4D4D4', activeforeground='#252525',
                        command=lambda: load_file(),
                        pady=0.02, fg='#373737', borderwidth='2', relief="groove")
    button1.place(relx=0.29, rely=0.22, relwidth=0.42, relheight=0.05)

    button2 = tk.Button(frame, text='Start Analysis',
                        bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0,
                        bd=0, activebackground='#D4D4D4', activeforeground='#252525',
                        command=lambda: analysis(),
                        pady=0.02, fg='#373737', borderwidth='2', relief="groove")
    button2.place(relx=0.29, rely=0.28, relwidth=0.42, relheight=0.05)

    button3 = tk.Button(frame, text='Create Distance Bar Chart',
                        bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0,
                        bd=0, activebackground='#D4D4D4', activeforeground='#252525',
                        command=lambda: bar_chart(),
                        pady=0.02, fg='#373737', borderwidth='2', relief="groove")
    button3.place(relx=0.29, rely=0.34, relwidth=0.42, relheight=0.05)

    button4 = tk.Button(frame, text='Create Strain Heat Map',
                        bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0,
                        bd=0, activebackground='#D4D4D4', activeforeground='#252525',
                        command=lambda: begin_strain_map(),
                        pady=0.02, fg='#373737', borderwidth='2', relief="groove")
    button4.place(relx=0.29, rely=0.40, relwidth=0.42, relheight=0.05)

    # button5 = tk.Button(frame, text='Export Strain Data to .csv',
    #                     bg='#F3F3F3', font=('Calibri', 20), highlightthickness=0,
    #                     bd=0, activebackground='#D4D4D4', activeforeground='#252525',
    #                     command=lambda: to_csv(input_filename),
    #                     pady=0.02, fg='#373737', borderwidth='2', relief="groove")
    # button5.place(relx=0.29, rely=0.46, relwidth=0.42, relheight=0.05)

    root.mainloop()