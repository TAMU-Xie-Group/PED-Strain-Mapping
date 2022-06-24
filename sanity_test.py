# Sanity test
# lets user upload 1 file
# shows gaussian blurred version of image
# shows peak points in diffraction spots on blurred image


import tkinter as tk
from tkinter import filedialog
from hyperspy.api import load
from PIL import Image, ImageTk
from pixstem.api import PixelatedSTEM
from numpy import array, mean, sqrt
from skimage.feature import match_template
from scipy.ndimage import gaussian_filter
from skimage.restoration import estimate_sigma, denoise_nl_means
from scipy.signal import wiener


# prompts file dialog for user to select file
def load_file():
    global file
    label1['text'] = label1['text'] + "Loading file...\n"
    input_file = filedialog.askopenfilename()
    root.update()
    try:
        file = load(input_file)
        label1['text'] = label1['text'] + "File loaded.\n"
    except ValueError:
        label1['text'] = label1['text'] + "Please select a file and try again.\n"
    except OSError:
        label1['text'] = label1['text'] + "Error loading. Please check the file path and try again.\n"


# applies filter to first image in .blo file
# calls find_peaks() function
# calls display() function
def filter_and_analyze():
    if file is not None:
        s = PixelatedSTEM(file.inav[0, 0])
        original = array(s)

        # # FILTERS
        # Gaussian is the default filter
        # Uncomment nlm (non-local means) or wien (wiener) to use those filters
        # # patch_size for 580 - 1, patch_size for 144 = 3
        # nlm = denoise_nl_means(original, h=1.15*sigma_est, fast_mode=True, patch_size=1, patch_distance=6, )
        # wien = wiener(original, 5, 3)

        sigma_est = mean(estimate_sigma(original, ))
        gaussian = gaussian_filter(original, 1.15 * sigma_est)

        filtered = array(gaussian)

        print(original)
        print(filtered)
        f = find_peaks(filtered)
        display(filtered, f)


# calculates and returns distance between two points
def distance(x1, y1, x2, y2):
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2))


# finds the centers of diffraction spots in filtered image
def find_peaks(filtered):
    template = filtered[265:320, 265:320]  # CHANGE BASED ON IMAGE SIZE
    print(template)
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
            if distance(point[0], point[1], tempList[j][0], tempList[j][1]) < 15:
                # minimum center distance can be changed
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
    print(peaks_list)
    return peaks_list


# displays filtered image and identified diffraction spots in pop-up window
def display(im, p_list):
    img = Image.fromarray(im)

    r = tk.Toplevel(root)

    c = tk.Canvas(r, height=720, width=1080)
    c.pack()
    f = tk.Frame(r, bg='#FFFFFF')
    f.place(relwidth=1, relheight=1)
    c2 = tk.Canvas(r, width=img.size[1], height=img.size[0])
    c2.place(relx=0.5, rely=0.5, anchor="center")
    image = ImageTk.PhotoImage(img)
    c2.create_image(0, 0, anchor='nw', image=image)

    for p in p_list:
        c2.create_oval(p[1]-2, p[0]-2, p[1]+2, p[0]+2, fill='#05FF00', outline='#05FF00')

    r.mainloop()


# main UI
if __name__ == "__main__":
    file = None

    HEIGHT = 600
    WIDTH = 800

    root = tk.Tk()

    canvas = tk.Canvas(root, height=HEIGHT, width=WIDTH)
    canvas.pack()
    frame = tk.Frame(root, bg='#FFFFFF')
    frame.place(relwidth=1, relheight=1)

    # Menu Label
    label = tk.Label(frame, text='Sanity Test', bg='#FFFFFF', font=('Calibri', 40), fg='#373737')
    label.place(relx=0.30, rely=0.05, relwidth=0.4, relheight=0.1)

    # Text Output box
    label1 = tk.Message(frame, bg='#F0F0F0', font=('Calibri', 15), anchor='nw', justify='left', highlightthickness=0,
                        bd=0, width=1500, fg='#373737')
    label1.place(relx=0.2, rely=0.6, relwidth=0.6, relheight=0.30)

    # Buttons
    button1 = tk.Button(frame, text='Load File', bg='#F0F0F0', font=('Calibri', 30), highlightthickness=0, bd=0,
                        activebackground='#D4D4D4', activeforeground='#252525', command=lambda: load_file(),
                        pady=0.02, fg='#373737')
    button1.place(relx=0.3, rely=0.25, relwidth=0.4, relheight=0.1)

    button2 = tk.Button(frame, text='Filter and Analyze', bg='#F0F0F0', font=('Calibri', 30), highlightthickness=0,
                        bd=0, activebackground='#D4D4D4', activeforeground='#252525',
                        command=lambda: filter_and_analyze(),
                        pady=0.02, fg='#373737')
    button2.place(relx=0.3, rely=0.4, relwidth=0.4, relheight=0.1)

    root.mainloop()
