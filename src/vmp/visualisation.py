# import packages
import argparse
import time
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from PIL import Image
from scipy.optimize import curve_fit



# SHOW IMAGE : --------------------------------------------------------------------------------------------------
def show_image():

    # function to plot an edge
    def plot_edge(edge):
        plt.plot([edge[0], edge[2]], [edge[1], edge[3]], color = 'grey', zorder=1, linewidth = 1)

    # load files for first snapshot
    nr = 0
    seeds = np.loadtxt('files/seed_list' + str(nr) +  '.csv', delimiter=',', skiprows=1)
    verticies = np.loadtxt('files/vertex_list' + str(nr) +  '.csv', delimiter=',', skiprows=1)
    edges = np.loadtxt('files/edge_list' + str(nr) +  '.csv', delimiter=',', skiprows=1)

    # optional style settings
    #plt.style.use('dark_background')
    #plt.figure(figsize=(10,10), facecolor='black')
    plt.figure(figsize=(6,6))



    # plotting all the edges
    print("plotting all the edges:")
    progress_bar = tqdm(total=len(edges), position=0, leave=True)

    for edge in edges:
        plot_edge(edge)
        progress_bar.update(1)

    progress_bar.close()

    # optional: scatter seeds and verticies
    if (len(edges) < 4000):
        plt.scatter(seeds[:, 0], seeds[:, 1], s=10,  zorder = 2)
        plt.scatter(verticies[:, 0], verticies[:, 1], s=10, zorder =3)

    # set x and y limits
    plt.xlim(0, 1)
    plt.ylim(0, 1)

    print("saving image ...")
    plt.savefig("../figures/single_picture.pdf")
    if (len(edges) < 60000):
        print("done! show image ...")
        plt.show()
    else:
        print("done saving image. you can finde the file under ../figures/single_picture.pdf")



# BENCHMARK : --------------------------------------------------------------------------------------------------
def benchmark():
    # load needed benchmark data
    data = np.loadtxt('benchmarks/time_benchmark.csv', delimiter=',', skiprows=1)

    # fitfunction
    def line(x, m, b):
        return m*x + b

    # fit for O(n^3) algorithm
    popt, pcov = curve_fit(line, np.log(data[-4:, 0]), np.log(data[-4:,1]/1000000))
    m = popt[0]
    b = popt[1]
    x_fit1= np.linspace(min(data[:, 0]), max(data[:, 0])*5, 100)
    y_fit1 = np.exp(line(np.log(x_fit1), m, b))

    plt.figure(figsize=(8,6))
    plt.title('performance benchmark')
    plt.xlabel('number of seeds')
    plt.ylabel('time in [s]')
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(data[:, 0], data[:, 1]/1000000, label = "benchmark", marker = '+', zorder = 2)
    plt.plot(x_fit1, y_fit1, label=f'fit: y = {np.exp(b):.2e} * x^{m:.2f}', color='grey', linestyle = '--', zorder = 1)
    plt.legend()
    plt.savefig("../figures/time_benchmark.png")
    plt.show()

    # load data
    memory = np.loadtxt('benchmarks/memory_benchmark.csv', delimiter=',', skiprows=1)

    # fit linear function 
    popt1, pcov1 = curve_fit(line, memory[:, 0], memory[:,1]/1024/1024)
    print((popt1[0]*1024), ' KB per Seed')
    m1 = popt1[0]
    b1 = popt1[1]


    # plot the linear fit
    x_fit = np.linspace(min(memory[1:, 0]), max(memory[:, 0])*1.2, 100)
    y_fit = line(x_fit, m1, b1)
    plt.title('memory benchmark')
    plt.xlabel('number of seeds')
    plt.ylabel('max_rss_memory in [MB]')
    plt.plot(memory[:, 0], memory[:, 1]/1024/1024, label = 'point insertion', marker = '+', color = "orange")
    plt.plot(x_fit, y_fit, label=f'fit: y = {m1*1024:.2f} KB * x + {b1*1024:.2f} KB', color='grey', linestyle = '--')
    plt.legend(loc = 'best')
    plt.savefig("../figures/memory_benchmark.png")
    plt.show()



# MOVING MESH ANIMATION : ----------------------------------------------------------------------------
def mm_anim(num_frames, frames_per_second):
    print("num frames:", num_frames, "fps", frames_per_second)

        # function to plot an edge
    def plot_edge(ax, edge):
        ax.plot([edge[0], edge[2]], [edge[1], edge[3]], color='grey', zorder=1)

    # function to update the frame
    def update(frame):

        # clear the frame
        plt.clf()
    
        # load data from the current file
        seeds = np.loadtxt(f'files/seed_list{frame}.csv', delimiter=',', skiprows=1)
        verticies = np.loadtxt(f'files/vertex_list{frame}.csv', delimiter=',', skiprows=1)
        edges = np.loadtxt(f'files/edge_list{frame}.csv', delimiter=',', skiprows=1)

        # plot edges
        for edge in edges:
            plot_edge(plt.gca(), edge)
    
        # optional : scatter seeds and verticies
        plt.scatter(seeds[:, 0], seeds[:, 1], s=25, zorder=2)
        #plt.scatter(verticies[:, 0], verticies[:, 1], s=10, zorder=3)

        # title and stuff
        plt.title(f'Frame {frame}')
        plt.axis('equal')  # Keep the aspect ratio equal for better visualization
        progress_bar.update(1)
        #print(frame)


    # Number of frames (assuming you have files numbered from 0 to num_frames)
    num_frames = num_frames - 1

    # Create the animation
    print('create animation...')
    animation = FuncAnimation(plt.figure(figsize=(10, 10)), update, frames=num_frames, interval=200, repeat=False)

    # Save the animation as a GIF using Pillow writer
    animation_file = '../figures/voronoi_animation.gif'
    progress_bar = tqdm(total=num_frames, position=0, leave=True)
    animation.save(animation_file, writer='pillow', fps=frames_per_second)
    progress_bar.close()
    print('done')



# GRID GENERATION ANIMATION : ----------------------------------------------------------------------------
def gg_anim(num_frames, frames_per_second):
    
    # function to plot an edge
    def plot_edge(ax, edge):
        ax.plot([edge[0], edge[2]], [edge[1], edge[3]], color='grey', zorder=1)

    # function to update the frame
    def update(frame):

        # clear the frame
        plt.clf()
        plt.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color = 'grey', zorder =1)
        plt.xlim(0,1)
        plt.ylim(0,1)

        # load data from the current file
        seeds = np.loadtxt(f'files/seed_list{frame+1}.csv', delimiter=',', skiprows=1)
        verticies = np.loadtxt(f'files/vertex_list{frame+1}.csv', delimiter=',', skiprows=1)
        edges = np.loadtxt(f'files/edge_list{frame+1}.csv', delimiter=',', skiprows=1)

        #seeds_end = np.loadtxt(f'build/files/seed_list98.csv', delimiter=',', skiprows=1)
        #plt.scatter(seeds_end[:, 0], seeds_end[:, 1], s=25, zorder = 2)

        # plot edges
        for edge in edges:
            plot_edge(plt.gca(), edge)
    
        # optional : scatter seeds and verticies
        if frame == 0:
            plt.scatter(seeds[0], seeds[1], s = 25, zorder = 2)
        else:
            plt.scatter(seeds[:, 0], seeds[:, 1], s=25, zorder=2)
        #plt.scatter(verticies[:, 0], verticies[:, 1], s=10, zorder=3)


        # title and stuff
        plt.title(f'generate grid animation, Frame {frame}')
        plt.axis('equal')  # Keep the aspect ratio equal for better visualization
        progress_bar.update(1)
        #print(frame)


    # Number of frames (assuming you have files numbered from 0 to num_frames)
    num_frames = num_frames - 2

    # Create the animation
    print('create animation...')
    animation = FuncAnimation(plt.figure(figsize=(10, 10)), update, frames=num_frames, interval=200, repeat=False)

    # Save the animation as a GIF using Pillow writer
    animation_file = '../figures/voronoi_animation.gif'
    progress_bar = tqdm(total=num_frames, position=0, leave=True)
    animation.save(animation_file, writer='pillow', fps=frames_per_second)
    progress_bar.close()
    print('done')



# MAIN FUNCTION : ---------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='Script to visualize outputs from vornoi mesh project')
    parser.add_argument('-program', type=int, help='which visualisation to run (0: show image, 1: benchmark, 2: moving mesh animation, 3: grid generation animation)')
    parser.add_argument('-num_frames', type=int, help='number of frames for the animations')
    parser.add_argument('-fps', type=int, help='fps for the animations')

    args = parser.parse_args()

    # Access the value of the option
    program = args.program
    num_frames = args.num_frames
    fps = args.fps

    # start the specified program
    print("starting python visualisation...")
    if (program == 0):
        show_image()
    elif (program == 1):
        benchmark()
    elif (program == 2):
        mm_anim(num_frames, fps)
    elif (program == 3):
        gg_anim(num_frames, fps)

if __name__ == '__main__':
    main()