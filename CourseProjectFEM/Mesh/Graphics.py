import numpy as np
import matplotlib.pyplot as plt

def main():
    
    x = []
    y = []

    figure, axes = plt.subplots()
    plt.xlabel("x")
    plt.ylabel("y")

    # plt.xlim(-3,3)
    # plt.ylim(-3,3)
    
    #plt.xticks(np.arange(-4, 4, 0.2))
    #plt.yticks(np.arange(-4, 4, 0.2))
    
    with open("AllRibs.dat") as file:
        for line in file:
            x1C, y1C, x2C, y2C = line.split()
            x.append(float(x1C))
            y.append(float(y1C))
            x.append(float(x2C))
            y.append(float(y2C))
            axes.plot(x, y, '-s', markersize=2, color='black')
            x.clear()
            y.clear()

    #plt.grid()
    plt.show()

if __name__ == "__main__":
    main()
